"""
MPI-based parallel polymerisation runner that works directly with
classifier output.

.. important::

    This runner **must** be launched via ``mpirun`` (or ``mpiexec``).
    It requires at least **2 MPI ranks** (1 master + 1 worker) and
    will raise ``RuntimeError`` if run as a plain Python script.

    .. code-block:: bash

        # Minimum: 2 ranks (1 master, 1 worker)
        mpirun -np 2 python -m polyscript.polymerizer.mpi_df

        # Production: 4 ranks (1 master, 3 workers)
        mpirun -np 4 python my_runner_script.py

Accepts a single DataFrame — the result of
:meth:`PolyScriptClassifier.classify` — and automatically extracts
monomer subsets for each polymer class / monomer-type rule.

Usage::

    from polyscript.classifier import PolyScriptClassifier
    from polyscript.polymerizer.mpi_df import MpiRunnerDf

    classifier = PolyScriptClassifier()
    classified = classifier.classify(my_df, col_name="smiles")

    runner = MpiRunnerDf(exp_name="run42", df=classified, smiles_col="smiles")
    runner.run()

    # Restrict to specific monomer types
    runner = MpiRunnerDf(
        exp_name="run42", df=classified, smiles_col="smiles",
        monomer_classes=["vinyl", "epo", "lactone"],
    )
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
from polyscript.polymerizer.mpi import MpiPolymerizer, SKIP_CLS

# ---------------------------------------------------------------------------
# Package root for default paths
# ---------------------------------------------------------------------------
_PKG_DIR = Path(__file__).resolve().parent


class MpiRunnerDf(MpiPolymerizer):
    """Distributed polymerisation runner fed by a classifier DataFrame.

    .. important::

        Must be launched via ``mpirun`` with at least **2 ranks**::

            mpirun -np 4 python my_script.py

        Running as a plain ``python`` script will raise ``RuntimeError``.

    Extends ``MpiPolymerizer`` so that instead of requiring pre-split
    parquet files on disk, it accepts a single DataFrame — the output of
    ``PolyScriptClassifier.classify()`` — and filters it on-the-fly for
    each monomer-type rule.

    Column format expected in *df*:

        ============  =====  ======  ========
        smiles        vinyl  epo     lactone
        ============  =====  ======  ========
        C=C           True   False   False
        C1OC1         False  True    False
        O=C1OCCC1     False  False   True
        ============  =====  ======  ========

    The boolean columns correspond to monomer functional-group types
    and are the same names returned by the classifier.
    """

    def __init__(
        self,
        exp_name: str,
        df: pd.DataFrame,
        smiles_col: str = "psmiles",
        *,
        monomer_classes: Optional[List[str]] = None,
        rule_dir: Optional[str] = None,
        output_path: Optional[str] = None,
        ckpt_path: str = "ckpt.db",
        log_filepath: str = "polymerizer.log",
        col_name: str = "psmiles",
        df_chunk_size: int = 1000,
        task_chunk_size: Optional[int] = None,
        log_interval: int = 100,
        use_rich: bool = True,
    ):
        """Initialise the runner from a classifier DataFrame.

        Extracts monomer subsets from boolean columns in *df*, creates
        chunked parquet files on disk, and populates the checkpoint
        database.

        Args:
            exp_name (str): Experiment name for the output subdirectory.
            df (pandas.DataFrame): Classifier output with a SMILES column
                plus one boolean column per monomer type.
            smiles_col (str): Name of the column containing SMILES strings.
            monomer_classes (list of str or None): Optional whitelist of
                monomer type columns to include.
            rule_dir (str or None): Path to the rule directory.
            output_path (str or None): Base output directory.
            ckpt_path (str): SQLite checkpoint database filename.
            log_filepath (str): Log file path.
            col_name (str): Column name to use in generated parquet chunks.
            df_chunk_size (int): Number of rows per chunk file.
            task_chunk_size (int or None): MPI task dispatch chunk size.
            log_interval (int): Progress log interval in completed tasks.
            use_rich (bool): Enable Rich progress bars on master.
        """
        # ── Default paths ────────────────────────────────────────────
        if rule_dir is None:
            rule_dir = os.path.join(_PKG_DIR, "rules")
        if output_path is None:
            output_path = os.path.join(_PKG_DIR, "outputs")

        self._df_chunk_size = df_chunk_size
        self._smiles_col = smiles_col

        # ── Build monomer-type → DataFrame mapping ───────────────────
        # Each value is a single-column DataFrame with the given col_name
        self._monomer_dfs: Dict[str, pd.DataFrame] = {}
        for col in df.columns:
            if col == smiles_col:
                continue
            if df[col].dtype != bool:
                continue
            if monomer_classes is not None and col not in monomer_classes:
                continue
            subset = df.loc[df[col], [smiles_col]].copy()
            subset = subset.rename(columns={smiles_col: col_name})
            if len(subset) > 0:
                self._monomer_dfs[col] = subset

        # ── Delegate MPI init + logger + checkpoint to base ──────────
        # We pass a dummy input_data_split_dir because the base requires
        # it, but we override _populate_checkpoint and _build_config so
        # it's never actually used.
        super().__init__(
            exp_name=exp_name,
            input_data_split_dir=_PKG_DIR,  # unused by this subclass
            output_path=output_path,
            ckpt_path=ckpt_path,
            rule_dir=rule_dir,
            log_filepath=log_filepath,
            col_name=col_name,
            chunk_size=task_chunk_size,
            log_interval=log_interval,
            use_rich=use_rich,
        )

        # ── Post-init: create chunk files + populate checkpoint ──────
        if self._is_master:
            self.logger.info(
                "Loaded %d monomer types from classifier DataFrame "
                "(smiles_col=%r, %d total rows)",
                len(self._monomer_dfs), smiles_col, len(df),
            )
            if monomer_classes:
                self.logger.info(
                    "Monomer class filter active: %s", monomer_classes
                )
            self._chunks_dir = os.path.join(self.exp_path, "_chunks")
            os.makedirs(self._chunks_dir, exist_ok=True)
            self._chunk_map: Dict[str, Dict[int, str]] = {}
            self._create_chunk_files()

            # Replace the base's file-system populate with chunk-based
            self._populate_checkpoint_df()

    # ──────────────────────────────────────────────────────────────────
    # Chunk file management
    # ──────────────────────────────────────────────────────────────────

    def _create_chunk_files(self):
        """Split each monomer-type DataFrame into parquet chunk files on disk.

        Each monomer type gets its own set of numbered chunk files stored
        under ``exp_path/_chunks/``.  The chunk-to-path mapping is saved
        in ``self._chunk_map``.
        """
        for mon_type, df in self._monomer_dfs.items():
            n_rows = len(df)
            n_chunks = max(1, (n_rows + self._df_chunk_size - 1) // self._df_chunk_size)
            self._chunk_map[mon_type] = {}

            for i in range(n_chunks):
                start = i * self._df_chunk_size
                end = min(start + self._df_chunk_size, n_rows)
                chunk_df = df.iloc[start:end].copy()
                chunk_path = os.path.join(
                    self._chunks_dir, f"{mon_type}_chunk_{i:05d}.parquet"
                )
                chunk_df.to_parquet(chunk_path, index=False)
                self._chunk_map[mon_type][i] = chunk_path

            self.logger.info(
                "Monomer '%s': %d rows -> %d chunks (%d rows/chunk)",
                mon_type, n_rows, n_chunks, self._df_chunk_size,
            )

    def _get_chunk_count(self, mon_type: str) -> int:
        """Return the number of chunk files for a monomer type.

        Args:
            mon_type (str): Monomer type name.

        Returns:
            int: Number of chunk files (0 if the type is unknown).
        """
        return len(self._chunk_map.get(mon_type, {}))

    def _get_chunk_path(self, mon_type: str, chunk_idx: int) -> str:
        """Return the file path for a specific chunk.

        Args:
            mon_type (str): Monomer type name.
            chunk_idx (int): Zero-based chunk index.

        Returns:
            str: Path to the chunk parquet file.
        """
        return self._chunk_map[mon_type][chunk_idx]

    # ──────────────────────────────────────────────────────────────────
    # Override: populate checkpoint from DataFrame chunks
    # ──────────────────────────────────────────────────────────────────

    def _populate_checkpoint_df(self):
        """Populate the checkpoint DB from DataFrame-derived chunks.

        Overrides the base-class file-system-based population.  Loads
        rules from ``ps_gen.pkl`` and enumerates all (class, mon_type1,
        mon_type2, chunk_pair) tasks using the chunk count metadata.
        """
        # Load rules the same way the base class does
        import joblib
        ps_cls = joblib.load(os.path.join(self.rule_dir, "ps_gen.pkl"))

        rows: List[Tuple] = []
        for p_class, rules_list in ps_cls.items():
            if p_class in SKIP_CLS:
                continue
            for values in rules_list:
                mon_type1 = values[0]
                mon_type2 = values[1]
                n1 = self._get_chunk_count(mon_type1)
                if n1 == 0:
                    self.logger.warning(
                        "No chunks for '%s' — skipping %s/%s_%s",
                        mon_type1, p_class, mon_type1, mon_type2,
                    )
                    continue

                mon_pair_key = f"{mon_type1}_{mon_type2}"
                if mon_type2 == "none":
                    n_chunks_total = n1
                    permuted = [str(i) for i in range(n1)]
                else:
                    n2 = self._get_chunk_count(mon_type2)
                    if n2 == 0:
                        self.logger.warning(
                            "No chunks for '%s' — skipping %s/%s",
                            mon_type2, p_class, mon_pair_key,
                        )
                        continue
                    permuted = [f"{i}+{j}" for i in range(n1) for j in range(n2)]
                    n_chunks_total = n1 * n2

                self.logger.info(
                    "  %s/%s: %d tasks (mon1_chunks=%d%s)",
                    p_class, mon_pair_key, n_chunks_total, n1,
                    f", mon2_chunks={n2}" if mon_type2 != "none" else "",
                )

                for batch_key in permuted:
                    rows.append(
                        (p_class, mon_pair_key, mon_type1, mon_type2, batch_key, 0)
                    )

        if not rows:
            self.logger.warning("No tasks generated — check rules and DataFrame.")
            return

        self.checkpoint.populate(rows)
        self.logger.info("Checkpoint DB populated with %d tasks", len(rows))

    # ──────────────────────────────────────────────────────────────────
    # Override: build config from chunk indices
    # ──────────────────────────────────────────────────────────────────

    def _build_config(
        self, p_class: str, mon_type1: str, mon_type2: str, batch_key: str
    ) -> Optional[Dict]:
        """Build a configuration dict from chunk indices.

        Overrides the base-class method.  *batch_key* is a chunk index
        (or ``"idx1+idx2"`` for pairs), resolved via ``_get_chunk_path``.

        Args:
            p_class (str): Polymer class label.
            mon_type1 (str): Monomer type 1.
            mon_type2 (str): Monomer type 2.
            batch_key (str): Chunk index or pair expression.

        Returns:
            dict or None: Configuration dict, or ``None`` if resolution fails.
        """
        output_dir = os.path.join(self.exp_path, p_class)
        os.makedirs(output_dir, exist_ok=True)
        output_filename = (
            f"{mon_type1}_{mon_type2}_{batch_key.replace('+', '_x_')}.parquet"
        )
        output_path = os.path.join(output_dir, output_filename)

        if "+" in batch_key:
            idx1_str, idx2_str = batch_key.split("+")
            mon1_path = self._get_chunk_path(mon_type1, int(idx1_str))
            mon2_path = self._get_chunk_path(mon_type2, int(idx2_str))
        else:
            mon1_path = self._get_chunk_path(mon_type1, int(batch_key))
            mon2_path = None

        return {
            "mon1_path": mon1_path,
            "mon2_path": mon2_path,
            "output_path": output_path,
            "col_name": self.col_name,
            "P_class": p_class,
            "mon_type1": mon_type1,
            "mon_type2": mon_type2 if mon_type2 != "none" else None,
        }


# ---------------------------------------------------------------------------
# Convenience function
# ---------------------------------------------------------------------------


def run_from_dataframes(
    exp_name: str,
    df: pd.DataFrame,
    smiles_col: str = "psmiles",
    **kwargs,
) -> Dict:
    """One-shot: create ``MpiRunnerDf`` from a classifier DataFrame and run it.

    .. important::

        Must be called from a script launched via ``mpirun``::

            mpirun -np 4 python my_script.py

    .. code-block:: python

        from polyscript.classifier import PolyScriptClassifier
        from polyscript.polymerizer.mpi_df import run_from_dataframes

        classified = PolyScriptClassifier().classify(my_df, col_name="smiles")
        run_from_dataframes("exp01", df=classified, smiles_col="smiles")
    """
    runner = MpiRunnerDf(exp_name=exp_name, df=df, smiles_col=smiles_col, **kwargs)
    return runner.run()
