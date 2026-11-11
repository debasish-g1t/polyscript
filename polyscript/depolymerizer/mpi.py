"""
MPI-parallel depolymerizer.

``MpiDepolymerizer`` is a distributed batch-depolymerization runner that
composes from ``polyscript.utils.mpi`` components (MpiLogger, MpiProtocol,
WorkSplitter, ProgressTracker) to scatter sequences across MPI ranks,
run ``PolyScriptDepolymerizer``, and gather per-class statistics.

Usage (launch with mpirun)::

    from polyscript.depolymerizer.mpi import MpiDepolymerizer

    runner = MpiDepolymerizer(
        input_path="benchmark.parquet",
        mode="validate",
    )
    stats = runner.run()
"""

from __future__ import annotations

import time
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd

from polyscript.depolymerizer import PolyScriptDepolymerizer
from polyscript.utils.mpi import (
    MpiLogger,
    MpiProtocol,
    ProgressTracker,
    WorkSplitter,
)


class MpiDepolymerizer:
    """Distributed depolymerization runner.

    Components:
    * ``self.logger``   — MpiLogger (direct-stderr + file)
    * ``self.mpi``      — MpiProtocol (scatter/gather)
    * ``self.splitter`` — WorkSplitter (round-robin)
    * ``self.progress`` — ProgressTracker (elapsed, rate, ETA)

    Parameters
    ----------
    input_path : str
        Path to a parquet or CSV file containing polymer sequences.
    mode : str
        ``"validate"`` (PolyScript sequences) or ``"convert"`` (raw PSMILES).
    smiles_col : str, optional
        Column name for SMILES strings (auto-detected when ``None``).
    output_dir : str
        Directory for results and logs.
    log_filepath : str
        Log filename written inside *output_dir*.
    """

    def __init__(
        self,
        input_path: str,
        *,
        mode: str = "validate",
        smiles_col: Optional[str] = None,
        output_dir: str = ".",
        log_filepath: str = "depolymerizer.log",
        log_interval: int = 100,
    ):
        # ── MPI init ──────────────────────────────────────────────────
        from mpi4py import MPI as _MPI  # noqa: N813

        comm = _MPI.COMM_WORLD
        self._rank = comm.Get_rank()
        self._size = comm.Get_size()
        self._is_master = self._rank == 0

        if self._size < 2:
            raise RuntimeError(
                "MpiDepolymerizer needs at least 2 MPI ranks. "
                f"Got {self._size}. Launch with mpirun."
            )

        # ── Components ────────────────────────────────────────────────
        self.mpi = MpiProtocol(comm, self._rank, self._size, self._is_master)
        self.splitter = WorkSplitter()
        self._log_interval = max(log_interval, 1)

        self._input_path = str(Path(input_path).resolve())
        self._mode = mode
        self._smiles_col = smiles_col
        self._output_dir = Path(output_dir).resolve()

        if self._is_master:
            self.logger = MpiLogger(
                name="depolymerizer",
                log_filepath=str(self._output_dir / log_filepath),
                progress_file=str(self._output_dir / "progress.txt"),
                use_rich=False,
                use_direct_stderr=True,
            )
        else:
            import logging
            self.logger = logging.getLogger(__name__ + f".rank{self._rank}")
            self.logger.setLevel(logging.WARNING)

        self.mpi.barrier()

    # ------------------------------------------------------------------
    # Lifecycle
    # ------------------------------------------------------------------

    def close(self):
        """Flush and close the logger if present."""
        if hasattr(self, "logger"):
            try:
                self.logger.flush()
                self.logger.close()
            except AttributeError:
                pass

    def __enter__(self):
        """Support context-manager usage (``with MpiDepolymerizer(...) as``)."""
        return self

    def __exit__(self, *args):
        """Clean up on context-manager exit."""
        self.close()

    # ------------------------------------------------------------------
    # Master: read input & scatter
    # ------------------------------------------------------------------

    def _load_sequences(self) -> tuple:
        """Read the input file and return (sequences, indices, classes)."""
        ip = Path(self._input_path)
        if ip.suffix == ".parquet":
            df = pd.read_parquet(ip)
        else:
            df = pd.read_csv(ip)

        if self._mode == "convert":
            col = self._smiles_col
            if col is None:
                for candidate in ["smiles", "canonical_smiles", "psmiles", "polym"]:
                    if candidate in df.columns:
                        col = candidate
                        break
                if col is None:
                    raise KeyError("No SMILES column found. Pass --smiles_col.")
            sequences = df[col].dropna().tolist()
            classes = ["unknown"] * len(sequences)
        else:
            sequences = df["polym"].tolist()
            classes = (
                df["polymer_class"].tolist()
                if "polymer_class" in df.columns
                else ["unknown"] * len(sequences)
            )

        indices = list(range(len(sequences)))
        return sequences, indices, classes

    # ------------------------------------------------------------------
    # Worker: process one sequence
    # ------------------------------------------------------------------

    @staticmethod
    def _process_sequence(
        dp: PolyScriptDepolymerizer,
        seq: str,
        mode: str,
    ) -> list:
        """Run depolymerization on a single sequence. Returns pathway list."""
        if mode == "convert":
            return dp.convert_psmiles(seq)
        return dp.validate(seq)

    # ------------------------------------------------------------------
    # run
    # ------------------------------------------------------------------

    def run(self) -> Dict:
        """Run the distributed depolymerization job.

        On the master rank (0), loads sequences, scatters them to
        workers, processes its own chunk, gathers results, and returns
        summary statistics.  Worker ranks process their chunk and return
        an empty dict.

        Returns:
            dict: Summary stats with keys ``total``, ``success``,
            ``failed``, ``pathways`` (master only; workers return ``{}``).
        """
        if self._is_master:
            return self._run_master()
        else:
            self._run_worker()
            return {}

    # ------------------------------------------------------------------
    # Master orchestration
    # ------------------------------------------------------------------

    def _run_master(self) -> Dict:
        """Master orchestration: load, scatter, process, gather, summarize.

        Loads the input file, distributes sequences across MPI ranks
        via round-robin, processes the master's own chunk, gathers
        results from all ranks, logs per-class statistics, and saves
        results/failures to parquet files.

        Returns:
            dict: Summary stats (``total``, ``success``, ``failed``,
            ``pathways``).
        """
        comm = self.mpi.comm
        t_start = time.monotonic()

        # ── Load & scatter ────────────────────────────────────────────
        sequences, indices, classes = self._load_sequences()
        total = len(sequences)

        seq_chunks = self.splitter.round_robin(sequences, self._size)
        idx_chunks = self.splitter.round_robin(indices, self._size)
        cls_chunks = self.splitter.round_robin(classes, self._size)

        self.logger.info(
            "Loaded %d sequences, distributing across %d ranks", total, self._size,
        )

        comm.scatter(seq_chunks, root=0)
        comm.scatter(idx_chunks, root=0)
        comm.scatter(cls_chunks, root=0)

        # ── Process master's own chunk ────────────────────────────────
        my_seqs = seq_chunks[0]
        my_idxs = idx_chunks[0]
        my_clss = cls_chunks[0]

        results, fails, ok, bad, per_class = self._process_chunk(
            my_seqs, my_idxs, my_clss, t_start,
        )

        # ── Gather from workers ───────────────────────────────────────
        all_results = comm.gather(results, root=0)
        all_fails = comm.gather(fails, root=0)
        all_per_class = comm.gather(per_class, root=0)

        # ── Summarise ─────────────────────────────────────────────────
        aggregated: dict[str, int] = {}
        for pc_dict in all_per_class:
            for pclass, count in pc_dict.items():
                aggregated[pclass] = aggregated.get(pclass, 0) + count

        # Include master's own results
        flat_results = [r for chunk in all_results for r in chunk]
        flat_fails = [f for chunk in all_fails for f in chunk]

        seq_ok = {r.seq_index for r in flat_results}
        seq_fail_idx = {f["seq_index"] for f in flat_fails}
        n_total = len(seq_ok) + len(seq_fail_idx)

        self.logger.info("=" * 60)
        self.logger.info("DEPOLYMERIZATION SUMMARY")
        self.logger.info("=" * 60)
        self.logger.info("  MPI ranks:       %d", self._size)
        self.logger.info("  Total sequences: %d", n_total)
        self.logger.info("  Success:         %d", len(seq_ok))
        self.logger.info("  Failed:          %d", len(seq_fail_idx))
        self.logger.info("  Total pathways:  %d", len(flat_results))

        for pclass in sorted(aggregated):
            pclass_fails = sum(
                1 for f in flat_fails if f["polymer_class"] == pclass
            )
            total_class = aggregated[pclass] + pclass_fails
            pct = 100 * aggregated[pclass] / max(1, total_class)
            self.logger.info(
                "    %-20s %4d/%d  (%.0f%%)",
                pclass, aggregated[pclass], total_class, pct,
            )

        # ── Save results ──────────────────────────────────────────────
        if flat_results:
            out_df = pd.DataFrame([
                {
                    "seq_index": r.seq_index,
                    "input_monomers": r.input_monomers,
                    "output_monomers": r.monomers,
                    "exact_match": r.is_exact_match,
                    "reaction": r.reaction,
                    "polymer": r.polymer,
                    "polymer_pattern": r.polymer_pattern,
                }
                for r in flat_results
            ])
            results_path = self._output_dir / "depolymerize_results.parquet"
            out_df.to_parquet(results_path)
            self.logger.info("  Saved %d pathways to %s", len(out_df), results_path)

        if flat_fails:
            fail_df = pd.DataFrame(flat_fails)
            fail_path = self._output_dir / "depolymerize_fails.parquet"
            fail_df.to_parquet(fail_path)
            self.logger.info("  Saved %d fail cases to %s", len(fail_df), fail_path)

        self.mpi.shutdown_workers()
        self.logger.flush()
        return {
            "total": n_total,
            "success": len(seq_ok),
            "failed": len(seq_fail_idx),
            "pathways": len(flat_results),
        }

    # ------------------------------------------------------------------
    # Shared: process a chunk of sequences
    # ------------------------------------------------------------------

    def _process_chunk(
        self,
        sequences: list,
        indices: list,
        classes: list,
        t_start: float,
    ) -> tuple:
        """Run depolymerization on a local chunk. Returns results tuple."""
        dp = PolyScriptDepolymerizer()
        results = []
        fails = []
        ok = 0
        bad = 0
        per_class: dict[str, int] = {}
        n = len(sequences)

        progress = ProgressTracker(total=n, log_interval=self._log_interval)

        for i, (seq, gidx, pclass) in enumerate(zip(sequences, indices, classes)):
            seq_results = self._process_sequence(dp, seq, self._mode)

            if seq_results:
                ok += 1
                per_class[pclass] = per_class.get(pclass, 0) + 1
                for r in seq_results:
                    r.seq_index = gidx
                    results.append(r)
            else:
                bad += 1
                fails.append(
                    {"seq_index": gidx, "sequence": seq, "polymer_class": pclass}
                )

            progress.tick(ok=bool(seq_results))

            if progress.should_log() or i == n - 1:
                self.logger.info(
                    "rank %2d: %s",
                    self._rank, progress.format_line(),
                )
                progress.mark_logged()

        self.logger.info(
            "rank %2d: DONE  ok=%d fail=%d out of %d  (%.1fs)",
            self._rank, ok, bad, n, progress.elapsed,
        )
        return results, fails, ok, bad, per_class

    # ------------------------------------------------------------------
    # Worker loop
    # ------------------------------------------------------------------

    def _run_worker(self):
        """Worker loop: receive chunk from master, process, send back results."""
        comm = self.mpi.comm

        # Receive scattered data from master
        my_sequences = comm.scatter(None, root=0)
        my_indices = comm.scatter(None, root=0)
        my_classes = comm.scatter(None, root=0)

        t_start = time.monotonic()
        results, fails, ok, bad, per_class = self._process_chunk(
            my_sequences, my_indices, my_classes, t_start,
        )

        # Send results back
        comm.gather(results, root=0)
        comm.gather(fails, root=0)
        comm.gather(per_class, root=0)
