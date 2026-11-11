"""Core polymerizer implementation: the ``Polymerizer`` class and helpers."""

import json
import os
import pickle
import sys
from functools import lru_cache
from itertools import islice, product
from pathlib import Path
from typing import Optional

# Ensure the package directory is on sys.path so that sibling modules
# (smip_utils, etc.) can be imported from any working directory.
_PKG_DIR = str(Path(__file__).resolve().parent)
if _PKG_DIR not in sys.path:
    sys.path.insert(0, _PKG_DIR)


def get_rules_dir() -> Optional[str]:
    """Return the path to the ``rules/`` directory shipped with this package.

    Returns ``None`` if the directory does not exist (e.g. in a partial
    installation).
    """
    rules_dir = os.path.join(_PKG_DIR, "rules")
    return rules_dir if os.path.isdir(rules_dir) else None


import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from canonicalize_psmiles.canonicalize import canonicalize
from rdkit.Chem import AllChem
from polyscript.polymerizer.utils import *

# ---------------------------------------------------------------------------
# Worker-local parquet cache (3.2)
# ---------------------------------------------------------------------------


@lru_cache(maxsize=8)
def _cached_read_parquet(path: str):
    """Read a parquet file, caching the result in this process."""
    return pd.read_parquet(path)


# ---------------------------------------------------------------------------
# Module-level cache — the 8 rule files are read-only, so load them once
# and reuse across all Polymerizer instances (and across all MPI tasks).
# Keyed by the resolved ``db_file`` directory path.
# ---------------------------------------------------------------------------
_shared_data_cache: dict = {}


def _load_data_dict(db_file: str) -> dict:
    """Load and cache the polymerisation rule files from *db_file*.

    The returned dict is immutable in practice — callers must not mutate it.
    Subsequent calls with the same *db_file* hit the cache and return
    instantly (zero disk I/O).
    """
    cached = _shared_data_cache.get(db_file)
    if cached is not None:
        return cached

    with open(os.path.join(db_file, "mon_vals.json"), "r") as f:
        mon_vals = json.load(f)
    with open(os.path.join(db_file, "mon_dic.json"), "r") as f:
        mon_dic = json.load(f)
    with open(os.path.join(db_file, "mon_dic_inv.json"), "r") as f:
        mon_dic_inv = json.load(f)
    with open(os.path.join(db_file, "mon_lst.json"), "r") as f:
        monL = json.load(f)
    with open(os.path.join(db_file, "excl_lst.json"), "r") as f:
        exclL = json.load(f)
    with open(os.path.join(db_file, "ps_rxn.pkl"), "rb") as f:
        Ps_rxnL = pickle.load(f)
    with open(os.path.join(db_file, "ps_class.json"), "r") as f:
        Ps_classL = json.load(f)
    with open(os.path.join(db_file, "ps_gen.pkl"), "rb") as f:
        Ps_GenL = pickle.load(f)

    monLg = {int(k): v for k, v in monL.items()}
    exclLg = {int(k): v for k, v in exclL.items()}
    mon_dic_inv = {int(k): v for k, v in mon_dic_inv.items()}

    data_dict = {
        "monLg": monLg,
        "exclLg": exclLg,
        "mon_dic_inv": mon_dic_inv,
        "mon_dic": mon_dic,
        "mon_vals": mon_vals,
        "Ps_rxnL": Ps_rxnL,
        "Ps_classL": Ps_classL,
        "Ps_GenL": Ps_GenL,
    }

    _shared_data_cache[db_file] = data_dict
    return data_dict


# ---------------------------------------------------------------------------
# Canonicalization helpers for dedup (post-1.8)
# ---------------------------------------------------------------------------


def _extract_psmiles(polym_seq: str) -> str:
    """Extract the raw pSMILES from a PolyScript sequence.

    PolyScript format:
        {mon1}+{mon2}=>rxn_smarts>>atom_mapping=>psmiles
    Multi-step (||| separated) — take the pSMILES from the final segment.
    """
    # Take the last segment (after final |||) and extract after final =>
    last_segment = polym_seq.rsplit("|||", 1)[-1]
    if "=>" in last_segment:
        return last_segment.rsplit("=>", 1)[-1]
    return last_segment


@lru_cache(maxsize=65536)
def _canonicalize_safe(ps: str) -> str:
    """Canonicalize a pSMILES string; fall back to the input on failure.

    LRU-cached: within a batch, the same pSMILES extracted from different
    PolyScript sequences will hit the cache instead of re-constructing
    RDKit molecular graphs.
    """
    try:
        return canonicalize(ps)
    except Exception:
        return ps


# ---------------------------------------------------------------------------
# Thread-safe block processor (1.7)
# ---------------------------------------------------------------------------


def _process_bipolymer_block(
    block, mol1_map, mol2_map, reaction, monL, Ps_rxnL, P_class, seen=None
):
    """Process one block of (m1, m2) pairs.

    All arguments are read-only → safe to call concurrently from threads.
    *seen* is an optional set of ``(reactset_tuple, polym_smiles)`` keys;
    when provided, duplicates are filtered out during generation (1.8).
    Returns three parallel lists ready for DataFrame construction.
    """
    m1_col = []
    m2_col = []
    polym_col = []
    for m1, m2 in block:
        polym_list = seq_bipolymA(
            [mol1_map[m1], mol2_map[m2]],
            targ_rxn=reaction,
            monL=monL,
            Ps_rxnL=Ps_rxnL,
            P_class=P_class,
        )
        for polym in polym_list:
            if seen is not None:
                reactset = (m1, m2) if m1 < m2 else (m2, m1)
                # Dedup on canonicalized pSMILES — two different PolyScript
                # sequences can yield the same polymer backbone.
                canon_ps = _canonicalize_safe(_extract_psmiles(polym))
                key = (reactset, canon_ps)
                if key in seen:
                    continue
                seen.add(key)
            m1_col.append(m1)
            m2_col.append(m2)
            polym_col.append(polym)
    return m1_col, m2_col, polym_col


class Polymerizer:
    """Generate polymers from monomer pairs using pre-defined reaction rules.

    Loads polymerisation rule data (monomer dictionaries, reaction templates,
    exclusion lists, etc.) from a directory of rule files and provides methods
    to perform bi-polymerisation and homopolymerisation reactions on batches
    of monomer candidates.

    The rule data is cached at the module level, so multiple instances share
    the same in-memory structures.

    Attributes:
        db_file (str): Resolved path to the rule directory.
        data_dict (dict): Cached rule data with keys such as ``monLg``,
            ``exclLg``, ``mon_dic``, ``mon_vals``, ``Ps_rxnL``, ``Ps_classL``,
            and ``Ps_GenL``.
    """

    def __init__(self, db_file_path="rules"):
        """Initialize the polymerizer by loading rule data.

        Args:
            db_file_path (str): Path to the directory containing rule files
                (``mon_vals.json``, ``mon_dic.json``, ``ps_rxn.pkl``, etc.).
                Relative paths are resolved against the package directory.
        """
        db_file = os.path.join(str(Path(__file__).resolve().parent), db_file_path)
        self.db_file = db_file
        self.data_dict = self.load_data()
        # O(1) lookup: (P_class, mon_type1, mon_type2) → rdkit Reaction.
        # Built once at init from Ps_GenL; stored separately so
        # data_dict.values() still unpacks into exactly 8 items.
        # Uses indexing (rule[0], rule[1], rule[2]) instead of tuple
        # unpacking because some rule entries may carry extra metadata.
        self._reaction_index = {
            (str(pc), rule[0], rule[1]): rule[2]
            for pc, rules in self.data_dict["Ps_GenL"].items()
            for rule in rules
        }

    def load_data(self):
        """Return the cached rule data dict (loads from disk on first call only)."""
        return _load_data_dict(self.db_file)

    def get_reaction_match(self, P_class, mon_type1, mon_type2):
        """O(1) dict lookup — was a linear scan over Ps_GenL[P_class] before."""
        return self._reaction_index.get((str(P_class), mon_type1, mon_type2))

    def bipolymerize(
        self,
        mon_df1,
        mon_type1,
        P_class,
        mon_df2=None,
        mon_type2=None,
        candidate_col_name="smiles",
        output_path: Optional[str] = None,
    ):
        """Generate polymers from one or two monomer DataFrames.

        When *mon_df2* is provided, performs binary (copolymer) reactions
        between all pairs of monomers from the two DataFrames.  When
        *mon_df2* is ``None``, performs homopolymerisation on *mon_df1*.

        For large batches (estimated total pairs >= 500,000), results are
        streamed directly to a Parquet file via ``pyarrow`` — the function
        returns ``None`` and never materialises the full result DataFrame
        in memory.  Smaller batches are returned as a pandas DataFrame.

        Deduplication is performed during generation (by canonical pSMILES)
        rather than as a post-processing step.

        Args:
            mon_df1 (pandas.DataFrame): First monomer DataFrame.  Must
                contain a column named by *candidate_col_name* with SMILES
                strings.
            mon_type1 (str): Monomer functional-group type for *mon_df1*
                (e.g. ``"vinyl"``, ``"epo"``).
            P_class (str): Polymer class label (e.g. ``"polyolefin"``).
            mon_df2 (pandas.DataFrame or None): Optional second monomer
                DataFrame for copolymerisation.
            mon_type2 (str or None): Monomer type for *mon_df2*, or
                ``None`` for homopolymerisation.
            candidate_col_name (str): Column name containing SMILES
                strings in the DataFrames.  Defaults to ``"smiles"``.
            output_path (str or None): If provided and the total pair
                count exceeds the streaming threshold, results are written
                directly to this Parquet path instead of being returned.

        Returns:
            pandas.DataFrame or None: A DataFrame with columns ``mon1``,
            ``mon2``, ``polym``, ``polymer_class``, or ``None`` when
            streaming output was used.
        """
        monLg, exclLg, mon_dic_inv, mon_dic, mon_vals, Ps_rxnL, Ps_classL, Ps_GenL = (
            self.data_dict.values()
        )

        monL = {
            k: v
            for k, v in monLg.items()
            if k in mon_vals[0] + mon_vals[1] + mon_vals[2]
        }
        exclL = {
            k: v
            for k, v in exclLg.items()
            if k in mon_vals[0] + mon_vals[1] + mon_vals[2]
        }

        if "ROMol" in mon_df1.columns:
            DF = mon_df1.drop("ROMol", axis=1).dropna(subset=[candidate_col_name])
        else:
            DF = mon_df1.dropna(subset=[candidate_col_name])
        if mon_df2 is not None and mon_type2 is not None:
            if "ROMol" in mon_df2.columns:
                DF2 = mon_df2.drop("ROMol", axis=1).dropna(subset=[candidate_col_name])
            else:
                DF2 = mon_df2.dropna(subset=[candidate_col_name])
        else:
            DF2 = None

        DF_Pgen = pd.DataFrame(
            columns=["mon1", "mon2", "polym", "polymer_class"]
        )

        # ── parquet streaming (1.9) ────────────────────────────────────
        # Streaming via pyarrow avoids materialising the full result
        # DataFrame in memory, but adds per-block overhead.  For small
        # batches that fit comfortably in RAM, skip streaming and use
        # the faster single-shot pandas path.
        STREAM_THRESHOLD = 500_000  # pairs — below this, don't stream
        writer = None
        _output_path = None
        _writer_schema = None
        if output_path is not None:
            total_pairs_est = len(DF) * (len(DF2) if DF2 is not None else 1)
            if total_pairs_est >= STREAM_THRESHOLD:
                os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
                _output_path = output_path
                _writer_schema = pa.schema(
                    [
                        ("mon1", pa.string()),
                        ("mon2", pa.string()),
                        ("polym", pa.string()),
                        ("polymer_class", pa.string()),
                    ]
                )

        def _ensure_writer():
            """Return the ParquetWriter, creating it lazily on first call.

            Returns:
                pyarrow.parquet.ParquetWriter or None: The writer instance,
                or ``None`` if streaming output is not enabled.
            """
            nonlocal writer
            if writer is None and _output_path is not None:
                writer = pq.ParquetWriter(
                    _output_path,
                    _writer_schema,
                    compression="zstd",
                    compression_level=3,
                )
            return writer

        reaction = self.get_reaction_match(P_class, mon_type1, mon_type2)
        if reaction is None:
            return DF_Pgen
        Ps_rxnL_key = []
        Ps_rxnL_key = [
            k
            for k, v in Ps_rxnL.items()
            if AllChem.ReactionToSmarts(v) == AllChem.ReactionToSmarts(reaction)
        ]
        temp1 = list(DF[candidate_col_name])
        temp2 = []
        if len(temp1) != 0:
            if mon_type2 is not None and mon_df2 is not None:
                temp2 = list(DF2[candidate_col_name])
                if len(temp2) != 0:
                    total_pairs = len(temp1) * len(temp2)
                    p_class_str = str(P_class)
                    ps_rxn_key = int(Ps_rxnL_key[0])

                    # ── pre-convert SMILES → Mol (1.5) ────────────────────
                    mol1_map = {s: genmol(s) for s in set(temp1)}
                    mol2_map = {s: genmol(s) for s in set(temp2)}
                    seen = set()  # (1.8) dedup during generation

                    if total_pairs <= STREAM_THRESHOLD:
                        # ── fast path: small batch, process all at once ──
                        blocks = [
                            (m1, m2) for m1, m2 in product(temp1, temp2) if m1 != m2
                        ]
                        m1_col, m2_col, polym_col = _process_bipolymer_block(
                            blocks,
                            mol1_map,
                            mol2_map,
                            reaction,
                            monL,
                            Ps_rxnL,
                            P_class,
                            seen,
                        )
                        if m1_col:
                            DF_temp = pd.DataFrame(
                                {
                                    "mon1": m1_col,
                                    "mon2": m2_col,
                                    "polym": polym_col,
                                    "polymer_class": p_class_str,
                                    # "Ps_rxnL": ps_rxn_key,
                                }
                            )
                            DF_Pgen = pd.concat([DF_Pgen, DF_temp], ignore_index=True)
                    else:
                        # ── large batch: chunked lazy processing (1.3) ──
                        BLOCK_SIZE = max(5000, min(100_000, total_pairs // 50))
                        pair_iter = (
                            (m1, m2) for m1, m2 in product(temp1, temp2) if m1 != m2
                        )
                        while True:
                            block = list(islice(pair_iter, BLOCK_SIZE))
                            if not block:
                                break
                            m1_col, m2_col, polym_col = _process_bipolymer_block(
                                block,
                                mol1_map,
                                mol2_map,
                                reaction,
                                monL,
                                Ps_rxnL,
                                P_class,
                                seen,
                            )
                            if m1_col:
                                w = _ensure_writer()
                                if w is not None:
                                    # (1.9) stream directly to Parquet —
                                    # zero DataFrame materialisation
                                    table = pa.table(
                                        {
                                            "mon1": pa.array(m1_col, type=pa.string()),
                                            "mon2": pa.array(m2_col, type=pa.string()),
                                            "polym": pa.array(
                                                polym_col, type=pa.string()
                                            ),
                                            "polymer_class": pa.array(
                                                [p_class_str] * len(m1_col),
                                                type=pa.string(),
                                            ),
                                            # "Ps_rxnL": pa.array(
                                            #     [ps_rxn_key] * len(m1_col),
                                            #     type=pa.int64(),
                                            # ),
                                        }
                                    )
                                    w.write_table(table)
                                else:
                                    DF_temp = pd.DataFrame(
                                        {
                                            "mon1": m1_col,
                                            "mon2": m2_col,
                                            "polym": polym_col,
                                            "polymer_class": p_class_str,
                                            # "Ps_rxnL": ps_rxn_key,
                                        }
                                    )
                                    DF_Pgen = pd.concat(
                                        [DF_Pgen, DF_temp], ignore_index=True
                                    )
            else:
                temp2 = ["" for i in range(len(temp1))]
                mons = monL[mon_dic[mon_type1]]
                excls = exclL[mon_dic[mon_type1]]
                # ── pre-convert SMILES → Mol (1.5) ────────────────────────
                mol_map = {s: genmol(s) for s in set(temp1)}

                # ── direct loop (1.6): replaces two chained .apply(axis=1).
                # Columns are built as flat lists (no dict-per-row overhead).
                m1_col = []
                polym_col = []
                p_class_str = str(P_class)
                seen = set()  # (1.8) dedup during generation
                for m1 in temp1:
                    polym_list = seq_homopolymA(
                        mol_map[m1],
                        mons=mons,
                        excls=excls,
                        targ_mon1=mon_type1,
                        Ps_rxnL=Ps_rxnL,
                        mon_dic=mon_dic,
                        monL=monL,
                    )
                    for polym in polym_list:
                        # dedup key for homopolymers: ("", mon1) is the
                        # canonical reactset ("" < any SMILES string).
                        # Use canonicalized pSMILES to catch duplicates
                        # that differ only in star position.
                        canon_ps = _canonicalize_safe(_extract_psmiles(polym))
                        key = (("", m1), canon_ps)
                        if key in seen:
                            continue
                        seen.add(key)
                        m1_col.append(m1)
                        polym_col.append(polym)

                if m1_col:
                    w = _ensure_writer()
                    if w is not None:
                        # (1.9) stream directly to Parquet
                        table = pa.table(
                            {
                                "mon1": pa.array(m1_col, type=pa.string()),
                                "mon2": pa.array([""] * len(m1_col), type=pa.string()),
                                "polym": pa.array(polym_col, type=pa.string()),
                                "polymer_class": pa.array(
                                    [p_class_str] * len(m1_col), type=pa.string()
                                ),
                                # "Ps_rxnL": pa.array(
                                #     [ps_rxn_key] * len(m1_col), type=pa.int64()
                                # ),
                            }
                        )
                        w.write_table(table)
                    else:
                        DF_temp = pd.DataFrame(
                            {
                                "mon1": m1_col,
                                "mon2": "",
                                "polym": polym_col,
                                "polymer_class": p_class_str,
                                # "Ps_rxnL": ps_rxn_key,
                            }
                        )
                        DF_Pgen = pd.concat([DF_Pgen, DF_temp], ignore_index=True)

        # ── final cleanup ──────────────────────────────────────────────
        if writer is not None:
            writer.close()
            return None

        # Deduplication is now done during generation (1.8) — the
        # expensive explode/sort/drop_duplicates chain below is replaced
        # by a simple dropna + reset_index.
        DF_gendP = DF_Pgen.dropna(subset=["polym"])
        DF_gendP = DF_gendP.reset_index(drop=True)
        return DF_gendP

    def post_process(
        self,
        res_df,
        reduce_column=None,
        drop_outro=True,
    ):
        """Filter and clean polymer results.

        Uses Polars when available; falls back to pandas otherwise.
        Accepts either a pandas or Polars DataFrame.

        Args:
            res_df (pandas.DataFrame or polars.DataFrame): The raw
                polymerisation result DataFrame.
            reduce_column (str or None): Optional column name to drop
                from the result.
            drop_outro (bool): If ``True`` (default), strips any ``|||``
                -separated suffix from the ``polym`` column, keeping only
                the main polymerisation result.

        Returns:
            polars.DataFrame or pandas.DataFrame: The filtered result.
        """
        if isinstance(res_df, pd.DataFrame) and res_df.shape[0] == 0:
            return res_df

        try:
            import polars as pl

            _USE_POLARS = True
        except ImportError:
            _USE_POLARS = False

        if _USE_POLARS:
            df = pl.from_pandas(res_df) if isinstance(res_df, pd.DataFrame) else res_df

            if drop_outro:
                # Keep only the part before the first "|||"
                df = df.with_columns(pl.col("polym").str.split("|||").list.first())

            if reduce_column is not None:
                df = df.drop(reduce_column)

            return df

        # ── pandas fallback ────────────────────────────────────────────
        if drop_outro:
            res_df["polym"] = res_df["polym"].apply(lambda x: x.split("|||")[0])

        if reduce_column is not None:
            res_df = res_df.drop(columns=reduce_column)

        res_df = res_df.reset_index(drop=True)
        return res_df

    def write_parquet(self, res_df, path):
        """Write a polymer result DataFrame to a Parquet file.

        Args:
            res_df (pandas.DataFrame or polars.DataFrame): Result
                DataFrame to write.
            path (str): Output file path (``.parquet``).
        """
        res_df.to_parquet(path, index=False, compression="zstd")

    def execute(self, config):
        """Run a single polymerisation batch from a configuration dict.

        Reads monomer parquet files, runs ``bipolymerize``, and writes
        the result to the output path specified in *config*.

        Args:
            config (dict): A dictionary with the following keys:

                - ``mon1_path`` (str): Path to the first monomer parquet file.
                - ``mon2_path`` (str or None): Path to the optional second
                  monomer parquet file.
                - ``output_path`` (str): Path for the output parquet file.
                - ``col_name`` (str): Column name for SMILES strings.
                - ``P_class`` (str): Polymer class label.
                - ``mon_type1`` (str): Monomer type for the first file.
                - ``mon_type2`` (str): Monomer type for the second file.

        Returns:
            tuple: ``(success, error_msg)`` where *success* is a bool and
            *error_msg* is ``None`` on success or an error string on failure.
        """
        try:
            mon_df1 = _cached_read_parquet(config["mon1_path"])
            if config["mon2_path"] is not None:
                mon_df2 = _cached_read_parquet(config["mon2_path"])
            else:
                mon_df2 = None
            P_class = config["P_class"]
            mon_type1 = config["mon_type1"]
            mon_type2 = config["mon_type2"]
            if mon_type2 is None:
                mon_type2 = (
                    "none"  # bug fix as the reaction mapper does not handle None
                )

            res_df = self.bipolymerize(
                mon_df1,
                mon_type1,
                P_class,
                mon_df2,
                mon_type2,
                candidate_col_name=config["col_name"],
                output_path=config.get("output_path"),  # (1.9)
            )

            if res_df is None:
                # (1.9) bipolymerize already streamed to output_path
                return True, None

            if res_df.shape[0] == 0:
                ...
            else:
                self.write_parquet(res_df, config["output_path"])
            return True, None
        except Exception as e:
            return False, str(e)
