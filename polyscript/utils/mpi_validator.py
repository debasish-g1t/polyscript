"""
MPI-parallel sequence validator — ``MpiSeqValidator``.

Uses ``MpiLogger`` and ``MpiProtocol`` from ``polyscript.utils.mpi``.

Usage (launch with mpirun)::

    from polyscript.utils.mpi_validator import MpiSeqValidator

    validator = MpiSeqValidator(
        input_dir="./polymerizer_outputs",
        output_dir="./validated",
        rule_dir="./rules",
    )
    stats = validator.run()
"""

from __future__ import annotations

import json
import os
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import joblib
import pandas as pd
from canonicalize_psmiles.canonicalize import canonicalize
from rdkit import RDLogger
from polyscript.utils.parsers import NestedPolyScriptParser
from polyscript.utils.validators import SeqValidator

from polyscript.utils.mpi import MpiLogger, MpiProtocol

RDLogger.logger().setLevel(RDLogger.ERROR)

TAG_WORK = 0
TAG_RESULT = 1
TAG_STOP = 2

SKIP_CLS = ["COC", "rec:radi", "rec:cati", "rec:ani", "rec:coord"]


# ---------------------------------------------------------------------------
# Worker-side validation function  (top-level so MPI can call it everywhere)
# ---------------------------------------------------------------------------


def _validate_file(
    args: Tuple[str, str, str, str],
) -> Dict:
    """Validate every row in a single parquet file.

    Parameters
    ----------
    args : tuple
        (input_parquet_path, output_parquet_path, p_class, mon_pair_key)

    Returns
    -------
    dict
        ``p_class``, ``mon_pair_key``, ``batch_key``, ``total_rows``,
        ``valid``, ``invalid``, ``unique_canon``, ``nested_counts``,
        ``error_type_counts``, ``elapsed_s``.
    """
    import time

    t0 = time.monotonic()

    input_path, output_path, p_class, mon_pair_key = args

    parser = NestedPolyScriptParser()
    validator = SeqValidator()

    df = pd.read_parquet(input_path)
    total_rows = len(df)

    valid_count = 0
    invalid_count = 0

    seen_canonical: Dict[Tuple[str, str], int] = {}
    keep_indices: List[int] = []
    linear_indices: List[int] = []
    invalid_rows: List[Tuple[int, List[str]]] = []
    nested_counts: Dict[int, int] = defaultdict(int)
    error_counts: Dict[str, int] = defaultdict(int)

    batch_key = Path(input_path).stem

    for idx, row in df.iterrows():
        polym = str(row.get("polym", ""))

        parser.parse(polym)

        if not parser.sequences:
            invalid_count += 1
            error_counts["<E-parse-|no-sequences|>"] += 1
            invalid_rows.append((idx, ["<E-parse-|no-sequences|>"]))
            continue

        num_nested = len(parser.sequences)
        nested_counts[num_nested] += 1

        all_sub_valid = True
        sub_errors: List[str] = []
        for monomers, reaction_smarts, polymer_smiles in parser.sequences:
            if monomers is None or reaction_smarts is None or polymer_smiles is None:
                all_sub_valid = False
                sub_errors.append("<E-parse-|sub-sequence-parse-failed|>")
                continue
            errs = validator._validate(monomers, reaction_smarts, polymer_smiles)
            if errs:
                all_sub_valid = False
                sub_errors.extend(errs)

        if not all_sub_valid:
            invalid_count += 1
            for e in sub_errors:
                error_counts[e] += 1
            invalid_rows.append((idx, sub_errors))
            continue

        first_monomers = parser.sequences[0][0]
        first_mon = first_monomers[0] if first_monomers else ""
        last_polymer = parser.sequences[-1][2]

        try:
            canon_mon = (
                canonicalize(first_mon)
                if first_mon and first_mon != "none"
                else first_mon
            )
        except Exception:
            canon_mon = first_mon

        try:
            canon_poly = canonicalize(last_polymer) if last_polymer else last_polymer
        except Exception:
            canon_poly = last_polymer

        canon_key = (canon_mon, canon_poly)

        if canon_key not in seen_canonical:
            seen_canonical[canon_key] = idx
            keep_indices.append(idx)
            valid_count += 1
            if num_nested == 1:
                linear_indices.append(idx)

    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    if keep_indices:
        validated_df = df.iloc[keep_indices].reset_index(drop=True)
        validated_df.to_parquet(output_path, index=False, compression="zstd")
    else:
        df.iloc[:0].to_parquet(output_path, index=False, compression="zstd")

    invalid_csv: Optional[str] = None
    if invalid_rows:
        tmp_dir = os.path.join(os.path.dirname(output_path), "_tmp_invalid")
        os.makedirs(tmp_dir, exist_ok=True)
        invalid_csv = os.path.join(tmp_dir, f"{batch_key}_invalid.csv")
        invalid_indices, error_lists = zip(*invalid_rows)
        invalid_df = df.iloc[list(invalid_indices)].copy()
        invalid_df["_validation_errors"] = [" | ".join(errs) for errs in error_lists]
        invalid_df.to_csv(invalid_csv, index=False)

    linear_parquet: Optional[str] = None
    if linear_indices:
        tmp_dir = os.path.join(os.path.dirname(output_path), "_tmp_linear")
        os.makedirs(tmp_dir, exist_ok=True)
        linear_parquet = os.path.join(tmp_dir, f"{batch_key}_linear.parquet")
        linear_df = df.iloc[linear_indices].copy()
        linear_df.to_parquet(linear_parquet, index=False, compression="zstd")

    elapsed = time.monotonic() - t0

    return {
        "p_class": p_class,
        "mon_pair_key": mon_pair_key,
        "batch_key": batch_key,
        "total_rows": total_rows,
        "valid": valid_count,
        "invalid": invalid_count,
        "unique_canon": len(seen_canonical),
        "linear_count": len(linear_indices),
        "nested_counts": dict(nested_counts),
        "error_counts": dict(error_counts),
        "invalid_csv": invalid_csv,
        "linear_parquet": linear_parquet,
        "elapsed_s": round(elapsed, 3),
    }


# ---------------------------------------------------------------------------
# MpiSeqValidator
# ---------------------------------------------------------------------------


class MpiSeqValidator:
    """Distributed PolyScript sequence validator.

    Uses ``MpiLogger`` and ``MpiProtocol`` for MPI orchestration.
    Validates sequences via ``SeqValidator``, canonicalizes and
    deduplicates, and writes statistics.

    Parameters
    ----------
    input_dir : str
        Directory containing parquet files from polymerizer output.
    output_dir : str
        Directory for validated outputs and logs.
    rule_dir : str
        Path to rules directory (monomer-type definitions).
    ckpt_path : str, optional
        Path to SQLite checkpoint DB (auto-detected when ``None``).
    log_filepath : str
        Log filename written inside *output_dir*.
    log_interval : int
        Minimum work units between progress log lines.
    use_rich : bool
        Enable Rich console output (default ``True``).
    """

    def __init__(
        self,
        input_dir: str,
        output_dir: str,
        rule_dir: str,
        ckpt_path: Optional[str] = None,
        log_filepath: str = "validator.log",
        log_interval: int = 10,
        use_rich: bool = True,
    ):
        """Initialise the distributed validator.

        On the master rank this sets up logging, loads rules, and
        discovers input tasks.  Workers only initialise the MPI
        communicator.  All ranks synchronise with a barrier before
        returning.

        Args:
            input_dir: Directory containing input parquet files.
            output_dir: Directory for validated outputs and logs.
            rule_dir: Path to the rules directory (``ps_gen.pkl``).
            ckpt_path: Path to an optional SQLite checkpoint DB.
                Auto-detected from *input_dir* when ``None``.
            log_filepath: Log filename written inside *output_dir*.
            log_interval: Minimum work units between progress log lines.
            use_rich: Enable Rich console progress display.

        Raises:
            RuntimeError: If fewer than 2 MPI ranks are available.
            ValueError: If *input_dir* or *rule_dir* does not exist.
        """
        from mpi4py import MPI as _MPI  # noqa: N813

        comm = _MPI.COMM_WORLD
        self._comm = comm
        self._rank = comm.Get_rank()
        self._size = comm.Get_size()
        self._is_master = self._rank == 0

        if self._size < 2:
            raise RuntimeError("MpiSeqValidator needs at least 2 MPI ranks.")
        self._n_workers = self._size - 1

        self.mpi = MpiProtocol(comm, self._rank, self._size, self._is_master)

        if self._is_master:
            self.input_dir = input_dir
            if not os.path.isdir(self.input_dir):
                raise ValueError(f"Input directory not found: {self.input_dir}")
            self.output_dir = output_dir
            os.makedirs(self.output_dir, exist_ok=True)
            if not os.path.isdir(rule_dir):
                raise ValueError(f"Rule directory not found: {rule_dir}")
            self.rule_dir = rule_dir

            self._log_interval = max(log_interval, 1)
            self._use_rich = use_rich

            self.logger = MpiLogger(
                name="validator",
                log_filepath=os.path.join(self.output_dir, log_filepath),
                progress_file=os.path.join(self.output_dir, "progress.txt"),
                use_rich=use_rich,
                use_direct_stderr=True,
            )
            self._console = self.logger._console if use_rich else None

            self.logger.info("MPI validator: 1 master + %d workers", self._n_workers)

            self._ckpt_db: Any = None
            self._db_hit_count = 0
            self._db_miss_count = 0
            db_path = ckpt_path or os.path.join(self.input_dir, "ckpt.db")
            if os.path.exists(db_path):
                import sqlite3
                self._ckpt_db = sqlite3.connect(db_path)
                self.logger.info("Checkpoint DB: %s", db_path)
            else:
                self.logger.warning("Checkpoint DB not found at %s", db_path)

            self._load_rules()
            self._task_list = self._discover_tasks()
            self.logger.info("Discovered %d parquet files", len(self._task_list))

        self.mpi.barrier()

    # ------------------------------------------------------------------
    # Lifecycle
    # ------------------------------------------------------------------

    def close(self):
        """Release resources: close the checkpoint DB and flush/close the logger."""
        if self._is_master and self._ckpt_db is not None:
            self._ckpt_db.close()
            self._ckpt_db = None
        if hasattr(self, "logger"):
            try:
                self.logger.flush()
                self.logger.close()
            except AttributeError:
                pass

    def __enter__(self):
        return self

    def __exit__(self, *a):
        self.close()

    # ------------------------------------------------------------------
    # Rule loading
    # ------------------------------------------------------------------

    def _load_rules(self):
        """Load polymerisation rules from ``ps_gen.pkl``.

        Populates ``self._class_mon_pairs`` (polymer class → list of
        monomer-type pairs) and ``self._mon_pair_set`` (class → first
        pair string).  Polymer classes listed in ``SKIP_CLS`` are
        excluded.
        """
        ps_gen_path = os.path.join(self.rule_dir, "ps_gen.pkl")
        if not os.path.exists(ps_gen_path):
            self.logger.warning("Rule file %s not found.", ps_gen_path)
            self._class_mon_pairs: Dict[str, List[Tuple[str, str]]] = {}
            self._mon_pair_set: Dict[str, str] = {}
            return
        ps_cls = joblib.load(ps_gen_path)
        self._class_mon_pairs = {}
        self._mon_pair_set = {}
        for p_class, rules in ps_cls.items():
            if p_class in SKIP_CLS:
                continue
            pairs = []
            for values in rules:
                mt1 = values[0]
                mt2 = values[1] if len(values) > 1 else "none"
                pairs.append((mt1, mt2))
            self._class_mon_pairs[p_class] = pairs
            if pairs:
                self._mon_pair_set[p_class] = f"{pairs[0][0]}_{pairs[0][1]}"

    # ------------------------------------------------------------------
    # Task discovery
    # ------------------------------------------------------------------

    def _discover_tasks(self) -> List[Tuple[str, str, str, str]]:
        """Walk *input_dir* and build the task list.

        Each task is a 4-tuple of
        ``(input_path, output_path, p_class, mon_pair_key)``.

        Returns:
            A sorted list of task tuples for all discovered parquet files.
        """
        tasks: List[Tuple[str, str, str, str]] = []
        input_root = Path(self.input_dir)
        output_root = Path(self.output_dir)
        _mpk_cache: Dict[Tuple[str, str], str] = {}
        for p_class_dir in sorted(input_root.iterdir()):
            if not p_class_dir.is_dir() or p_class_dir.name.startswith("_"):
                continue
            p_class = p_class_dir.name
            for pf in sorted(p_class_dir.glob("*.parquet")):
                input_path = str(pf)
                rel = pf.relative_to(input_root)
                output_path = str(output_root / rel)
                mpk = self._resolve_mon_pair_key(p_class, pf.stem, _mpk_cache)
                tasks.append((input_path, output_path, p_class, mpk))
        from collections import Counter
        mpk_counts = Counter(t[3] for t in tasks)
        self.logger.info("Monomer pairs: %d unique across %d files",
                         len(mpk_counts), len(tasks))
        if self._ckpt_db is not None:
            self.logger.info("DB resolution: %d hits, %d fallbacks",
                             self._db_hit_count, self._db_miss_count)
        return tasks

    def _resolve_mon_pair_key(self, p_class, filename_stem, cache):
        """Resolve the monomer-pair key for a given file.

        First tries the checkpoint DB, then falls back to matching
        monomer-type pairs loaded from the rules file.

        Args:
            p_class: Polymer class name.
            filename_stem: Stem of the parquet filename.
            cache: Dict used to memoize resolved keys (mutated in-place).

        Returns:
            A monomer-pair key string, e.g. ``"mt1_mt2"``.
        """
        if self._ckpt_db is not None:
            candidates = [filename_stem]
            if "_x_" in filename_stem:
                candidates.append(filename_stem.replace("_x_", "+"))
            for cand in candidates:
                ck = (p_class, cand)
                if ck in cache:
                    self._db_hit_count += 1
                    return cache[ck]
                try:
                    row = self._ckpt_db.execute(
                        "SELECT mon_pair_key FROM checkpoints "
                        "WHERE p_class=? AND batch_key=?",
                        (p_class, cand),
                    ).fetchone()
                except Exception:
                    row = None
                if row:
                    cache[ck] = row[0]
                    self._db_hit_count += 1
                    return row[0]
        self._db_miss_count += 1
        pairs = self._class_mon_pairs.get(p_class, [])
        for mt1, mt2 in pairs:
            mpk = f"{mt1}_{mt2}"
            if "_x_" in filename_stem:
                parts = filename_stem.split("_x_")
                if len(parts) == 2 and parts[0].startswith(mt1) and parts[1].startswith(mt2):
                    return mpk
            elif filename_stem.startswith(mt1 + "_"):
                return mpk
        return self._mon_pair_set.get(p_class, p_class)

    # ------------------------------------------------------------------
    # Progress helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _fmt_duration(seconds):
        """Format a duration in seconds as a human-readable string.

        Args:
            seconds: Duration in seconds (may be negative).

        Returns:
            A formatted string like ``"2h30m"``, ``"5m12s"``, or
            ``"45s"``.  Returns ``"--"`` for negative values.
        """
        if seconds < 0:
            return "--"
        m, s = divmod(int(seconds), 60)
        h, m = divmod(m, 60)
        if h:
            return f"{h}h{m:02d}m"
        if m:
            return f"{m}m{s:02d}s"
        return f"{s}s"

    def _log_progress_line(self, completed, total, valid, invalid, rows, elapsed,
                           label=""):
        """Log a one-line progress summary at the INFO level.

        Args:
            completed: Number of tasks completed so far.
            total: Total number of tasks.
            valid: Cumulative valid rows.
            invalid: Cumulative invalid rows.
            rows: Cumulative rows processed.
            elapsed: Elapsed wall-clock seconds.
            label: Optional trailing label (e.g. heartbeat indicator).
        """
        pct = (completed / total * 100) if total > 0 else 0
        if completed > 0 and elapsed > 0:
            row_rate = rows / elapsed
            remaining = (total - completed) / (completed / elapsed)
        else:
            row_rate = 0
            remaining = -1
        parts = [
            f"[{completed}/{total} {pct:.1f}%]",
            f"rows={rows:,}", f"ok={valid:,}",
        ]
        if invalid:
            parts.append(f"bad={invalid:,}")
        parts.append(f"{row_rate:,.0f} r/s")
        parts.append(f"elapsed={self._fmt_duration(elapsed)}")
        parts.append(f"eta={self._fmt_duration(remaining)}")
        if label:
            parts.append(label)
        self.logger.info(" | ".join(parts))

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def run(self, p_classes=None):
        """Execute validation across all MPI ranks.

        Args:
            p_classes: Optional iterable of polymer class names to
                restrict validation to.  When ``None`` all discovered
                tasks are processed.

        Returns:
            A dict with summary counters on the master rank;
            ``None`` on worker ranks.
        """
        if self._is_master:
            if p_classes:
                p_set = set(p_classes)
                self._task_list = [
                    t for t in self._task_list if t[2] in p_set
                ]
                self.logger.info("Filtered to %d tasks in: %s",
                                 len(self._task_list), p_classes)
            return self._run_master()
        else:
            self._run_worker()
            return None

    # ------------------------------------------------------------------
    # Master orchestration
    # ------------------------------------------------------------------

    def _run_master(self):
        """Master orchestration loop.

        Distributes tasks to workers, collects and aggregates results,
        runs post-processing merges, writes statistics, and returns a
        summary dictionary.

        Returns:
            A dict with keys ``completed``, ``total``, ``total_rows``,
            ``valid``, ``invalid``, and ``unique_canon``.
        """
        import time as _time
        from mpi4py import MPI  # noqa: F811

        t_start = _time.monotonic()
        self.logger.info("Starting MPI parallel validation run")
        total_tasks = len(self._task_list)
        if total_tasks == 0:
            self.logger.info("No files to validate.")
            self.mpi.shutdown_workers()
            return {"completed": 0, "total": 0}
        self.logger.info("Files: %d  |  Workers: %d",
                         total_tasks, self._n_workers)

        completed = total_valid = total_invalid = 0
        total_unique_canon = total_rows_processed = 0

        class_stats = defaultdict(lambda: {
            "valid": 0, "invalid": 0, "unique_canon": 0,
            "total_rows": 0, "linear": 0,
        })
        mon_pair_stats = defaultdict(lambda: {
            "valid": 0, "invalid": 0, "unique_canon": 0,
            "total_rows": 0, "linear": 0,
        })
        global_nested: Dict[int, int] = defaultdict(int)
        class_nested = defaultdict(lambda: defaultdict(int))
        mon_pair_nested = defaultdict(lambda: defaultdict(int))
        global_errors: Dict[str, int] = defaultdict(int)
        class_errors = defaultdict(lambda: defaultdict(int))
        _invalid_csv_paths: List[Tuple[str, str]] = []
        _linear_csv_paths: List[Tuple[str, str]] = []

        active_workers: Dict[int, Tuple] = {}
        task_idx = 0

        use_rich = self._use_rich
        progress = None
        overall_task = None
        stall_timeout = 120.0
        next_log_at = 1
        last_log_time = t_start

        if use_rich:
            try:
                if self._console is None:
                    from rich.console import Console
                    self._console = Console(
                        stderr=True, force_terminal=True, force_interactive=True,
                    )
                from rich.progress import (
                    BarColumn, MofNCompleteColumn, Progress,
                    TaskProgressColumn, TextColumn,
                    TimeElapsedColumn, TimeRemainingColumn,
                )
                progress = Progress(
                    TextColumn("[progress.description]{task.description:<20}"),
                    BarColumn(), TaskProgressColumn(), TextColumn("•"),
                    MofNCompleteColumn(), TextColumn("•"),
                    TextColumn("[green]✓ {task.fields[ok]}"),
                    TextColumn("[red]✗ {task.fields[bad]}"),
                    TextColumn("•"),
                    TextColumn("[cyan]{task.fields[rows]:,} rows"),
                    TextColumn("•"),
                    TimeElapsedColumn(), TextColumn("<"), TimeRemainingColumn(),
                    console=self._console, expand=True,
                    refresh_per_second=4, transient=True,
                )
                progress.start()
                overall_task = progress.add_task(
                    "[bold green]Validating",
                    total=total_tasks, ok=0, bad=0, rows=0,
                )
            except Exception:
                use_rich = False

        comm = self._comm
        for worker in range(1, self._size):
            if task_idx < total_tasks:
                t = self._task_list[task_idx]
                comm.send(t, dest=worker, tag=TAG_WORK)
                active_workers[worker] = (task_idx, t)
                task_idx += 1

        while active_workers:
            status_obj = MPI.Status()
            result = comm.recv(
                source=MPI.ANY_SOURCE, tag=TAG_RESULT, status=status_obj,
            )
            source = status_obj.Get_source()
            _, completed_task = active_workers.pop(source)
            completed += 1

            if isinstance(result, dict):
                r = result
                total_rows_processed += r.get("total_rows", 0)
                v = r.get("valid", 0)
                iv = r.get("invalid", 0)
                uc = r.get("unique_canon", 0)
                p_class = r.get("p_class", "unknown")
                mpk = r.get("mon_pair_key", "unknown")

                total_valid += v
                total_invalid += iv
                total_unique_canon += uc

                class_stats[p_class]["valid"] += v
                class_stats[p_class]["invalid"] += iv
                class_stats[p_class]["unique_canon"] += uc
                class_stats[p_class]["total_rows"] += r.get("total_rows", 0)
                class_stats[p_class]["linear"] += r.get("linear_count", 0)

                mon_pair_stats[mpk]["valid"] += v
                mon_pair_stats[mpk]["invalid"] += iv
                mon_pair_stats[mpk]["unique_canon"] += uc
                mon_pair_stats[mpk]["total_rows"] += r.get("total_rows", 0)
                mon_pair_stats[mpk]["linear"] += r.get("linear_count", 0)

                for k, c in r.get("nested_counts", {}).items():
                    global_nested[k] += c
                    class_nested[p_class][k] += c
                    mon_pair_nested[mpk][k] += c
                for k, c in r.get("error_counts", {}).items():
                    global_errors[k] += c
                    class_errors[p_class][k] += c

                ic = r.get("invalid_csv")
                if ic:
                    _invalid_csv_paths.append((p_class, ic))
                lc = r.get("linear_parquet")
                if lc:
                    _linear_csv_paths.append((p_class, lc))
            else:
                total_invalid += 1

            elapsed = _time.monotonic() - t_start
            now = _time.monotonic()

            if progress is not None and overall_task is not None:
                progress.update(
                    overall_task, advance=1,
                    ok=total_valid, bad=total_invalid,
                    rows=total_rows_processed,
                )
            else:
                if completed >= next_log_at:
                    self._log_progress_line(
                        completed, total_tasks, total_valid, total_invalid,
                        total_rows_processed, elapsed,
                    )
                    next_log_at = completed + self._log_interval
                    last_log_time = now
                elif now - last_log_time > stall_timeout:
                    self._log_progress_line(
                        completed, total_tasks, total_valid, total_invalid,
                        total_rows_processed, elapsed,
                        label=f"(heartbeat — {len(active_workers)} busy)",
                    )
                    last_log_time = now

            if task_idx < total_tasks:
                nt = self._task_list[task_idx]
                comm.send(nt, dest=source, tag=TAG_WORK)
                active_workers[source] = (task_idx, nt)
                task_idx += 1

        self.mpi.shutdown_workers()
        if progress is not None:
            progress.stop()

        elapsed = _time.monotonic() - t_start
        self._merge_invalid_csvs(_invalid_csv_paths)
        self._merge_linear_parquets(_linear_csv_paths)
        self._write_statistics(
            total_tasks, total_rows_processed, total_valid, total_invalid,
            total_unique_canon,
            dict(class_stats), dict(mon_pair_stats),
            dict(global_nested), dict(class_nested), dict(mon_pair_nested),
            dict(global_errors), dict(class_errors), round(elapsed, 1),
        )

        self.logger.info(
            "Validation finished. files=%d/%d rows=%s valid=%s "
            "invalid=%s unique_canon=%s elapsed=%s",
            completed, total_tasks,
            f"{total_rows_processed:,}", f"{total_valid:,}",
            f"{total_invalid:,}", f"{total_unique_canon:,}",
            self._fmt_duration(elapsed),
        )

        return {
            "completed": completed, "total": total_tasks,
            "total_rows": total_rows_processed,
            "valid": total_valid, "invalid": total_invalid,
            "unique_canon": total_unique_canon,
        }

    # ------------------------------------------------------------------
    # Worker loop
    # ------------------------------------------------------------------

    def _run_worker(self):
        """Worker receive-validate-send loop.

        Blocks on ``comm.recv`` for work items from the master.  Each
        item is validated via ``_validate_file`` and the result dict is
        sent back.  Exits when it receives a ``None`` (stop) message.
        """
        from mpi4py import MPI  # noqa: F811

        comm = self._comm
        while True:
            data = comm.recv(source=0, tag=MPI.ANY_TAG)
            if data is None:
                break
            try:
                result = _validate_file(data)
            except Exception as exc:
                input_path, output_path, p_class, mon_pair_key = data
                result = {
                    "p_class": p_class,
                    "mon_pair_key": mon_pair_key,
                    "batch_key": Path(input_path).stem,
                    "total_rows": 0, "valid": 0, "invalid": 0,
                    "unique_canon": 0, "nested_counts": {},
                    "error_counts": {f"<E-worker-crash|>": 1},
                    "invalid_csv": None, "linear_parquet": None,
                    "linear_count": 0, "elapsed_s": 0.0,
                }
            comm.send(result, dest=0, tag=TAG_RESULT)

    # ------------------------------------------------------------------
    # Post-processing
    # ------------------------------------------------------------------

    def _merge_invalid_csvs(self, csv_paths):
        """Merge per-batch invalid-row CSVs into per-class CSV files.

        Reads all temporary invalid CSVs, concatenates them grouped by
        polymer class, and writes ``<class>_invalid.csv`` under
        ``<output_dir>/invalid/``.  Temporary files and directories are
        cleaned up afterward.

        Args:
            csv_paths: List of ``(p_class, path)`` tuples.
        """
        if not csv_paths:
            return
        by_class: Dict[str, List[str]] = defaultdict(list)
        for p_class, pth in csv_paths:
            by_class[p_class].append(pth)
        for p_class, paths in by_class.items():
            dfs = []
            for p in paths:
                try:
                    dfs.append(pd.read_csv(p))
                except Exception:
                    pass
            if dfs:
                merged = pd.concat(dfs, ignore_index=True)
                out_dir = os.path.join(self.output_dir, "invalid")
                os.makedirs(out_dir, exist_ok=True)
                merged.to_csv(
                    os.path.join(out_dir, f"{p_class}_invalid.csv"), index=False,
                )
                for p in paths:
                    try:
                        os.remove(p)
                    except OSError:
                        pass
                try:
                    os.rmdir(os.path.dirname(paths[0]))
                except OSError:
                    pass

    def _merge_linear_parquets(self, pq_paths):
        """Merge per-batch linear-only parquet files into per-class files.

        Reads all temporary linear parquet files, concatenates them
        grouped by polymer class, and writes ``<class>_linear.parquet``
        under ``<output_dir>/linear_sets/``.  Temporary files and
        directories are cleaned up afterward.

        Args:
            pq_paths: List of ``(p_class, path)`` tuples.
        """
        if not pq_paths:
            return
        by_class: Dict[str, List[str]] = defaultdict(list)
        for p_class, pth in pq_paths:
            by_class[p_class].append(pth)
        for p_class, paths in by_class.items():
            dfs = []
            for p in paths:
                try:
                    dfs.append(pd.read_parquet(p))
                except Exception:
                    pass
            if dfs:
                merged = pd.concat(dfs, ignore_index=True)
                out_dir = os.path.join(self.output_dir, "linear_sets")
                os.makedirs(out_dir, exist_ok=True)
                merged.to_parquet(
                    os.path.join(out_dir, f"{p_class}_linear.parquet"),
                    index=False,
                )
                for p in paths:
                    try:
                        os.remove(p)
                    except OSError:
                        pass
                try:
                    os.rmdir(os.path.dirname(paths[0]))
                except OSError:
                    pass

    def _write_statistics(
        self, total_tasks, total_rows, total_valid, total_invalid,
        total_unique_canon, class_stats, mon_pair_stats,
        nested_counts, class_nested, mon_pair_nested,
        error_counts, class_errors, elapsed_s,
    ):
        """Compute and write per-class and per-monomer-pair statistics.

        Writes ``validation_stats.json`` under
        ``<output_dir>/statistics/`` with global, per-class, and
        per-monomer-pair breakdowns including validity rates, nested
        depth distributions, and the top 30 error types.

        Args:
            total_tasks: Total number of task files processed.
            total_rows: Total rows across all tasks.
            total_valid: Cumulative valid rows.
            total_invalid: Cumulative invalid rows.
            total_unique_canon: Cumulative unique canonical pairs.
            class_stats: Per-class aggregated counters.
            mon_pair_stats: Per-monomer-pair aggregated counters.
            nested_counts: Global nested depth distribution.
            class_nested: Per-class nested depth distribution.
            mon_pair_nested: Per-monomer-pair nested depth distribution.
            error_counts: Global error-type counters.
            class_errors: Per-class error-type counters.
            elapsed_s: Total wall-clock duration in seconds.
        """
        stats_dir = os.path.join(self.output_dir, "statistics")
        os.makedirs(stats_dir, exist_ok=True)

        stats: Dict[str, Any] = {
            "total_tasks": total_tasks,
            "total_rows": total_rows,
            "valid": total_valid,
            "invalid": total_invalid,
            "unique_canonical": total_unique_canon,
            "elapsed_s": elapsed_s,
            "per_class": {},
            "per_monomer_pair": {},
            "nested_depth_distribution": {},
            "top_errors": {},
        }

        for p_class, cs in sorted(class_stats.items()):
            t = cs["valid"] + cs["invalid"]
            stats["per_class"][p_class] = {
                "total_rows": cs["total_rows"],
                "valid": cs["valid"],
                "invalid": cs["invalid"],
                "unique_canonical": cs["unique_canon"],
                "linear": cs["linear"],
                "valid_pct": round(cs["valid"] / max(t, 1) * 100, 2),
            }

        for mpk, ms in sorted(mon_pair_stats.items()):
            t = ms["valid"] + ms["invalid"]
            stats["per_monomer_pair"][mpk] = {
                "total_rows": ms["total_rows"],
                "valid": ms["valid"],
                "invalid": ms["invalid"],
                "unique_canonical": ms["unique_canon"],
                "valid_pct": round(ms["valid"] / max(t, 1) * 100, 2),
            }

        for d, c in sorted(nested_counts.items()):
            stats["nested_depth_distribution"][f"depth_{d}"] = c

        top = sorted(error_counts.items(), key=lambda x: -x[1])[:30]
        stats["top_errors"] = {e: c for e, c in top}

        with open(os.path.join(stats_dir, "validation_stats.json"), "w") as fh:
            json.dump(stats, fh, indent=2)

        self.logger.info("Statistics written to %s/", stats_dir)
