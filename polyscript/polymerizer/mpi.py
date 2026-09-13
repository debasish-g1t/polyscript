"""
MPI-parallel polymerizer.

``MpiPolymerizer`` is a distributed batch-polymerisation runner that
composes from ``polyscript.utils.mpi`` components (MpiLogger, CheckpointDB,
ProgressFile, ResourceMonitor) with Rich progress bars and per-class
statistics.

Usage (launch with mpirun)::

    from polyscript.polymerizer.mpi import MpiPolymerizer

    runner = MpiPolymerizer(
        exp_name="my_run",
        input_data_split_dir="./splits",
    )
    stats = runner.run()
"""

from __future__ import annotations

import os
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import joblib

from polyscript.utils.mpi import (
    CheckpointDB,
    MpiLogger,
    MpiProtocol,
    ProgressFile,
    ProgressTracker,
    ResourceMonitor,
    WorkSplitter,
)

SKIP_CLS = ["COC", "rec:radi", "rec:cati", "rec:ani", "rec:coord"]

# Auto-resolve the rules directory shipped with this package
_HERE = Path(__file__).resolve().parent
_DEFAULT_RULE_DIR = str(_HERE / "rules") if (_HERE / "rules").is_dir() else "./rules"


# ---------------------------------------------------------------------------
# Worker function  (top-level, required by MPI/spawn)
# ---------------------------------------------------------------------------


def _polymerize_batch(
    args: Tuple,
) -> Tuple[str, str, str, bool, Optional[str], float]:
    """Process a single polymerisation batch (top-level worker function).

    Instantiates a ``Polymerizer``, executes one config, and returns the
    result along with peak memory usage.

    Args:
        args (tuple): A ``(config, p_class, mon_pair_key, batch_key)``
            tuple where *config* is a dict as accepted by
            ``Polymerizer.execute``.

    Returns:
        tuple: ``(p_class, mon_pair_key, batch_key, success, error, peak_rss_mb)``
        where *success* is a bool, *error* is ``None`` or an error string,
        and *peak_rss_mb* is the peak RSS memory in MiB (or ``-1.0`` if
        ``psutil`` is not available).
    """
    from polyscript.polymerizer import Polymerizer

    config, p_class, mon_pair_key, batch_key = args

    try:
        import psutil
        proc = psutil.Process()
        mem_before = proc.memory_info().rss
    except ImportError:
        psutil = None
        mem_before = 0

    polymerizer = Polymerizer()
    success, error = polymerizer.execute(config)

    if psutil is not None:
        mem_after = proc.memory_info().rss
        peak_rss_mb = max(mem_before, mem_after) / (1024 * 1024)
    else:
        peak_rss_mb = -1.0

    return (p_class, mon_pair_key, batch_key, success, error, peak_rss_mb)


# ---------------------------------------------------------------------------
# MpiPolymerizer
# ---------------------------------------------------------------------------


class MpiPolymerizer:
    """Distributed batch-polymerisation runner.

    .. important::

        Must be launched via ``mpirun`` with at least **2 MPI ranks**
        (1 master + 1 or more workers).  Running as a plain Python
        script raises ``RuntimeError``::

            mpirun -np 4 python my_script.py

    Components (all from ``polyscript.utils.mpi``):
    * ``self.logger``      — MpiLogger (direct-stderr + file + Rich)
    * ``self.mpi``         — MpiProtocol (scatter/gather/messaging)
    * ``self.checkpoint``  — CheckpointDB (SQLite task tracking)
    * ``self.progress_file`` — ProgressFile (atomic JSON snapshot)
    * ``self.resources``   — ResourceMonitor (background RSS sampler)
    * ``self.splitter``    — WorkSplitter (round-robin / chunks)
    """

    def __init__(
        self,
        exp_name: str,
        input_data_split_dir: str,
        *,
        output_path: str = "./outputs",
        ckpt_path: str = "ckpt.db",
        rule_dir: Optional[str] = None,
        log_filepath: str = "polymerizer.log",
        col_name: str = "psmiles",
        chunk_size: Optional[int] = None,
        log_interval: int = 100,
        use_rich: bool = True,
    ):
        """Initialize the MPI polymeriser.

        Args:
            exp_name (str): Experiment name; used to create the output
                subdirectory ``output_path/exp_name``.
            input_data_split_dir (str): Directory containing pre-split
                monomer parquet files.
            output_path (str): Base directory for outputs.
            ckpt_path (str): Filename for the SQLite checkpoint database.
            rule_dir (str or None): Path to the rule directory.  Defaults
                to the bundled ``rules/`` directory.
            log_filepath (str): Path for the log file.
            col_name (str): Column name used for SMILES strings in the
                monomer DataFrames.
            chunk_size (int or None): Number of tasks to dispatch per
                fetch round.  Defaults to ``max(n_workers * 4, 50)``.
            log_interval (int): How often (in completed tasks) to emit
                progress log lines.
            use_rich (bool): Enable Rich progress bars on the master rank.

        Raises:
            RuntimeError: If there are fewer than 2 MPI ranks.
            ValueError: If the rule directory or input split directory
                does not exist.
        """
        # ── MPI init ──────────────────────────────────────────────────
        from mpi4py import MPI as _MPI  # noqa: N813

        comm = _MPI.COMM_WORLD
        self._rank = comm.Get_rank()
        self._size = comm.Get_size()
        self._is_master = self._rank == 0

        if self._size < 2:
            raise RuntimeError(
                "MpiPolymerizer needs at least 2 MPI ranks (1 master + >=1 worker). "
                f"Got {self._size} rank(s).  Launch with mpirun."
            )
        self._n_workers = self._size - 1

        if rule_dir is None:
            rule_dir = _DEFAULT_RULE_DIR

        # ── Components ────────────────────────────────────────────────
        self.mpi = MpiProtocol(comm, self._rank, self._size, self._is_master)
        self.splitter = WorkSplitter()

        if self._is_master:
            self.exp_name = exp_name
            self.exp_path = os.path.join(output_path, self.exp_name)
            os.makedirs(self.exp_path, exist_ok=True)

            self.logger = MpiLogger(
                name="polymerizer",
                log_filepath=os.path.join(self.exp_path, log_filepath),
                progress_file=os.path.join(self.exp_path, "progress.txt"),
                use_rich=use_rich,
                use_direct_stderr=True,
            )
            self._use_rich = use_rich
            self._console = self.logger._console if use_rich else None
            self._log_interval = max(log_interval, 1)

            if not os.path.exists(rule_dir):
                raise ValueError(f"Rule directory {rule_dir} does not exist.")
            if not os.path.exists(input_data_split_dir):
                raise ValueError(
                    f"Input data split directory {input_data_split_dir} does not exist."
                )
            self.rule_dir = rule_dir
            self.input_data_split_dir = input_data_split_dir
            self.col_name = col_name

            self._chunk_size = (
                chunk_size if chunk_size is not None else max(self._n_workers * 4, 50)
            )

            # ── Checkpoint DB ─────────────────────────────────────────
            ckpt_db_path = os.path.join(self.exp_path, ckpt_path)
            self.checkpoint = CheckpointDB(ckpt_db_path)
            self.checkpoint.init_schema()
            self.progress_file = ProgressFile(
                os.path.join(self.exp_path, "progress.json")
            )

            # ── Migrate legacy JSON checkpoint ────────────────────────
            legacy_json = os.path.join(self.exp_path, "ckpt.json")
            if os.path.exists(legacy_json) and self.checkpoint.count_total() == 0:
                self.logger.info("Migrating legacy JSON checkpoint …")
                self.checkpoint.migrate_json(legacy_json)
                self.logger.info("Migration complete")
            elif self.checkpoint.count_total() == 0:
                self.logger.info("Populating checkpoint DB …")
                self._populate_checkpoint()
            else:
                self.logger.info(
                    "Checkpoint DB loaded: %d pending / %d total",
                    self.checkpoint.count_pending(),
                    self.checkpoint.count_total(),
                )

            self.logger.info("Using %d MPI worker ranks", self._n_workers)
            self.logger.info("Task-fetch chunk size = %d", self._chunk_size)

            self.resources = ResourceMonitor(enabled=True)
            self._worker_peaks: List[float] = []
        else:
            import logging
            self.logger = logging.getLogger(__name__ + f".rank{self._rank}")
            self.logger.setLevel(logging.WARNING)

        self.mpi.barrier()

    # ------------------------------------------------------------------
    # Lifecycle
    # ------------------------------------------------------------------

    def close(self):
        """Close the checkpoint database and flush/close the logger."""
        if self._is_master:
            if hasattr(self, "checkpoint"):
                self.checkpoint.close()
            if hasattr(self, "logger"):
                try:
                    self.logger.flush()
                    self.logger.close()
                except AttributeError:
                    pass

    def __enter__(self):
        """Enter the runtime context (returns ``self``)."""
        return self

    def __exit__(self, *args):
        """Exit the runtime context and close resources."""
        self.close()

    # ──────────────────────────────────────────────────────────────────
    # Monomer / split helpers
    # ──────────────────────────────────────────────────────────────────

    def _get_cls_split(self, cls_name: str) -> List[str]:
        """Find parquet split files for a given monomer class.

        Args:
            cls_name (str): Monomer class name (e.g. ``"vinyl"``).

        Returns:
            list of str: The first element is the directory name, and
            subsequent elements are split file basenames (without
            ``.parquet`` extension).  Returns an empty list if no
            matching directory is found.
        """
        cls_split: List[str] = []
        valid_dirs = os.listdir(self.input_data_split_dir)
        for dirs in valid_dirs:
            if "_".join(dirs.split("_")[:-1]) == cls_name:
                cls_split.append(dirs)
        if cls_split == []:
            return cls_split
        assert len(cls_split) == 1, (
            f"Expected exactly one directory starting with {cls_name}, "
            f"found {len(cls_split)}"
        )
        for files in os.listdir(
            os.path.join(self.input_data_split_dir, cls_split[0])
        ):
            if files.endswith(".parquet"):
                cls_split.append(files.split(".")[0])
        return cls_split

    def _permute_mon1_mon2(
        self, mon1_splits: List[str], mon2_splits: List[str]
    ) -> List[str]:
        """Generate all cross-product pair keys from two split lists.

        Args:
            mon1_splits (list of str): Split IDs for monomer type 1.
            mon2_splits (list of str): Split IDs for monomer type 2.
                If empty, returns *mon1_splits* unchanged.

        Returns:
            list of str: Strings of the form ``"mon1+mon2"`` for every
            combination.
        """
        if len(mon2_splits) == 0:
            return mon1_splits
        permuted: List[str] = []
        for mon1 in mon1_splits:
            for mon2 in mon2_splits:
                permuted.append(f"{mon1}+{mon2}")
        return permuted

    # ──────────────────────────────────────────────────────────────────
    # Checkpoint population
    # ──────────────────────────────────────────────────────────────────

    def _populate_checkpoint(self):
        """Scan rules and input splits, then populate the checkpoint DB.

        Loads ``ps_gen.pkl``, enumerates all (class, mon_type1, mon_type2)
        rules, discovers their corresponding parquet splits, and inserts
        tasks into the checkpoint database.
        """
        ps_cls = joblib.load(os.path.join(self.rule_dir, "ps_gen.pkl"))
        rows: List[Tuple] = []
        for k, v in ps_cls.items():
            if k in SKIP_CLS:
                continue
            for values in v:
                mon_type1 = values[0]
                mon_type2 = values[1]
                mon1_splits = self._get_cls_split(mon_type1)
                if len(mon1_splits) in [0, 1]:
                    continue
                mon2_splits: list = []
                if mon_type2 != "none":
                    mon2_splits = self._get_cls_split(mon_type2)
                permuted_list = self._permute_mon1_mon2(
                    mon1_splits[1:], mon2_splits[1:]
                )
                mon_pair_key = f"{mon_type1}_{mon_type2}"
                self.logger.info(
                    "Permutations for %s: %d", mon_pair_key, len(permuted_list)
                )
                for batch_key in permuted_list:
                    rows.append((k, mon_pair_key, mon_type1, mon_type2, batch_key, 0))
        self.checkpoint.populate(rows)
        self.logger.info("Checkpoint DB populated with %d tasks", len(rows))

    # ──────────────────────────────────────────────────────────────────
    # Config / task collection
    # ──────────────────────────────────────────────────────────────────

    def _parse_batch_key(
        self, mon_type1: str, mon_type2: str, batch_key: str
    ) -> Tuple[Optional[str], Optional[str]]:
        """Resolve a batch key into concrete parquet file paths.

        Args:
            mon_type1 (str): Monomer type 1.
            mon_type2 (str): Monomer type 2.
            batch_key (str): Batch identifier, either a single split ID
                or ``"id1+id2"`` for a pair.

        Returns:
            tuple: ``(mon1_path, mon2_path)`` where each is a full path
            or ``None`` if the key could not be resolved.
        """
        if "+" in batch_key:
            parts = batch_key.split("+")
            if len(parts) != 2:
                return None, None
            b1, b2 = parts
        else:
            b1, b2 = batch_key, None

        cls1 = self._get_cls_split(mon_type1)
        if len(cls1) < 2:
            return None, None
        dir1 = cls1[0]
        try:
            cls1[1:].index(b1)
        except ValueError:
            return None, None
        mon1_path = os.path.join(self.input_data_split_dir, dir1, f"{b1}.parquet")

        if b2 is None:
            return mon1_path, None

        cls2 = self._get_cls_split(mon_type2)
        if len(cls2) < 2:
            return None, None
        dir2 = cls2[0]
        try:
            cls2[1:].index(b2)
        except ValueError:
            return None, None
        mon2_path = os.path.join(self.input_data_split_dir, dir2, f"{b2}.parquet")
        return mon1_path, mon2_path

    def _build_config(
        self, p_class: str, mon_type1: str, mon_type2: str, batch_key: str
    ) -> Optional[Dict]:
        """Build a configuration dict for a single task.

        Args:
            p_class (str): Polymer class label.
            mon_type1 (str): Monomer type 1.
            mon_type2 (str): Monomer type 2.
            batch_key (str): Batch identifier.

        Returns:
            dict or None: Configuration dict with keys ``mon1_path``,
            ``mon2_path``, ``output_path``, ``col_name``, ``P_class``,
            ``mon_type1``, ``mon_type2``, or ``None`` if paths can't be
            resolved.
        """
        mon1_path, mon2_path = self._parse_batch_key(mon_type1, mon_type2, batch_key)
        if mon1_path is None:
            return None

        output_dir = os.path.join(self.exp_path, p_class)
        os.makedirs(output_dir, exist_ok=True)
        output_filename = (
            f"{mon_type1}_{mon_type2}_{batch_key.replace('+', '_x_')}.parquet"
        )
        output_path = os.path.join(output_dir, output_filename)

        return {
            "mon1_path": mon1_path,
            "mon2_path": mon2_path,
            "output_path": output_path,
            "col_name": self.col_name,
            "P_class": p_class,
            "mon_type1": mon_type1,
            "mon_type2": mon_type2 if mon_type2 != "none" else None,
        }

    def _collect_pending_tasks(
        self, limit: int, p_classes: Optional[List[str]] = None
    ) -> List[Tuple]:
        """Fetch up to *limit* uncompleted tasks from the checkpoint DB.

        Args:
            limit (int): Maximum number of tasks to return.
            p_classes (list of str or None): Optional filter to specific
                polymer classes.

        Returns:
            list of tuple: Each element is
            ``(config, p_class, mon_pair_key, batch_key)`` where *config*
            is the dict returned by ``_build_config``.
        """
        query = (
            "SELECT p_class, mon_pair_key, mon_type1, mon_type2, batch_key "
            "FROM checkpoints WHERE completed = 0"
        )
        params: tuple = ()
        if p_classes:
            placeholders = ", ".join("?" for _ in p_classes)
            query += f" AND p_class IN ({placeholders})"
            params = tuple(p_classes)
        query += " LIMIT ?"
        params += (limit,)
        rows = self.checkpoint.conn.execute(query, params).fetchall()
        tasks = []
        for p_class, mpk, mt1, mt2, bk in rows:
            config = self._build_config(p_class, mt1, mt2, bk)
            if config is None:
                self.checkpoint.mark_complete(p_class, mpk, bk, success=False)
                continue
            tasks.append((config, p_class, mpk, bk))
        return tasks

    # ──────────────────────────────────────────────────────────────────
    # Public API
    # ──────────────────────────────────────────────────────────────────

    def list_classes(self) -> List[str]:
        """Return the distinct polymer classes present in the checkpoint DB.

        Returns:
            list of str: Sorted list of class names.
        """
        return self.checkpoint.list_classes()

    def get_progress(self) -> Dict:
        """Return a summary of overall task progress.

        Returns:
            dict: Keys are ``total``, ``completed``, ``pending``, and
            ``progress_pct`` (float).
        """
        total = self.checkpoint.count_total()
        completed = total - self.checkpoint.count_pending()
        return {
            "total": total,
            "completed": completed,
            "pending": total - completed,
            "progress_pct": (completed / total * 100) if total > 0 else 0,
        }

    def reset_checkpoint(self, p_class: Optional[str] = None):
        """Reset the checkpoint state for a class or all classes.

        Args:
            p_class (str or None): Specific class to reset, or ``None``
                to reset all.
        """
        self.checkpoint.reset(p_class)
        self.logger.info("Reset checkpoint: %s", p_class or "all classes")

    # ──────────────────────────────────────────────────────────────────
    # Main entry point
    # ──────────────────────────────────────────────────────────────────

    def run(
        self,
        show_class_progress: bool = True,
        p_classes: Optional[List[str]] = None,
        exclude_classes: Optional[List[str]] = None,
    ):
        """Start the polymerisation run (entry point for all ranks).

        On the master rank, orchestrates task dispatch and returns a
        summary dict.  On worker ranks, enters the receive/process loop
        and returns ``None``.

        Args:
            show_class_progress (bool): Whether to show per-class Rich
                progress bars.
            p_classes (list of str or None): Specific polymer classes to
                process.  Cannot be used with *exclude_classes*.
            exclude_classes (list of str or None): Polymer classes to
                skip.  Cannot be used with *p_classes*.

        Returns:
            dict or None: On master, a dict with keys ``completed``,
            ``failed``, ``total``.  On workers, ``None``.

        Raises:
            ValueError: If both *p_classes* and *exclude_classes* are
                specified, or if *p_classes* contains unknown class names.
        """
        if self._is_master:
            return self._run_master(
                show_class_progress=show_class_progress,
                p_classes=p_classes,
                exclude_classes=exclude_classes,
            )
        else:
            self._run_worker()
            return None

    # ──────────────────────────────────────────────────────────────────
    # Master orchestration
    # ──────────────────────────────────────────────────────────────────

    def _run_master(
        self,
        show_class_progress: bool = True,
        p_classes: Optional[List[str]] = None,
        exclude_classes: Optional[List[str]] = None,
    ) -> Dict:
        """Orchestrate task dispatch and result collection on the master rank.

        Seeds initial work to all workers, then enters an event loop that
        receives completed results, updates checkpoints and progress
        display, and dispatches new work until all tasks are done.

        Args:
            show_class_progress (bool): Whether to show per-class Rich
                progress bars.
            p_classes (list of str or None): Filter to specific classes.
            exclude_classes (list of str or None): Classes to skip.

        Returns:
            dict: Summary with keys ``completed``, ``failed``, ``total``.

        Raises:
            ValueError: If both *p_classes* and *exclude_classes* are set.
        """
        from mpi4py import MPI  # noqa: F811

        comm = self.mpi.comm
        t_start = time.monotonic()
        self.resources.start()
        self.logger.info("Starting MPI parallel batch polymerisation")

        self.checkpoint.conn.execute("PRAGMA synchronous=OFF")

        if p_classes and exclude_classes:
            raise ValueError("Cannot specify both p_classes and exclude_classes.")

        if exclude_classes:
            known = set(self.list_classes())
            unknown = [c for c in exclude_classes if c not in known]
            if unknown:
                self.logger.warning("Excluded unknown classes: %s", unknown)
            p_classes = sorted(known - set(exclude_classes))
            if not p_classes:
                self.logger.info("All classes excluded. Nothing to process.")
                self.mpi.shutdown_workers()
                return {"completed": 0, "failed": 0, "total": 0}
            self.logger.info("Excluding %s — processing: %s", exclude_classes, p_classes)

        if p_classes:
            known = set(self.list_classes())
            unknown = [c for c in p_classes if c not in known]
            if unknown:
                raise ValueError(
                    f"Unknown class(es): {unknown}. Available: {sorted(known)}"
                )
            self.logger.info("Filtering to classes: %s", p_classes)

        total_tasks = self.checkpoint.count_pending(p_classes=p_classes)
        if total_tasks == 0:
            self.logger.info("No pending tasks — all complete.")
            self.mpi.shutdown_workers()
            return {"completed": 0, "failed": 0, "total": 0}

        self.logger.info(
            "Pending: %d  chunk_size=%d  n_workers=%d",
            total_tasks, self._chunk_size, self._n_workers,
        )

        self.progress = ProgressTracker(total=total_tasks, log_interval=self._log_interval)

        # Per-class counters
        class_progress: Dict[str, Dict[str, int]] = {}
        query = "SELECT p_class, COUNT(*) FROM checkpoints WHERE completed = 0"
        params: tuple = ()
        if p_classes:
            placeholders = ", ".join("?" for _ in p_classes)
            query += f" AND p_class IN ({placeholders})"
            params = tuple(p_classes)
        query += " GROUP BY p_class"
        for p_class, count in self.checkpoint.conn.execute(query, params):
            class_progress[p_class] = {"completed": 0, "failed": 0, "total": count}

        self.progress_file.write(0, 0, total_tasks)

        # Rich progress
        use_rich = self._use_rich
        rich_progress = None
        overall_task = None
        class_tasks: Dict[str, int] = {}

        if use_rich:
            try:
                if self._console is None:
                    from rich.console import Console
                    self._console = Console(
                        force_terminal=True, force_interactive=True, stderr=True,
                    )
                from rich.progress import (
                    BarColumn, MofNCompleteColumn, Progress,
                    TaskProgressColumn, TextColumn,
                    TimeElapsedColumn, TimeRemainingColumn,
                )
                rich_progress = Progress(
                    TextColumn("[progress.description]{task.description:<30}"),
                    BarColumn(), TaskProgressColumn(), TextColumn("•"),
                    MofNCompleteColumn(), TextColumn("•"),
                    TextColumn("[green]✓ {task.fields[ok]:<5}"),
                    TextColumn("[red]✗ {task.fields[bad]:<5}"),
                    TextColumn("•"),
                    TimeElapsedColumn(), TextColumn("<"), TimeRemainingColumn(),
                    console=self._console, expand=True,
                    refresh_per_second=4, transient=True,
                )
                rich_progress.start()
                overall_task = rich_progress.add_task(
                    "[bold green]Overall", total=total_tasks, ok=0, bad=0,
                )
                if show_class_progress:
                    for p_class in sorted(class_progress):
                        class_tasks[p_class] = rich_progress.add_task(
                            f"[blue]  {p_class[:25]}[/]",
                            total=class_progress[p_class]["total"], ok=0, bad=0,
                        )
            except ImportError:
                use_rich = False

        # State
        active_workers: Dict[int, Tuple] = {}
        dispatched: set = set()
        completed = 0
        failed = 0
        tasks_submitted = 0

        # Seed initial tasks
        chunk = self._collect_pending_tasks(
            limit=min(self._chunk_size, total_tasks), p_classes=p_classes,
        )
        for worker in range(1, self._size):
            if chunk:
                task = chunk.pop(0)
                _, p_class, mpk, bk = task
                dispatched.add((p_class, mpk, bk))
                comm.send(task, dest=worker, tag=0)  # TAG_WORK
                active_workers[worker] = task
        tasks_submitted = len(active_workers)

        # Main event loop
        while active_workers:
            status_obj = MPI.Status()
            result = comm.recv(source=MPI.ANY_SOURCE, tag=1, status=status_obj)
            source = status_obj.Get_source()
            task = active_workers.pop(source)
            _, _pc, _mpk, _bk = task
            dispatched.discard((_pc, _mpk, _bk))

            _r_p_class, _r_mpk, _r_bk, success, error = result[:5]
            worker_peak_rss = result[5] if len(result) > 5 else -1.0

            if success:
                self.checkpoint.mark_complete(_r_p_class, _r_mpk, _r_bk, True)
                completed += 1
                class_progress[_r_p_class]["completed"] += 1
            else:
                failed += 1
                class_progress[_r_p_class]["failed"] += 1
                self.logger.warning(
                    "FAILED  %s/%s/%s  error=%s", _r_p_class, _r_mpk, _r_bk, error,
                )

            if worker_peak_rss > 0:
                self._worker_peaks.append(worker_peak_rss)

            self.progress.tick(ok=success)

            # Update display
            if use_rich and rich_progress is not None and overall_task is not None:
                rich_progress.update(
                    overall_task, advance=1,
                    ok=self.progress.completed, bad=failed,
                )
                if show_class_progress and _r_p_class in class_tasks:
                    rich_progress.update(
                        class_tasks[_r_p_class], advance=1,
                        ok=class_progress[_r_p_class]["completed"],
                        bad=class_progress[_r_p_class]["failed"],
                    )
            else:
                if self.progress.should_log():
                    self.logger.info(self.progress.format_line())
                    self.progress.mark_logged()
                elif self.progress.is_stalled():
                    self.logger.info("%s (heartbeat)", self.progress.format_line())
                    self.progress.mark_logged()

            if self.progress.processed % 10 == 0:
                self.progress_file.write(self.progress.completed, failed, total_tasks)

            # Fetch more work
            if not chunk:
                remaining = total_tasks - self.progress.processed
                fetch_limit = min(self._chunk_size, remaining)
                if fetch_limit > 0:
                    raw = self._collect_pending_tasks(limit=fetch_limit, p_classes=p_classes)
                    chunk = [t for t in raw if (t[1], t[2], t[3]) not in dispatched]

            if chunk:
                next_task = chunk.pop(0)
                _, pc2, mpk2, bk2 = next_task
                dispatched.add((pc2, mpk2, bk2))
                comm.send(next_task, dest=source, tag=0)
                active_workers[source] = next_task
                tasks_submitted += 1

        # Drain remaining
        while chunk:
            task = chunk.pop(0)
            for worker in range(1, self._size):
                if worker not in active_workers:
                    _, pc3, mpk3, bk3 = task
                    dispatched.add((pc3, mpk3, bk3))
                    comm.send(task, dest=worker, tag=0)
                    active_workers[worker] = task
                    tasks_submitted += 1
                    break
            status_obj = MPI.Status()
            result = comm.recv(source=MPI.ANY_SOURCE, tag=1, status=status_obj)
            source = status_obj.Get_source()
            task2 = active_workers.pop(source)
            _rpc, _rmpk, _rbk = task2[1], task2[2], task2[3]
            dispatched.discard((_rpc, _rmpk, _rbk))
            success = result[3]
            if success:
                self.checkpoint.mark_complete(_rpc, _rmpk, _rbk, True)
                completed += 1
                class_progress[_rpc]["completed"] += 1
            else:
                failed += 1
                class_progress[_rpc]["failed"] += 1
            self.progress.tick(ok=success)

        # Shutdown
        self.mpi.shutdown_workers()
        if rich_progress is not None:
            rich_progress.stop()

        elapsed = time.monotonic() - t_start
        self.progress_file.write(completed, failed, total_tasks)
        self.resources.stop()

        self.logger.info(
            "Run finished.  ok=%d  fail=%d  total=%d  submitted=%d  elapsed=%s",
            completed, failed, total_tasks, tasks_submitted,
            self.progress._fmt_duration(elapsed),
        )
        self.logger.info("Per-class summary:")
        for p_class, cp in sorted(class_progress.items()):
            status = "OK" if cp["failed"] == 0 else "WARN"
            self.logger.info(
                "  [%s] %s: %d/%d  failed=%d",
                status, p_class, cp["completed"], cp["total"], cp["failed"],
            )

        self.checkpoint.conn.execute("PRAGMA synchronous=NORMAL")
        return {"completed": completed, "failed": failed, "total": total_tasks}

    # ──────────────────────────────────────────────────────────────────
    # Worker loop
    # ──────────────────────────────────────────────────────────────────

    def _run_worker(self):
        """Worker loop: receive tasks, process, and return results.

        Runs in an infinite loop receiving task tuples from the master,
        passing them to ``_polymerize_batch``, and sending results back.
        Exits when a ``None`` shutdown signal is received.
        """
        from mpi4py import MPI  # noqa: F811

        comm = self.mpi.comm
        while True:
            data = comm.recv(source=0, tag=MPI.ANY_TAG)
            if data is None:
                break
            try:
                result = _polymerize_batch(data)
            except Exception as exc:
                config, p_class, mon_pair_key, batch_key = data
                result = (p_class, mon_pair_key, batch_key, False, str(exc), -1.0)
            comm.send(result, dest=0, tag=1)  # TAG_RESULT
