"""
Base class for MPI-parallel executors.

Composes ``MpiProtocol``, ``MpiLogger``, ``ProgressTracker``, and
``ResourceMonitor`` so downstream MPI classes only need to implement
their domain-specific ``run()`` and ``_run_worker()`` methods.

Usage
-----
    from polyscript.utils.mpi_executor import BaseMPIExecutor

    class MyRunner(BaseMPIExecutor):
        def run(self) -> dict:
            if self._is_master:
                tasks = [...]
                chunks = self.splitter.round_robin(tasks, self._size)
                ...
            else:
                self._run_worker()
            return {}

        def _run_worker(self):
            while True:
                task = self.mpi.recv_work()
                if task is None:  # poison pill
                    break
                result = process(task)
                self.mpi.send_result(result)
"""

from __future__ import annotations

from abc import abstractmethod
from typing import Any, Dict

from polyscript.utils.executor import BaseExecutor
from polyscript.utils.logger import INFO
from polyscript.utils.mpi import (
    MpiLogger,
    MpiProtocol,
    ProgressTracker,
    ResourceMonitor,
    WorkSplitter,
)


class BaseMPIExecutor(BaseExecutor):
    """Base class for MPI-parallel runners — composition over inheritance.

    Provides ready-to-use components:

    * ``self.mpi`` — :class:`MpiProtocol` for scatter/gather/messaging
    * ``self.logger`` — :class:`MpiLogger` with Rich + direct-stderr
    * ``self.splitter`` — :class:`WorkSplitter` for round-robin/chunks
    * ``self.progress`` — :class:`ProgressTracker` (created per-run)
    * ``self.resources`` — :class:`ResourceMonitor` background sampler

    Subclasses implement :meth:`run` and :meth:`_run_worker`.

    Parameters
    ----------
    log_filepath : str, optional
        Detailed log file path (master only).
    progress_file : str, optional
        ``tail -f``-friendly progress file.
    log_interval : int
        Minimum work units between progress log lines (default 100).
    use_rich : bool
        Enable Rich console handler (default ``True``).
    enable_resource_monitor : bool
        Start RSS sampling thread (requires psutil).
    """

    _log_level = INFO

    def __init__(
        self,
        *,
        log_filepath: str | None = None,
        progress_file: str | None = None,
        log_interval: int = 100,
        use_rich: bool = True,
        enable_resource_monitor: bool = True,
        **kwargs: Any,
    ) -> None:
        # ── MPI init ──────────────────────────────────────────────────
        from mpi4py import MPI as _MPI  # noqa: N813

        comm = _MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()
        is_master = rank == 0

        if size < 2:
            cls_name = type(self).__name__
            raise RuntimeError(
                f"{cls_name} needs at least 2 MPI ranks (1 master + ≥1 worker). "
                f"Got {size} rank(s). "
                "Did you forget to launch with mpirun?"
            )

        self._comm = comm
        self._rank = rank
        self._size = size
        self._is_master = is_master

        # ── Components ─────────────────────────────────────────────────
        self.mpi = MpiProtocol(comm, rank, size, is_master)
        self.splitter = WorkSplitter()

        # Logger: master gets full setup, workers minimal
        if is_master:
            self.logger = MpiLogger(
                name=type(self).__name__.lower(),
                level=self._log_level,
                log_filepath=log_filepath,
                progress_file=progress_file,
                use_rich=use_rich,
                use_direct_stderr=True,
            )
        else:
            # Workers need a logger too (errors show on stderr)
            super().__init__(**kwargs)  # sets self.logger via BaseExecutor

        self._log_interval = max(log_interval, 1)
        self._use_rich = use_rich

        # Resource monitor (master only)
        self.resources = ResourceMonitor(
            enabled=enable_resource_monitor and is_master
        )

        # Progress tracker — created lazily in run()
        self.progress: ProgressTracker | None = None

        # Synchronise
        self.mpi.barrier()

    # ------------------------------------------------------------------
    # Lifecycle
    # ------------------------------------------------------------------

    def close(self) -> None:
        self.resources.stop()
        if hasattr(self, "logger"):
            self.logger.flush()
            self.logger.close()

    def __enter__(self) -> "BaseMPIExecutor":
        return self

    def __exit__(self, *args: Any) -> None:
        self.close()

    # ------------------------------------------------------------------
    # Abstract interface
    # ------------------------------------------------------------------

    @abstractmethod
    def run(self) -> Dict[str, Any]:
        """Execute the MPI-parallel workload on all ranks."""
        ...

    @abstractmethod
    def _run_worker(self) -> None:
        """Worker loop: receive work, process, send results."""
        ...
