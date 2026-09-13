"""
MPI master↔worker communication protocol.

Encapsulates scatter / gather, tagged messaging, and the poison-pill
shutdown pattern so runners don't need to manage MPI tags directly.
"""

from __future__ import annotations

from typing import Any, List, Optional

# Tags used for the master↔worker protocol
TAG_WORK = 0  # master → worker: a task tuple
TAG_RESULT = 1  # worker → master: result tuple
TAG_STOP = 2  # master → worker: poison pill (None)


class MpiProtocol:
    """Thin wrapper around an MPI communicator for master↔worker messaging.

    Parameters
    ----------
    comm : mpi4py.MPI.Comm
        The MPI communicator (typically ``MPI.COMM_WORLD``).
    rank : int
        This process's rank.
    size : int
        Total number of ranks.
    is_master : bool
        Whether this rank is the master (rank 0).
    """

    def __init__(self, comm: Any, rank: int, size: int, is_master: bool) -> None:
        self.comm = comm
        self.rank = rank
        self.size = size
        self.is_master = is_master
        self.n_workers = size - 1

    # ------------------------------------------------------------------
    # Scatter / gather
    # ------------------------------------------------------------------

    def scatter(
        self, items: List[Any]
    ) -> Optional[List[Any]]:
        """Scatter a list from master to all ranks.

        Master passes the full list; workers receive their chunk.
        Returns ``None`` on workers if master has no data.
        """
        if self.is_master:
            chunks = WorkSplitter.round_robin(items, self.size)
        else:
            chunks = None
        return self.comm.scatter(chunks, root=0)

    def gather(self, items: List[Any]) -> Optional[List[List[Any]]]:
        """Gather lists from all ranks to master.

        Returns ``None`` on workers, a list of per-rank lists on master.
        """
        return self.comm.gather(items, root=0)

    # ------------------------------------------------------------------
    # Tagged point-to-point
    # ------------------------------------------------------------------

    def send_work(self, dest: int, task: Any) -> None:
        """Send a work item to a worker."""
        self.comm.send(task, dest=dest, tag=TAG_WORK)

    def recv_work(self, source: int = 0) -> Any:
        """Receive a work item from the master."""
        return self.comm.recv(source=source, tag=TAG_WORK)

    def send_result(self, result: Any) -> None:
        """Send a result back to the master."""
        self.comm.send(result, dest=0, tag=TAG_RESULT)

    def recv_result(self, source: int) -> Any:
        """Receive a result from a worker."""
        return self.comm.recv(source=source, tag=TAG_RESULT)

    def send_stop(self, dest: int) -> None:
        """Send a poison pill to a worker."""
        self.comm.send(None, dest=dest, tag=TAG_STOP)

    def recv_any(self) -> Any:
        """Block until a message arrives from any source / any tag."""
        status = self.comm.Probe()
        return self.comm.recv(source=status.Get_source(), tag=status.Get_tag())

    # ------------------------------------------------------------------
    # Shutdown
    # ------------------------------------------------------------------

    def shutdown_workers(self) -> None:
        """Send poison pills to all workers."""
        for rank in range(1, self.size):
            self.send_stop(rank)

    # ------------------------------------------------------------------
    # Barrier
    # ------------------------------------------------------------------

    def barrier(self) -> None:
        """Synchronise all ranks."""
        self.comm.Barrier()


# Import at bottom to avoid circular dependency
from polyscript.utils.mpi.splitter import WorkSplitter
