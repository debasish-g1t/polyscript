"""
MPI-parallel runtime components.

Each component is independently usable — compose them into your runner
instead of inheriting from a monolithic base class.

Components
----------
* ``MpiLogger``        — Rich / direct-stderr / progress-file logging
* ``WorkSplitter``     — round-robin & chunk-based work distribution
* ``MpiProtocol``      — master↔worker scatter / gather / poison-pill
* ``ProgressTracker``  — elapsed time, rate, ETA, heartbeat detection
* ``ResourceMonitor``  — background RSS-memory sampler (psutil)
* ``CheckpointDB``     — SQLite task-completion tracker with JSON migration
* ``ProgressFile``     — atomic JSON progress snapshot for monitoring
"""

from polyscript.utils.mpi.checkpoint import CheckpointDB, ProgressFile
from polyscript.utils.mpi.logger import MpiLogger
from polyscript.utils.mpi.progress import ProgressTracker
from polyscript.utils.mpi.protocol import MpiProtocol
from polyscript.utils.mpi.resource import ResourceMonitor
from polyscript.utils.mpi.splitter import WorkSplitter

__all__ = [
    "CheckpointDB",
    "MpiLogger",
    "MpiProtocol",
    "ProgressFile",
    "ProgressTracker",
    "ResourceMonitor",
    "WorkSplitter",
]
