"""Polymerizer package.

Provides the :class:`Polymerizer` engine, MPI-parallel runners, and
polymerisation utility functions.
"""

from polyscript.polymerizer.core import Polymerizer, get_rules_dir
from polyscript.polymerizer.mpi import MpiPolymerizer, SKIP_CLS
from polyscript.polymerizer.mpi_df import MpiRunnerDf, run_from_dataframes
from polyscript.polymerizer.utils import *  # noqa: F403

__all__ = [
    # core
    "Polymerizer",
    "get_rules_dir",
    # mpi
    "MpiPolymerizer",
    "SKIP_CLS",
    # mpi_df
    "MpiRunnerDf",
    "run_from_dataframes",
]
