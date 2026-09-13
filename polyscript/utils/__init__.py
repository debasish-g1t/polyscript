"""
Utilities package — logging, parsing, validation, adapters, and MPI components.

Exports
-------
Logger
    ``PolyLogger``, ``get_logger``, level constants (``INFO``, ``DEBUG``, …).

Executor
    ``BaseExecutor`` — auto-wired logger mixin class.

MPI
    ``BaseMPIExecutor`` — MPI-aware executor.
    ``MpiSeqValidator`` — distributed sequence validator.
    ``mpi`` subpackage — ``MpiLogger``, ``MpiProtocol``, ``CheckpointDB``, …

Parsing & Validation
    ``PolyScriptParser``, ``NestedPolyScriptParser``, ``SeqValidator``.

Adapters
    ``Adapter``, ``PSMILESAdapter``, ``SMILESAdapter``, ``SELFIESAdapter``,
    ``BigSMILESAdapter``.

Analysis
    ``CanonicalizationAnalyzer``.
"""

from polyscript.utils.adapters import (
    Adapter,
    BigSMILESAdapter,
    PSMILESAdapter,
    SELFIESAdapter,
    SMILESAdapter,
)
from polyscript.utils.canon_analyze import CanonicalizationAnalyzer
from polyscript.utils.executor import BaseExecutor
from polyscript.utils.logger import (
    CRITICAL,
    DEBUG,
    ERROR,
    INFO,
    SILENT,
    WARNING,
    PolyLogger,
    get_logger,
)
from polyscript.utils.mpi_executor import BaseMPIExecutor
from polyscript.utils.mpi_validator import MpiSeqValidator
from polyscript.utils.parsers import NestedPolyScriptParser, PolyScriptParser
from polyscript.utils.validators import SeqValidator

__all__ = [
    # logger
    "PolyLogger",
    "get_logger",
    "INFO",
    "DEBUG",
    "WARNING",
    "ERROR",
    "CRITICAL",
    "SILENT",
    # executor
    "BaseExecutor",
    # mpi
    "BaseMPIExecutor",
    "MpiSeqValidator",
    # parsers
    "PolyScriptParser",
    "NestedPolyScriptParser",
    # validators
    "SeqValidator",
    # adapters
    "Adapter",
    "PSMILESAdapter",
    "SMILESAdapter",
    "SELFIESAdapter",
    "BigSMILESAdapter",
    # analysis
    "CanonicalizationAnalyzer",
]
