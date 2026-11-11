"""
Flexible multi-level logger with pluggable output sinks.

Supports logging at standard severity levels (DEBUG, INFO, WARNING, ERROR,
CRITICAL) and allows different "versions" of the output by attaching
independently configured handlers.  Each handler can have its own level
threshold and format, so a single log call can simultaneously produce:

- concise timestamps on the console,
- rich structured output in a file,
- machine-readable JSON lines to a separate stream.

Usage
-----
    from polyscript.utils.logger import PolyLogger

    logger = PolyLogger(name="my-module")
    logger.add_console(level="INFO")
    logger.add_file("app.log", level="DEBUG")
    logger.add_stream(sys.stderr, level="WARNING")

    logger.info("Starting polymerisation")          # console + file
    logger.debug("Monomer batch: ABC-123")           # file only
    logger.warning("High temperature detected!")     # all three sinks
"""

from __future__ import annotations

import logging
import sys
from typing import Optional, TextIO

# ---------------------------------------------------------------------------
# Public re-exports so callers can import level names from this module.
# ---------------------------------------------------------------------------
DEBUG = logging.DEBUG
INFO = logging.INFO
WARNING = logging.WARNING
ERROR = logging.ERROR
CRITICAL = logging.CRITICAL

# Above CRITICAL so that nothing passes through when set as the threshold.
SILENT = 100

LEVEL_NAME_MAP: dict[str, int] = {
    "DEBUG": DEBUG,
    "INFO": INFO,
    "WARNING": WARNING,
    "ERROR": ERROR,
    "CRITICAL": CRITICAL,
    "SILENT": SILENT,
}

# Convenience mapping: small integers → log levels.
#   0 = SILENT, 1 = WARNING, 2 = INFO, 3 = DEBUG.
# Integers outside this range are passed through as literal levels.
VERBOSITY_TO_LEVEL: dict[int, int] = {
    0: SILENT,
    1: WARNING,
    2: INFO,
    3: DEBUG,
}

# Default formats for convenience ------------------------------------------
FMT_CONSOLE = "%(asctime)s [%(levelname).1s] %(message)s"
FMT_FILE = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
FMT_DETAILED = (
    "%(asctime)s | %(name)s | %(levelname)s | %(filename)s:%(lineno)d | %(message)s"
)
FMT_JSON = (
    '{"ts":"%(asctime)s","lvl":"%(levelname)s","name":"%(name)s","msg":"%(message)s"}'
)

DATE_FMT = "%Y-%m-%d %H:%M:%S"


# ---------------------------------------------------------------------------
# Core logger
# ---------------------------------------------------------------------------


class PolyLogger:
    """A logger that routes messages to multiple outputs with per-output levels.

    Each output (console, file, stream) is backed by a standard-library
    ``logging.Handler`` and can be configured with its own *threshold* and
    *format string*.  The logger itself also carries a master threshold;
    messages below that threshold are discarded immediately.

    Parameters
    ----------
    name : str
        Logger name (typically ``__name__`` of the calling module).
    level : int or str, optional
        Master logging threshold.  Defaults to ``DEBUG``.
    propagate : bool
        Whether to propagate messages to ancestor loggers (default ``False``).
    """

    def __init__(
        self,
        name: str,
        level: int | str = DEBUG,
        propagate: bool = False,
    ) -> None:
        self._logger = logging.getLogger(name)
        self._logger.setLevel(self._normalise_level(level))
        self._logger.propagate = propagate

        # Avoid cascading handlers when re-creating a logger.
        self._logger.handlers.clear()

    # ------------------------------------------------------------------
    # Public helpers
    # ------------------------------------------------------------------

    @staticmethod
    def level_name(level: int | str) -> str:
        """Return the canonical level name for *level*."""
        return logging.getLevelName(PolyLogger._normalise_level(level))

    @staticmethod
    def _normalise_level(level: int | str) -> int:
        if isinstance(level, str):
            return LEVEL_NAME_MAP.get(level.upper(), DEBUG)
        # Resolve small verbosity integers (0-3) to named levels;
        # everything else passes through as a literal logging level.
        return VERBOSITY_TO_LEVEL.get(level, level)

    # ------------------------------------------------------------------
    # Handler factories
    # ------------------------------------------------------------------

    def add_console(
        self,
        level: int | str = INFO,
        fmt: str = FMT_CONSOLE,
        datefmt: str = DATE_FMT,
        stream: TextIO | None = None,
    ) -> PolyLogger:
        """Add a handler that prints to stderr (or *stream*).

        Returns *self* so calls can be chained.
        """
        handler = logging.StreamHandler(stream if stream is not None else sys.stderr)
        handler.setLevel(self._normalise_level(level))
        handler.setFormatter(logging.Formatter(fmt, datefmt=datefmt))
        self._logger.addHandler(handler)
        return self

    def add_file(
        self,
        path: str,
        level: int | str = DEBUG,
        fmt: str = FMT_FILE,
        datefmt: str = DATE_FMT,
        mode: str = "a",
    ) -> PolyLogger:
        """Add a handler that appends to *path*.

        Returns *self* so calls can be chained.
        """
        handler = logging.FileHandler(path, mode=mode)
        handler.setLevel(self._normalise_level(level))
        handler.setFormatter(logging.Formatter(fmt, datefmt=datefmt))
        self._logger.addHandler(handler)
        return self

    def add_stream(
        self,
        stream: TextIO,
        level: int | str = INFO,
        fmt: str = FMT_CONSOLE,
        datefmt: str = DATE_FMT,
    ) -> PolyLogger:
        """Add a handler that writes to an arbitrary *stream*.

        Returns *self* so calls can be chained.
        """
        handler = logging.StreamHandler(stream)
        handler.setLevel(self._normalise_level(level))
        handler.setFormatter(logging.Formatter(fmt, datefmt=datefmt))
        self._logger.addHandler(handler)
        return self

    def add_handler(
        self,
        handler: logging.Handler,
        level: Optional[int | str] = None,
    ) -> PolyLogger:
        """Attach a pre-configured :class:`logging.Handler`.

        If *level* is given the handler's threshold is updated.

        Returns *self* so calls can be chained.
        """
        if level is not None:
            handler.setLevel(self._normalise_level(level))
        self._logger.addHandler(handler)
        return self

    # ------------------------------------------------------------------
    # Level control
    # ------------------------------------------------------------------

    def set_level(self, level: int | str) -> None:
        """Set the master logging threshold."""
        self._logger.setLevel(self._normalise_level(level))

    def set_verbosity(self, verbosity: int) -> None:
        """Set the master threshold from a small verbosity integer.

        ``0`` = silent, ``1`` = warnings, ``2`` = info, ``3`` = debug.
        """
        self._logger.setLevel(self._normalise_level(verbosity))

    @property
    def level(self) -> int:
        return self._logger.level

    @property
    def is_debug(self) -> bool:
        return self._logger.isEnabledFor(DEBUG)

    @property
    def is_info(self) -> bool:
        return self._logger.isEnabledFor(INFO)

    @property
    def is_warning(self) -> bool:
        return self._logger.isEnabledFor(WARNING)

    @property
    def is_error(self) -> bool:
        return self._logger.isEnabledFor(ERROR)

    # ------------------------------------------------------------------
    # Logging methods
    # ------------------------------------------------------------------

    def debug(self, msg: str, *args, **kwargs) -> None:
        self._logger.debug(msg, *args, **kwargs)

    def info(self, msg: str, *args, **kwargs) -> None:
        self._logger.info(msg, *args, **kwargs)

    def warning(self, msg: str, *args, **kwargs) -> None:
        self._logger.warning(msg, *args, **kwargs)

    def error(self, msg: str, *args, **kwargs) -> None:
        self._logger.error(msg, *args, **kwargs)

    def critical(self, msg: str, *args, **kwargs) -> None:
        self._logger.critical(msg, *args, **kwargs)

    # Convenience aliases
    warn = warning

    # ------------------------------------------------------------------
    # Compatibility with the standard library
    # ------------------------------------------------------------------

    @property
    def std_logger(self) -> logging.Logger:
        """Expose the underlying ``logging.Logger``.

        Useful when a function expects a standard-library logger (e.g.
        third-party libraries that accept a ``logger=`` kwarg).
        """
        return self._logger

    # ------------------------------------------------------------------
    # Clean-up
    # ------------------------------------------------------------------

    def flush(self) -> None:
        """Flush all handlers without closing them.

        Useful in HPC / batch contexts where you want to force buffered
        log output to disk periodically.
        """
        for handler in self._logger.handlers:
            handler.flush()

    def close(self) -> None:
        """Flush and close all handlers."""
        for handler in tuple(self._logger.handlers):
            handler.flush()
            handler.close()
            self._logger.removeHandler(handler)

    def __repr__(self) -> str:
        handlers = ", ".join(type(h).__name__ for h in self._logger.handlers)
        return (
            f"PolyLogger(name={self._logger.name!r}, "
            f"level={self.level_name(self.level)}, "
            f"handlers=[{handlers}])"
        )


# ---------------------------------------------------------------------------
# Quick-start convenience
# ---------------------------------------------------------------------------


def get_logger(
    name: str,
    level: int | str = INFO,
    *,
    log_file: Optional[str] = None,
    console: bool = True,
) -> PolyLogger:
    """Create a pre-configured ``PolyLogger`` in one call.

    Parameters
    ----------
    name : str
        Logger name.
    level : int or str
        Master threshold.
    log_file : str, optional
        If provided, a file handler is added at this path.
    console : bool
        If ``True`` (the default), a console (stderr) handler is added.
    """
    logger = PolyLogger(name=name, level=level)
    if console:
        logger.add_console(level=level)
    if log_file is not None:
        logger.add_file(log_file, level=DEBUG)
    return logger
