"""
Base executor class that provides an auto-configured ``PolyLogger`` to
every subclass, so individual modules don't need to declare or set up a
logger in their ``__init__``.

Usage
-----
    from polyscript.utils.executor import BaseExecutor

    class MyRunner(BaseExecutor):
        _log_level = "DEBUG"   # override the class default (optional)

        def run(self):
            self.logger.info("running ...")   # self.logger is already wired
"""

from __future__ import annotations

from polyscript.utils.logger import INFO, PolyLogger


class BaseExecutor:
    """Mixin / base class that wires up a ``PolyLogger`` automatically.

    Subclasses get ``self.logger`` pre-configured at construction time.
    The log level is determined by the **class attribute** ``_log_level``
    (default ``INFO``).  Override it on a subclass or call
    :meth:`set_log_level` to change the level for all future instances.

    If an explicit *logger* is passed to ``__init__``, it is used
    as-is and ``_log_level`` has no effect.

    Parameters
    ----------
    logger : PolyLogger, optional
        Pre-configured logger.  When ``None``, one is created
        automatically using ``_log_level``.
    """

    _log_level: int | str = INFO

    def __init__(self, *, logger: PolyLogger | None = None, **kwargs: object) -> None:
        super().__init__(**kwargs)
        if logger is not None:
            self.logger = logger
        else:
            # Avoid clobbering if a parent class already set self.logger.
            if not hasattr(self, "logger"):
                self.logger = PolyLogger(
                    name=self.__class__.__name__.lower(),
                    level=self._log_level,
                )
                self.logger.add_console(level=self._log_level)

    # ------------------------------------------------------------------
    # Class-level level control
    # ------------------------------------------------------------------

    @classmethod
    def set_log_level(cls, level: int | str) -> None:
        """Set ``_log_level`` on this class.

        Affects all **future** instances that don't pass an explicit
        *logger* at construction time.
        """
        cls._log_level = level
