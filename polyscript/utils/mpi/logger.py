"""
MPI-aware logger built on ``PolyLogger``.

Handles the quirks of MPI environments: direct-stderr for unbuffered
output, optional Rich progress bars, and a ``tail -f`` friendly
progress file.
"""

from __future__ import annotations

import logging
import os
from typing import Any, Optional

from polyscript.utils.logger import INFO, PolyLogger


# ---------------------------------------------------------------------------
# Direct-stderr handler (bypasses mpirun buffering)
# ---------------------------------------------------------------------------


class _DirectStderrHandler(logging.Handler):
    """Write log records directly to fd 2 via ``os.write``."""

    def emit(self, record: logging.LogRecord) -> None:
        try:
            msg = self.format(record) + "\n"
            os.write(2, msg.encode("utf-8"))
        except Exception:
            self.handleError(record)


# ---------------------------------------------------------------------------
# MpiLogger
# ---------------------------------------------------------------------------


class MpiLogger(PolyLogger):
    """``PolyLogger`` pre-configured for MPI environments.

    Adds a direct-stderr handler (mpirun-safe), an optional progress
    file (``tail -f`` friendly), and an optional Rich console handler
    on top of the standard file handler.

    Parameters
    ----------
    name : str
        Logger name.
    level : int or str
        Master logging threshold.
    log_filepath : str, optional
        Path to the detailed log file.
    progress_file : str, optional
        Path to a concise progress file (overwritten on each run).
    use_rich : bool
        Attach a ``RichHandler`` for progress-bar output.
    use_direct_stderr : bool
        Attach the unbuffered stderr handler.
    """

    def __init__(
        self,
        name: str,
        level: int | str = INFO,
        *,
        log_filepath: Optional[str] = None,
        progress_file: Optional[str] = None,
        use_rich: bool = True,
        use_direct_stderr: bool = True,
    ) -> None:
        super().__init__(name=name, level=level)
        # Clear the auto-added console handler from PolyLogger
        self._logger.handlers.clear()

        # ── Direct stderr (mpirun-safe) ───────────────────────────────
        if use_direct_stderr:
            h = _DirectStderrHandler()
            h.setLevel(logging.INFO)
            h.setFormatter(
                logging.Formatter(
                    "%(asctime)s [%(levelname).1s] %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S",
                )
            )
            self._logger.addHandler(h)

        # ── Progress file ─────────────────────────────────────────────
        if progress_file is not None:
            os.makedirs(os.path.dirname(progress_file) or ".", exist_ok=True)
            h = logging.FileHandler(progress_file, mode="w")
            h.setLevel(logging.INFO)
            h.setFormatter(
                logging.Formatter(
                    "%(asctime)s [%(levelname).1s] %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S",
                )
            )
            self._logger.addHandler(h)

        # ── Detailed log file ─────────────────────────────────────────
        if log_filepath is not None:
            os.makedirs(os.path.dirname(log_filepath) or ".", exist_ok=True)
            h = logging.FileHandler(log_filepath)
            h.setLevel(logging.INFO)
            h.setFormatter(
                logging.Formatter(
                    "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
                )
            )
            self._logger.addHandler(h)

        # ── Rich console ──────────────────────────────────────────────
        self._use_rich = use_rich
        self._console: Any = None
        if use_rich:
            try:
                from rich.logging import RichHandler

                if self._console is None:
                    from rich.console import Console

                    self._console = Console(
                        stderr=True,
                        force_terminal=True,
                        force_interactive=True,
                    )
                rh = RichHandler(
                    console=self._console,
                    show_time=True,
                    show_path=False,
                    markup=False,
                )
                rh.setLevel(logging.INFO)
                rh.setFormatter(logging.Formatter("%(message)s"))
                self._logger.addHandler(rh)
            except ImportError:
                self._use_rich = False
