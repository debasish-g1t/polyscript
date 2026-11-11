"""
Lightweight progress tracker for MPI batch workloads.

Tracks elapsed time, computes rate / ETA, and emits periodic log lines
with optional heartbeat detection for stalled runs.
"""

from __future__ import annotations

import time
from typing import Optional


class ProgressTracker:
    """Track and report progress of a batch workload.

    Parameters
    ----------
    total : int
        Total number of work units.
    log_interval : int
        Minimum work units between progress log lines.
    stall_timeout : float
        Seconds without progress before emitting a heartbeat line
        (default 300 = 5 min).
    """

    def __init__(
        self,
        total: int,
        log_interval: int = 100,
        stall_timeout: float = 300.0,
    ) -> None:
        self.total = total
        self.log_interval = max(log_interval, 1)
        self.stall_timeout = stall_timeout

        self.completed: int = 0
        self.failed: int = 0
        self.rows: int = 0  # optional row count for row-rate display

        self._t_start = time.monotonic()
        self._next_log_at = log_interval
        self._last_log_time = self._t_start

    # ------------------------------------------------------------------
    # Updates
    # ------------------------------------------------------------------

    def tick(self, ok: bool = True, rows: int = 0) -> None:
        """Record one completed (or failed) work unit."""
        self.completed += 1
        if not ok:
            self.failed += 1
        self.rows += rows

    @property
    def processed(self) -> int:
        """Total work units processed (completed + failed)."""
        return self.completed

    @property
    def elapsed(self) -> float:
        """Wall-clock seconds since the tracker was created."""
        return time.monotonic() - self._t_start

    @property
    def rate(self) -> float:
        """Work units per second."""
        e = self.elapsed
        return self.processed / e if e > 0 else 0.0

    @property
    def eta(self) -> float:
        """Estimated seconds remaining."""
        r = self.rate
        remaining = self.total - self.processed
        return remaining / r if r > 0 else -1.0

    @property
    def pct(self) -> float:
        """Percentage complete."""
        return (self.processed / self.total * 100) if self.total > 0 else 0.0

    # ------------------------------------------------------------------
    # Log-line helpers
    # ------------------------------------------------------------------

    def should_log(self) -> bool:
        """Return ``True`` if it's time to emit a progress line."""
        return self.processed >= self._next_log_at

    def mark_logged(self) -> None:
        """Advance the next-log threshold and update last-log timestamp."""
        self._next_log_at = self.processed + self.log_interval
        self._last_log_time = time.monotonic()

    def is_stalled(self) -> bool:
        """Return ``True`` if no progress has been logged for too long."""
        return (time.monotonic() - self._last_log_time) > self.stall_timeout

    # ------------------------------------------------------------------
    # Formatting
    # ------------------------------------------------------------------

    @staticmethod
    def _fmt_duration(seconds: float) -> str:
        if seconds < 0:
            return "--"
        m, s = divmod(int(seconds), 60)
        h, m = divmod(m, 60)
        if h > 0:
            return f"{h}h{m:02d}m"
        if m > 0:
            return f"{m}m{s:02d}s"
        return f"{s}s"

    def format_line(self, label: str = "") -> str:
        """Build a single-line progress string."""
        parts = [
            f"[{self.processed}/{self.total} {self.pct:.1f}%]",
            f"ok={self.completed}",
        ]
        if self.failed:
            parts.append(f"fail={self.failed}")
        if self.rows and self.elapsed > 0:
            row_rate = self.rows / self.elapsed
            parts.append(f"rows={self.rows:,}")
            parts.append(f"{row_rate:,.0f} r/s")
        parts.append(f"elapsed={self._fmt_duration(self.elapsed)}")
        parts.append(f"eta={self._fmt_duration(self.eta)}")
        if label:
            parts.append(label)
        return " | ".join(parts)
