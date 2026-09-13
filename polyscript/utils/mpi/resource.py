"""
Background resource (memory) monitor for MPI runs.

Samples RSS memory every 5 seconds in a daemon thread.  Requires
``psutil`` — falls back gracefully when it's not installed.
"""

from __future__ import annotations

import threading
import time
from typing import List, Optional, Tuple

try:
    import psutil

    _HAS_PSUTIL = True
except ImportError:
    _HAS_PSUTIL = False


class ResourceMonitor:
    """Sample RSS memory in a background daemon thread.

    Parameters
    ----------
    enabled : bool
        When ``False`` (or psutil unavailable), all methods are no-ops.
    """

    def __init__(self, enabled: bool = True) -> None:
        self._enabled = enabled and _HAS_PSUTIL
        self._samples: List[Tuple[float, float]] = []
        self._start_time: Optional[float] = None
        self._stop_event: Optional[threading.Event] = None
        self._thread: Optional[threading.Thread] = None

    # ------------------------------------------------------------------
    # Lifecycle
    # ------------------------------------------------------------------

    def start(self) -> None:
        if not self._enabled:
            return
        self._start_time = time.monotonic()
        self._samples = []
        self._stop_event = threading.Event()
        self._thread = threading.Thread(target=self._loop, daemon=True)
        self._thread.start()

    def stop(self) -> None:
        if self._stop_event is not None:
            self._stop_event.set()
        if self._thread is not None and self._thread.is_alive():
            self._thread.join(timeout=6.0)

    # ------------------------------------------------------------------
    # Internals
    # ------------------------------------------------------------------

    def _sample(self) -> None:
        if self._start_time is None:
            return
        try:
            proc = psutil.Process()  # type: ignore[possibly-undefined]
            mem = proc.memory_info()
            rss_mb = mem.rss / (1024 * 1024)
            elapsed = time.monotonic() - self._start_time
            self._samples.append((elapsed, rss_mb))
        except Exception:
            pass

    def _loop(self) -> None:
        while self._stop_event is not None and not self._stop_event.wait(5.0):
            self._sample()

    # ------------------------------------------------------------------
    # Query
    # ------------------------------------------------------------------

    @property
    def samples(self) -> List[Tuple[float, float]]:
        return list(self._samples)

    def summary(self) -> dict:
        if not self._samples:
            return {"peak_rss_mb": 0, "avg_rss_mb": 0, "wall_seconds": 0}
        rss = [v[1] for v in self._samples]
        wall = self._samples[-1][0] if self._samples else 0
        return {
            "peak_rss_mb": round(max(rss), 1),
            "avg_rss_mb": round(sum(rss) / len(rss), 1),
            "wall_seconds": round(wall, 1),
        }
