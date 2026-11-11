"""
SQLite-backed task checkpoint for MPI batch workloads.

Tracks which (p_class, mon_pair_key, batch_key) tasks have been
completed, supports legacy JSON migration, and provides lightweight
query helpers for progress tracking.

Also includes ``ProgressFile`` — an atomic JSON progress snapshot
for external monitoring tools.
"""

from __future__ import annotations

import json
import os
import sqlite3
import time
from typing import Any, Dict, List, Optional, Tuple


class ProgressFile:
    """Write an atomic progress snapshot as JSON.

    Parameters
    ----------
    path : str
        Full path to the progress JSON file (e.g. ``progress.json``).
    """

    def __init__(self, path: str) -> None:
        self._path = path

    def write(self, completed: int, failed: int, total: int) -> None:
        """Atomically write the current progress snapshot."""
        data: Dict[str, Any] = {
            "completed": completed,
            "failed": failed,
            "total": total,
            "pending": total - completed - failed,
            "progress_pct": (
                round((completed + failed) / total * 100, 2) if total > 0 else 0
            ),
            "success_rate": (
                round(completed / (completed + failed) * 100, 2)
                if (completed + failed) > 0
                else 0
            ),
        }
        tmp = self._path + ".tmp"
        with open(tmp, "w") as fh:
            json.dump(data, fh, indent=2)
        try:
            for attempt in range(3):
                try:
                    os.replace(tmp, self._path)
                    return
                except PermissionError:
                    time.sleep(0.05 * (attempt + 1))
            os.replace(tmp, self._path)
        except PermissionError:
            try:
                with open(self._path, "w") as fh:
                    json.dump(data, fh, indent=2)
            except Exception:
                pass
            try:
                os.remove(tmp)
            except OSError:
                pass


# ---------------------------------------------------------------------------
# CheckpointDB
# ---------------------------------------------------------------------------


class CheckpointDB:
    """SQLite checkpoint database for MPI task tracking.

    Schema::

        CREATE TABLE checkpoints (
            p_class      TEXT NOT NULL,
            mon_pair_key TEXT NOT NULL,
            mon_type1    TEXT NOT NULL,
            mon_type2    TEXT NOT NULL,
            batch_key    TEXT NOT NULL,
            completed    INTEGER NOT NULL DEFAULT 0,
            PRIMARY KEY (p_class, mon_pair_key, batch_key)
        )

    Parameters
    ----------
    db_path : str
        Path to the SQLite database file.
    """

    def __init__(self, db_path: str) -> None:
        self._db_path = db_path
        self._conn: Optional[sqlite3.Connection] = None

    # ------------------------------------------------------------------
    # Connection management
    # ------------------------------------------------------------------

    @property
    def conn(self) -> sqlite3.Connection:
        """Lazy-connect (or return existing connection)."""
        if self._conn is None:
            self._conn = sqlite3.connect(self._db_path)
            self._conn.execute("PRAGMA journal_mode=WAL")
            self._conn.execute("PRAGMA synchronous=NORMAL")
        return self._conn

    def close(self) -> None:
        if self._conn is not None:
            self._conn.close()
            self._conn = None

    # ------------------------------------------------------------------
    # Schema
    # ------------------------------------------------------------------

    def init_schema(self) -> None:
        """Ensure the checkpoints table and index exist."""
        db = self.conn
        db.execute(
            """
            CREATE TABLE IF NOT EXISTS checkpoints (
                p_class      TEXT NOT NULL,
                mon_pair_key TEXT NOT NULL,
                mon_type1    TEXT NOT NULL,
                mon_type2    TEXT NOT NULL,
                batch_key    TEXT NOT NULL,
                completed    INTEGER NOT NULL DEFAULT 0,
                PRIMARY KEY (p_class, mon_pair_key, batch_key)
            )
            """
        )
        db.execute(
            """
            CREATE INDEX IF NOT EXISTS idx_pending
            ON checkpoints(completed) WHERE completed = 0
            """
        )
        db.commit()

    # ------------------------------------------------------------------
    # Populate
    # ------------------------------------------------------------------

    def populate(
        self, rows: List[Tuple[str, str, str, str, str, int]]
    ) -> None:
        """Insert or ignore a batch of task rows.

        Each row is ``(p_class, mon_pair_key, mon_type1, mon_type2,
        batch_key, completed)``.
        """
        db = self.conn
        db.executemany(
            "INSERT OR IGNORE INTO checkpoints "
            "(p_class, mon_pair_key, mon_type1, mon_type2, batch_key, completed) "
            "VALUES (?, ?, ?, ?, ?, ?)",
            rows,
        )
        db.commit()

    def migrate_json(self, json_path: str) -> None:
        """Import a legacy ``ckpt.json`` file into the SQLite DB."""
        with open(json_path) as fh:
            ckpt = json.load(fh)

        rows: List[Tuple[str, str, str, str, str, int]] = []
        for p_class, mon_pair_list in ckpt.items():
            for mon_pair_dict in mon_pair_list:
                for mon_pair_key, batch_dict in mon_pair_dict.items():
                    parts = mon_pair_key.rsplit("_", 1)
                    if len(parts) == 2:
                        mon_type1, mon_type2 = parts
                    else:
                        mon_type1 = mon_pair_key
                        mon_type2 = "none"

                    for batch_key, is_complete in batch_dict.items():
                        rows.append(
                            (
                                p_class,
                                mon_pair_key,
                                mon_type1,
                                mon_type2,
                                batch_key,
                                1 if is_complete else 0,
                            )
                        )

        self.populate(rows)
        backup = json_path + ".bak"
        os.rename(json_path, backup)

    # ------------------------------------------------------------------
    # Queries
    # ------------------------------------------------------------------

    def count_total(self) -> int:
        return self.conn.execute(
            "SELECT COUNT(*) FROM checkpoints"
        ).fetchone()[0]

    def count_pending(self, p_classes: Optional[List[str]] = None) -> int:
        query = "SELECT COUNT(*) FROM checkpoints WHERE completed = 0"
        params: tuple = ()
        if p_classes:
            placeholders = ", ".join("?" for _ in p_classes)
            query += f" AND p_class IN ({placeholders})"
            params = tuple(p_classes)
        return self.conn.execute(query, params).fetchone()[0]

    def list_classes(self) -> List[str]:
        rows = self.conn.execute(
            "SELECT DISTINCT p_class FROM checkpoints ORDER BY p_class"
        ).fetchall()
        return [r[0] for r in rows]

    def mark_complete(
        self,
        p_class: str,
        mon_pair_key: str,
        batch_key: str,
        success: bool = True,
    ) -> None:
        """Mark (or un-mark) a single task as completed."""
        self.conn.execute(
            "UPDATE checkpoints SET completed = ? "
            "WHERE p_class = ? AND mon_pair_key = ? AND batch_key = ?",
            (1 if success else 0, p_class, mon_pair_key, batch_key),
        )
        self.conn.commit()

    def reset(self, p_class: Optional[str] = None) -> None:
        """Reset all tasks (or tasks for a specific polymer class)."""
        if p_class:
            self.conn.execute(
                "UPDATE checkpoints SET completed = 0 WHERE p_class = ?",
                (p_class,),
            )
        else:
            self.conn.execute("UPDATE checkpoints SET completed = 0")
        self.conn.commit()

    def lookup_mon_pair_key(
        self, p_class: str, batch_key: str
    ) -> Optional[str]:
        """Look up the ``mon_pair_key`` for a given class + batch key.

        Returns ``None`` if not found.
        """
        row = self.conn.execute(
            "SELECT mon_pair_key FROM checkpoints "
            "WHERE p_class = ? AND batch_key = ?",
            (p_class, batch_key),
        ).fetchone()
        return row[0] if row else None
