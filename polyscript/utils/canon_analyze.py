"""
Canonicalization duplicate analyser.

``CanonicalizationAnalyzer`` scans parquet files under a base directory,
canonicalizes PSMILES, and reports duplication statistics.
"""

from __future__ import annotations

import sys
from collections import defaultdict
from pathlib import Path
from typing import List, Optional

import pandas as pd
from canonicalize_psmiles.canonicalize import canonicalize

from polyscript.utils.logger import INFO, PolyLogger, get_logger
from polyscript.utils.parsers import PolyScriptParser

DEFAULT_CLASSES = [
    "polyamide", "polyester", "polyether", "polyimide",
    "polyolefin", "polyoxazolidone", "polyurethane",
]


class CanonicalizationAnalyzer:
    """Analyse canonicalized polymer output for remaining duplicates.

    Parameters
    ----------
    base_dir : str
        Base directory with per-class subdirectories of parquet files.
    classes : list[str], optional
        Polymer classes to analyse (auto-detects all subdirs when ``None``).
    max_rows : int
        Max rows to collect per class (default -1 to include all the rows).
        This arg helps limit the rows for analysis in case of large datasets.
    logger : PolyLogger, optional
    """

    def __init__(
        self,
        base_dir: str,
        classes: Optional[List[str]] = None,
        max_rows: int = -1,
        logger: Optional[PolyLogger] = None,
    ):
        self._base = Path(base_dir)
        if not self._base.is_dir():
            raise FileNotFoundError(f"Base directory not found: {base_dir}")

        self._classes = classes
        self._max_rows = max_rows
        self.logger = logger or get_logger("canon_analyze", level=INFO)

    def collect(self) -> pd.DataFrame:
        """Scan parquet files and return a reference DataFrame."""
        if self._classes:
            classes = self._classes
        else:
            classes = sorted(
                d.name for d in self._base.iterdir()
                if d.is_dir() and not d.name.startswith(".")
            )
        if not classes:
            raise FileNotFoundError(f"No class subdirectories under {self._base}")

        self.logger.info("Base: %s  Classes: %s  Max rows/class: %d",
                         self._base, classes, self._max_rows)

        parser = PolyScriptParser()
        all_rows: list[dict] = []

        for pclass in classes:
            class_dir = self._base / pclass
            if not class_dir.is_dir():
                self.logger.warning("Skipping missing: %s", class_dir)
                continue

            collected = 0
            for pq in sorted(class_dir.glob("*.parquet")):
                if collected >= self._max_rows and self._max_rows > 0:
                    break
                try:
                    df = pd.read_parquet(pq)
                except Exception:
                    self.logger.debug("Failed to read %s", pq)
                    continue

                for _, row in df.iterrows():
                    if collected >= self._max_rows and self._max_rows > 0:
                        break
                    ps = str(row["polym"])
                    if "|||" in ps:
                        continue
                    errs = parser.parse(ps)
                    if errs:
                        continue
                    mon_pair = tuple(sorted(parser.monomers))
                    psmiles = parser.polymer
                    try:
                        canon = canonicalize(psmiles)
                    except Exception:
                        canon = f"ERR:{psmiles}"
                    all_rows.append({
                        "polymer_class": pclass,
                        "monomer_pair": str(mon_pair),
                        "original_psmiles": psmiles,
                        "canonicalized_psmiles": canon,
                        "source_file": pq.name,
                    })
                    collected += 1

            self.logger.info("%s: %d rows collected", pclass, collected)

        ref = pd.DataFrame(all_rows)
        self.logger.info("Total rows: %d", len(ref))
        return ref

    def analyse(self, ref: pd.DataFrame) -> dict:
        """Run duplication analysis on *ref* and return summary dict."""
        if ref.empty:
            self.logger.warning("No data — nothing to analyse.")
            return {"unique": 0, "total": 0, "dup_groups": 0}

        groups: dict[str, list] = defaultdict(list)
        for _, r in ref.iterrows():
            groups[r["canonicalized_psmiles"]].append(r)

        total_unique = len(groups)
        dups = {k: v for k, v in groups.items() if len(v) > 1}
        same_mon = sum(
            1 for rows in dups.values()
            if len({r["monomer_pair"] for r in rows}) == 1
        )
        diff_mon = len(dups) - same_mon

        self.logger.info("Unique canonical pSMILES: %d", total_unique)
        self.logger.info("Total rows: %d", len(ref))
        self.logger.info("Duplication: %d dups (%.1f%%)",
                         len(ref) - total_unique,
                         (1 - total_unique / len(ref)) * 100)
        self.logger.info("Duplicate groups: %d", len(dups))
        self.logger.info("  Same monomer pair:     %d", same_mon)
        self.logger.info("  Different monomer pair: %d", diff_mon)

        # Per-class
        self.logger.info("--- Per class ---")
        classes = ref["polymer_class"].unique()
        for pclass in sorted(classes):
            sub = ref[ref["polymer_class"] == pclass]
            n = len(sub)
            nu = sub["canonicalized_psmiles"].nunique()
            self.logger.info("  %s: %d rows -> %d unique (%d dups)",
                             pclass, n, nu, n - nu)

        # Examples
        self._show_examples(dups, same_mon, diff_mon)

        return {
            "unique": total_unique,
            "total": len(ref),
            "dup_groups": len(dups),
            "same_monomer": same_mon,
            "diff_monomer": diff_mon,
        }

    def _show_examples(self, dups: dict, same_mon: int, diff_mon: int):
        if same_mon > 0:
            self.logger.info("--- Same-monomer-pair duplicates ---")
            shown = 0
            for canon, rows in dups.items():
                if len({r["monomer_pair"] for r in rows}) == 1:
                    self.logger.info("  Canonical: %s", canon)
                    for r in rows:
                        self.logger.info("    %s", r["original_psmiles"])
                    shown += 1
                    if shown >= 5:
                        break

        if diff_mon > 0:
            self.logger.info("--- Different-monomer-pair duplicates ---")
            shown = 0
            for canon, rows in dups.items():
                mpairs = {r["monomer_pair"] for r in rows}
                if len(mpairs) > 1:
                    self.logger.info("  Canonical: %s", canon)
                    for mp in mpairs:
                        self.logger.info("    monomers: %s", mp)
                    shown += 1
                    if shown >= 5:
                        break
