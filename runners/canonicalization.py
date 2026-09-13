#!/usr/bin/env python3
"""
CLI entry point for canonicalization duplicate analysis.

Delegates to ``CanonicalizationAnalyzer`` (``polyscript.utils.canon_analyze``).

Usage:
    python canonicalization.py <BASE_DIR> [--classes C ...] [--max-rows N] [-o out.parquet]
"""

from __future__ import annotations

import argparse
from pathlib import Path

from polyscript.utils.canon_analyze import CanonicalizationAnalyzer
from polyscript.utils.logger import PolyLogger


def main() -> None:
    p = argparse.ArgumentParser(description="Analyse canonicalized polymer duplicates.")
    p.add_argument("base_dir", type=Path)
    p.add_argument("--classes", "-c", nargs="*", default=None)
    p.add_argument("--max-rows", "-n", type=int, default=200)
    p.add_argument("--output", "-o", type=Path, default=None)
    p.add_argument("--log-level", "-l", default="INFO")
    args = p.parse_args()

    analyzer = CanonicalizationAnalyzer(
        base_dir=str(args.base_dir),
        classes=args.classes,
        max_rows=args.max_rows,
    )
    analyzer.logger.set_level(args.log_level)

    ref = analyzer.collect()
    if args.output and not ref.empty:
        ref.to_parquet(args.output, index=False)
        print(f"Saved {len(ref)} rows to {args.output}")

    analyzer.analyse(ref)
    analyzer.logger.flush()


if __name__ == "__main__":
    main()
