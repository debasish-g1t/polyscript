#!/usr/bin/env python3
"""
CLI entry point for candidate filtering.

Delegates to ``CandidateFilterRunner`` (``polyscript.utils.filter_runner``).

Usage:
    python filter_candidates.py --input file.parquet [-o out.parquet]
    python filter_candidates.py --dir ./outputs/run_0 [-o ./filtered]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from polyscript.utils.filter_runner import CandidateFilterRunner
from polyscript.utils.logger import PolyLogger


def main() -> None:
    p = argparse.ArgumentParser(description="Filter polymer candidates.")
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--input", type=Path)
    g.add_argument("--dir", type=Path)
    p.add_argument("--output", "-o", type=Path, default=None)
    p.add_argument("--output-dir", type=Path, default=None)
    p.add_argument("--log-level", "-l", default="INFO")
    args = p.parse_args()

    logger = PolyLogger(name="filter_candidates", level=args.log_level)
    logger.add_console(level=args.log_level)
    runner = CandidateFilterRunner(logger=logger)

    if args.input:
        try:
            runner.run_single(args.input, args.output)
        except FileNotFoundError as e:
            logger.error("%s", e)
            sys.exit(1)
    else:
        try:
            runner.run_dir(args.dir, args.output_dir)
        except FileNotFoundError as e:
            logger.error("%s", e)
            sys.exit(1)

    logger.info("Done.")
    logger.flush()


if __name__ == "__main__":
    main()
