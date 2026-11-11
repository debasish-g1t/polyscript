#!/usr/bin/env python3
"""
Adapt a parquet file containing a ``polym`` column with alternative
representations (PSMILES, SMILES, SELFIES, BigSMILES) using the adapters.

Works on any file with a ``polym`` column — generator outputs, reference
training data, etc.

Usage
-----
    python adapt_dataframe.py <input.parquet> [--output OUTPUT] [--data-dir DIR]
                              [--progress-batch N] [--log-level LEVEL]

Examples
--------
    # Reference training data (fast)
    python adapt_dataframe.py ref_cryst_tendency_train.parquet
    python adapt_dataframe.py ref_IP_train.parquet

    # Merged generator data (slow — 2.9M rows with SELFIES encoding)
    python adapt_dataframe.py merged_polymer_data_filtered.parquet

    # Custom output, silent mode
    python adapt_dataframe.py input.parquet --output out.parquet --log-level 0
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from polyscript.utils.adapters import enrich
from polyscript.utils.logger import PolyLogger


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Adapt a dataframe with PSMILES, SMILES, SELFIES, BigSMILES columns.",
    )
    p.add_argument(
        "input",
        type=Path,
        help="Input parquet file (must contain a 'polym' column).",
    )
    p.add_argument(
        "--output", "-o",
        type=Path,
        default=None,
        help="Output parquet path (default: <input>_enriched.parquet in --data-dir).",
    )
    p.add_argument(
        "--data-dir", "-d",
        type=Path,
        default=Path("."),
        help="Default output directory when --output is not given (default: cwd).",
    )
    p.add_argument(
        "--progress-batch", "-b",
        type=int,
        default=5_000,
        help="Log progress every N rows (default: 5000).",
    )
    p.add_argument(
        "--log-level", "-l",
        type=str,
        default="INFO",
        help="Log level: SILENT, WARNING, INFO, DEBUG (default: INFO).",
    )
    return p


def main() -> None:
    args = _build_parser().parse_args()

    # Resolve output path
    if args.output is not None:
        output_path = args.output
    else:
        stem = args.input.stem
        output_path = args.data_dir / f"{stem}_enriched.parquet"

    logger = PolyLogger(name="adapt_dataframe", level=args.log_level)
    logger.add_console(level=args.log_level)

    logger.info("=" * 60)
    logger.info("  Enrich: %s", args.input.name)
    logger.info("  Output: %s", output_path.name)
    logger.info("=" * 60)

    logger.info("Loading %s …", args.input)
    df = pd.read_parquet(args.input)
    logger.info("Loaded %d rows, columns: %s", len(df), list(df.columns))

    logger.info("Converting representations …")
    df = enrich(df, progress_batch=args.progress_batch, logger=logger)

    logger.info("Writing → %s", output_path)
    df.to_parquet(output_path, index=False)
    logger.info("Done: %s  (%d rows, %s)", output_path, len(df), list(df.columns))
    logger.flush()


if __name__ == "__main__":
    main()
