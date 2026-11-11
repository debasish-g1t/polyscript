#!/usr/bin/env python3
"""Monomer classification runner for the BRICS dataset (Workflow step 1.4).

Reads the output of the BRICS pipeline (``1.2.0_decomp.csv`` fragments or
``1.3.0_all_rxn_data_BRICS.parquet`` compose output), classifies the SMILES
column with ``PolyScriptClassifier``, and writes the enriched table.

Usage
-----
    python 1_4_monomer_classification.py INPUT [-o OUTPUT] [-c COLUMN]
                                     [--min-fg N] [--max-fg N]
                                     [--no-carbonate] [-l LEVEL]

Examples
--------
    # Compose output; SMILES column is auto-detected (default column order:
    # product_smiles -> decomp_SMILES -> first column containing "smiles")
    python 1_4_monomer_classification.py 1.3.0_all_rxn_data_BRICS.parquet

    # Decomposed fragments
    python 1_4_monomer_classification.py 1.2.0_decomp.csv

    # Explicit column and output path
    python 1_4_monomer_classification.py 1.3.0_all_rxn_data_BRICS.parquet \
        -c product_smiles -o 1.4.0_classified.parquet
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

from polyscript.classifier import PolyScriptClassifier
from polyscript.utils.logger import PolyLogger

_READERS = {".csv": pd.read_csv, ".parquet": pd.read_parquet}
_WRITERS = {".csv": "to_csv", ".parquet": "to_parquet"}

# Preferred columns when --column is not given, in BRICS workflow order.
_PREFERRED_COLUMNS = ("product_smiles", "decomp_SMILES")


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Classify BRICS output monomers by functional-group type.",
    )
    p.add_argument(
        "input",
        type=Path,
        help="BRICS output file (.csv or .parquet).",
    )
    p.add_argument(
        "--output", "-o",
        type=Path,
        default=None,
        help="Output path (default: <input>_classified<same suffix>).",
    )
    p.add_argument(
        "--column", "-c",
        type=str,
        default=None,
        help="SMILES column to classify (default: auto-detect).",
    )
    p.add_argument(
        "--min-fg", type=int, default=2,
        help="Minimum functional groups for poly-functional classes (default: 2).",
    )
    p.add_argument(
        "--max-fg", type=int, default=4,
        help="Maximum functional groups for poly-functional classes (default: 4).",
    )
    p.add_argument(
        "--no-carbonate", action="store_true",
        help="Do not append the CO/HCHO reference rows used for carbonate classes.",
    )
    p.add_argument(
        "--log-level", "-l", type=str, default="INFO",
        help="Log level: SILENT, WARNING, INFO, DEBUG (default: INFO).",
    )
    return p


def _read_table(path: Path) -> pd.DataFrame:
    reader = _READERS.get(path.suffix.lower())
    if reader is None:
        raise ValueError(
            f"Unsupported input format '{path.suffix}' (expected .csv or .parquet)"
        )
    return reader(path)


def _resolve_column(df: pd.DataFrame, requested: str | None) -> str:
    if requested is not None:
        if requested not in df.columns:
            raise KeyError(
                f"Column '{requested}' not found (available: {list(df.columns)})"
            )
        return requested

    for candidate in _PREFERRED_COLUMNS:
        if candidate in df.columns:
            return candidate

    matches = [c for c in df.columns if "smiles" in c.lower()]
    if not matches:
        raise KeyError(
            f"No SMILES column found (available: {list(df.columns)})"
        )
    return matches[0]


def _write_table(df: pd.DataFrame, path: Path) -> None:
    writer = _WRITERS.get(path.suffix.lower())
    if writer is None:
        raise ValueError(
            f"Unsupported output format '{path.suffix}' (expected .csv or .parquet)"
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    getattr(df, writer)(path, index=False)


def main() -> None:
    args = _build_parser().parse_args()

    logger = PolyLogger(name="classify_monomers", level=args.log_level)
    logger.add_console(level=args.log_level)

    if not args.input.is_file():
        logger.error("Input file not found: %s", args.input)
        sys.exit(1)

    output_path = args.output
    if output_path is None:
        output_path = args.input.with_name(
            f"{args.input.stem}_classified{args.input.suffix}"
        )

    try:
        df = _read_table(args.input)
        column = _resolve_column(df, args.column)
    except (ValueError, KeyError) as exc:
        logger.error("%s", exc)
        sys.exit(1)

    logger.info("=" * 60)
    logger.info("  Classify: %s", args.input.name)
    logger.info("  Rows:     %d", len(df))
    logger.info("  Column:   %s", column)
    logger.info("  Output:   %s", output_path)
    logger.info("=" * 60)

    if not args.no_carbonate:
        logger.info("Carbonate mode: appending CO and HCHO reference rows")

    classifier = PolyScriptClassifier(
        MinFG=args.min_fg,
        MaxFG=args.max_fg,
        include_carbonate=not args.no_carbonate,
        logger=logger,
    )
    result = classifier.classify(df, col_name=column)

    added = [c for c in result.columns if c not in df.columns]
    class_cols = [c for c in added if pd.api.types.is_bool_dtype(result[c])]
    totals = result[class_cols].sum() if class_cols else {}
    detected = {str(name): int(n) for name, n in totals.items() if n}
    if detected:
        logger.info("Detected monomer types:")
        for name, n in sorted(detected.items(), key=lambda kv: kv[1], reverse=True):
            logger.info("  %-18s %d", name, n)

    try:
        _write_table(result, output_path)
    except ValueError as exc:
        logger.error("%s", exc)
        sys.exit(1)

    logger.info("Done: %s  (%d rows, %d class columns)",
                output_path, len(result), len(class_cols))
    logger.flush()


if __name__ == "__main__":
    main()
