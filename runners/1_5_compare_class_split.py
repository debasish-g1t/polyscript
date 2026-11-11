#!/usr/bin/env python3
"""Comprehensive classification report + per-class dataset split (Workflow 1.5).

Merges the classification-report and class-split scripts into a single run:

1. Scan classified parquet files, count every monomer class and the number of
   "unknown" rows (rows where no class column is true), and write
   ``classification_report.{csv,png,txt}``.
2. For every class with at least one hit, collect the matching rows (keeping
   only the SMILES column and the class flag) and write them in batches::

       <output>/<class>_<batch_size>/<class>_<start_index>.parquet

Usage
-----
    python 1_5_compare_class_split.py INPUT [-o OUTPUT] [--report-dir DIR]
        [--batch-size N] [--smiles-column psmiles] [-l LEVEL]

INPUT may be a single classified parquet file or a directory of them (the
output of ``1_4_monomer_classification.py`` works either way).

Examples
--------
    python 1_5_compare_class_split.py outputs/classified/
    python 1_5_compare_class_split.py classified.parquet \
        -o outputs/cls_split_batch_5K
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd
import tqdm

from polyscript.classifier.class_reference import mon_type_index
from polyscript.utils.logger import PolyLogger

# Stable, de-duplicated class order from the classifier reference table.
CLASS_NAMES = list(dict.fromkeys(mon_type_index.values()))
UNKNOWN = "unknown"


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Report classification counts and split the dataset by class.",
    )
    p.add_argument(
        "input",
        type=Path,
        help="Classified parquet file, or a directory of parquet files.",
    )
    p.add_argument(
        "--output", "-o",
        type=Path,
        default=None,
        help="Root output directory (default: ./outputs/cls_split_batch_<N>).",
    )
    p.add_argument(
        "--report-dir",
        type=Path,
        default=None,
        help="Where to write the report files (default: the output directory).",
    )
    p.add_argument(
        "--batch-size",
        type=int,
        default=5000,
        help="Rows per output parquet chunk (default: 5000).",
    )
    p.add_argument(
        "--smiles-column",
        type=str,
        default="psmiles",
        help="SMILES column to keep in the split files (default: psmiles).",
    )
    p.add_argument(
        "--log-level", "-l",
        type=str,
        default="INFO",
        help="Log level: SILENT, WARNING, INFO, DEBUG (default: INFO).",
    )
    return p


def _list_parquet_files(path: Path) -> list[Path]:
    if path.is_file():
        if path.suffix.lower() != ".parquet":
            raise ValueError(f"Input file is not parquet: {path}")
        return [path]
    if path.is_dir():
        files = sorted(p for p in path.iterdir() if p.suffix.lower() == ".parquet")
        if not files:
            raise ValueError(f"No .parquet files found in {path}")
        return files
    raise FileNotFoundError(f"Input path does not exist: {path}")


def _build_report(
    files: list[Path], smiles_column: str, logger: PolyLogger
) -> dict[str, int]:
    """Count class hits and unknown rows across all classified batches."""
    totals = {name: 0 for name in CLASS_NAMES}
    totals[UNKNOWN] = 0
    smiles_found = False

    for path in tqdm.tqdm(files, desc="Scanning batches"):
        df = pd.read_parquet(path)
        if smiles_column in df.columns:
            smiles_found = True
        present = [c for c in CLASS_NAMES if c in df.columns]
        if present:
            sums = df[present].sum().to_dict()
            for name in present:
                totals[name] += int(sums[name])
            totals[UNKNOWN] += int((df[present].sum(axis=1) == 0).sum())
        else:
            totals[UNKNOWN] += len(df)
            logger.warning("%s has no monomer-class columns", path.name)

    if not smiles_found:
        raise ValueError(
            f"SMILES column '{smiles_column}' not found in any input file"
        )
    return totals


def _write_report(totals: dict[str, int], report_dir: Path, logger: PolyLogger) -> None:
    report_dir.mkdir(parents=True, exist_ok=True)

    pd.DataFrame(totals, index=[0]).to_csv(
        report_dir / "classification_report.csv", index=False
    )

    # Headless backend so the runner works on servers without a display.
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    names = list(totals)
    values = list(totals.values())
    fig, ax = plt.subplots(figsize=(max(10, len(names) * 0.55), 6))
    ax.bar(names, values)
    ax.set_xlabel("Classification Type")
    ax.set_ylabel("Number of Samples")
    ax.set_title("Classification Report")
    ax.tick_params(axis="x", rotation=70)
    fig.tight_layout()
    fig.savefig(report_dir / "classification_report.png", dpi=150)
    plt.close(fig)

    total = sum(values)
    with open(report_dir / "classification_report.txt", "w") as fh:
        fh.write("Classification Report\n")
        fh.write("----------------------------\n")
        for name, count in totals.items():
            pct = count / total * 100 if total else 0.0
            fh.write(f"{name}: {count} ({pct:.2f}%)\n")
        fh.write(f"Total: {total}\n")

    logger.info("Wrote classification report to %s", report_dir)


def _extract_class(
    type_name: str,
    files: list[Path],
    smiles_column: str,
    batch_size: int,
    output_dir: Path,
    logger: PolyLogger,
) -> int:
    """Collect rows where *type_name* is true and write them in batches."""
    type_dir = output_dir / f"{type_name}_{batch_size}"
    type_dir.mkdir(parents=True, exist_ok=True)

    frames = []
    for path in tqdm.tqdm(
        files, desc=f"Collecting {type_name}", leave=False
    ):
        df = pd.read_parquet(path)
        if type_name not in df.columns:
            continue
        if smiles_column not in df.columns:
            logger.warning(
                "Skipping %s: no '%s' column", path.name, smiles_column
            )
            continue
        mask = df[type_name].fillna(False).astype(bool)
        subset = df.loc[mask, [smiles_column, type_name]]
        if not subset.empty:
            frames.append(subset)

    if not frames:
        logger.warning("No rows found for class '%s'", type_name)
        return 0

    ref_df = pd.concat(frames, ignore_index=True)
    for start in range(0, len(ref_df), batch_size):
        ref_df.iloc[start:start + batch_size].to_parquet(
            type_dir / f"{type_name}_{start}.parquet", index=False
        )
    return len(ref_df)


def main() -> None:
    args = _build_parser().parse_args()

    output_dir = args.output
    if output_dir is None:
        output_dir = Path("./outputs") / f"cls_split_batch_{args.batch_size}"
    report_dir = args.report_dir or output_dir

    logger = PolyLogger(name="class_split", level=args.log_level)
    logger.add_console(level=args.log_level)

    try:
        files = _list_parquet_files(args.input)
    except (ValueError, FileNotFoundError) as exc:
        logger.error("%s", exc)
        sys.exit(1)

    logger.info("=" * 60)
    logger.info("  Input:   %d parquet file(s)", len(files))
    logger.info("  Output:  %s", output_dir)
    logger.info("  Batch:   %d rows", args.batch_size)
    logger.info("  SMILES:  %s", args.smiles_column)
    logger.info("=" * 60)

    try:
        totals = _build_report(files, args.smiles_column, logger)
    except ValueError as exc:
        logger.error("%s", exc)
        sys.exit(1)
    _write_report(totals, report_dir, logger)

    valid = [name for name in CLASS_NAMES if totals.get(name, 0)]
    logger.info("Extracting %d class(es): %s", len(valid), ", ".join(valid))
    output_dir.mkdir(parents=True, exist_ok=True)

    summary = {}
    for name in valid:
        summary[name] = _extract_class(
            name, files, args.smiles_column, args.batch_size, output_dir, logger
        )
        logger.info("  %-18s %d rows", name, summary[name])

    logger.info(
        "Done. Wrote %d class(es) to %s", len(summary), output_dir
    )
    logger.flush()


if __name__ == "__main__":
    main()
