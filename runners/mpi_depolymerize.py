"""
CLI entry point for MPI-parallel depolymerization.

Delegates to ``MpiDepolymerizer`` (``polyscript.depolymerizer.mpi``).

Usage:
    mpirun -np 8 python mpi_depolymerize.py --input benchmark.parquet
    mpirun -np 8 python mpi_depolymerize.py --input benchmark.parquet --mode convert -l DEBUG
"""

from __future__ import annotations

import argparse
from pathlib import Path

from polyscript.depolymerizer.mpi import MpiDepolymerizer


def main() -> None:
    p = argparse.ArgumentParser(
        description="MPI-parallel depolymerization benchmark."
    )
    p.add_argument("--input", "-i", type=Path, required=True,
                   help="Parquet or CSV file with polymer sequences.")
    p.add_argument("--mode", "-m", choices=["validate", "convert"],
                   default="validate",
                   help="'validate' for PolyScript, 'convert' for raw PSMILES.")
    p.add_argument("--smiles-col", type=str, default=None,
                   help="SMILES column name (auto-detected in convert mode).")
    p.add_argument("--output-dir", "-o", type=Path, default=Path("."),
                   help="Output directory for results and logs.")
    p.add_argument("--log-level", "-l", type=str, default="INFO",
                   help="Log level: SILENT, WARNING, INFO, DEBUG.")
    args = p.parse_args()

    runner = MpiDepolymerizer(
        input_path=str(args.input),
        mode=args.mode,
        smiles_col=args.smiles_col,
        output_dir=str(args.output_dir),
        log_filepath="depolymerizer.log",
    )
    runner.logger.set_level(args.log_level)

    stats = runner.run()
    if runner._is_master:
        print(f"\nDone. success={stats.get('success')} failed={stats.get('failed')}")


if __name__ == "__main__":
    main()
