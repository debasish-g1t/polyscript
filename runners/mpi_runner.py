"""
CLI entry point for MPI-parallel polymerisation.

Delegates to ``MpiPolymerizer`` (``polyscript.polymerizer.mpi``).

Usage:
    mpirun -np 8 python mpi_runner.py <exp_name> <input_data_split_dir>
    mpirun -np 8 python mpi_runner.py my_run ./splits -o ./outputs -l DEBUG
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def main() -> None:
    p = argparse.ArgumentParser(
        description="MPI-parallel batch polymerisation runner."
    )
    p.add_argument("exp_name", help="Experiment name (subdirectory under output).")
    p.add_argument("input_data_split_dir", help="Directory with per-class data split files.")
    p.add_argument("--output-path", "-o", default="./outputs",
                   help="Root output directory (default: ./outputs).")
    p.add_argument("--ckpt-path", default="ckpt.db",
                   help="SQLite checkpoint filename (default: ckpt.db).")
    p.add_argument("--rule-dir", "-r", default=None,
                   help="Rules directory (default: auto-detected from polymerizer package).")
    p.add_argument("--log-filepath", default="polymerizer.log",
                   help="Log filename (default: polymerizer.log).")
    p.add_argument("--col-name", default="psmiles",
                   help="SMILES column name (default: psmiles).")
    p.add_argument("--chunk-size", type=int, default=None,
                   help="Max tasks per worker chunk (auto-tuned when omitted).")
    p.add_argument("--log-interval", type=int, default=100,
                   help="Min tasks between progress logs (default: 100).")
    p.add_argument("--no-rich", action="store_true",
                   help="Disable Rich progress bars.")
    p.add_argument("--log-level", "-l", default="INFO",
                   help="Log level: SILENT, WARNING, INFO, DEBUG.")
    args = p.parse_args()

    from polyscript.polymerizer.mpi import MpiPolymerizer

    runner = MpiPolymerizer(
        exp_name=args.exp_name,
        input_data_split_dir=args.input_data_split_dir,
        output_path=args.output_path,
        ckpt_path=args.ckpt_path,
        rule_dir=args.rule_dir,
        log_filepath=args.log_filepath,
        col_name=args.col_name,
        chunk_size=args.chunk_size,
        log_interval=args.log_interval,
        use_rich=not args.no_rich,
    )
    runner.logger.set_level(args.log_level)

    stats = runner.run()
    if runner._is_master:
        print(f"Done. completed={stats['completed']} failed={stats['failed']}")


if __name__ == "__main__":
    main()
