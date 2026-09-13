"""
CLI entry point for DataFrame-based MPI polymerisation.

Delegates to ``MpiRunnerDf`` (``polyscript.polymerizer.mpi_df``).

Usage:
    mpirun -np 8 python mpi_runner_df.py <exp_name> --monomer-dfs monomer_dfs.pkl
"""

from __future__ import annotations

import argparse
import pickle
from pathlib import Path


def main() -> None:
    p = argparse.ArgumentParser(description="MPI polymerisation from DataFrames.")
    p.add_argument("exp_name", help="Experiment name.")
    p.add_argument("--monomer-dfs", required=True, type=Path,
                   help="Pickle file containing dict[str, DataFrame].")
    p.add_argument("--output-path", "-o", default="./outputs")
    p.add_argument("--rule-dir", "-r", default=None)
    p.add_argument("--df-chunk-size", type=int, default=1000)
    p.add_argument("--log-interval", type=int, default=100)
    p.add_argument("--no-rich", action="store_true")
    args = p.parse_args()

    with open(args.monomer_dfs, "rb") as fh:
        monomer_dfs = pickle.load(fh)

    from polyscript.polymerizer.mpi_df import MpiRunnerDf

    runner = MpiRunnerDf(
        exp_name=args.exp_name,
        monomer_dfs=monomer_dfs,
        output_path=args.output_path,
        rule_dir=args.rule_dir,
        df_chunk_size=args.df_chunk_size,
        log_interval=args.log_interval,
        use_rich=not args.no_rich,
    )
    stats = runner.run()
    if runner._is_master:
        print(f"Done. completed={stats['completed']} failed={stats['failed']}")


if __name__ == "__main__":
    main()
