"""
CLI entry point for MPI-parallel sequence validation.

Delegates to ``MpiSeqValidator`` (``polyscript.utils.mpi_validator``).

Usage:
    mpirun -np 8 python mpi_validator.py --input-dir ./outputs --output-dir ./validated --rule-dir ./rules
"""

from __future__ import annotations

import argparse


def main() -> None:
    ap = argparse.ArgumentParser(description="MPI-parallel PolyScript validator")
    ap.add_argument("--input-dir", required=True)
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--rule-dir", required=True)
    ap.add_argument("--ckpt-path", default=None)
    ap.add_argument("--log-file", default="validator.log")
    ap.add_argument("--log-interval", type=int, default=10)
    ap.add_argument("--no-rich", action="store_true")
    ap.add_argument("--p-classes", nargs="*", default=None)
    args = ap.parse_args()

    from polyscript.utils.mpi_validator import MpiSeqValidator

    validator = MpiSeqValidator(
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        rule_dir=args.rule_dir,
        ckpt_path=args.ckpt_path,
        log_filepath=args.log_file,
        log_interval=args.log_interval,
        use_rich=not args.no_rich,
    )
    result = validator.run(p_classes=args.p_classes)
    if result is not None:
        print(f"\nDone. {result['completed']} files processed.")
        print(f"Valid: {result['valid']:,}  Invalid: {result['invalid']:,}")


if __name__ == "__main__":
    main()
