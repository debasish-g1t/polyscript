"""MPI sequence validator integration test.

Launches ``MpiSeqValidator`` via ``mpirun -np 2``, validating the
sample polymerizer outputs in ``test_data/sample_outputs/``.

Standalone:
    python distrib/tests/utils/test_mpi_validator.py

pytest:
    pytest distrib/tests/utils/test_mpi_validator.py -v
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pytest

# ---------------------------------------------------------------------------
# Resolve paths
# ---------------------------------------------------------------------------
_HERE = Path(__file__).resolve().parent
_PROJECT = _HERE.parent.parent  # distrib/
_TEST_DATA = _HERE.parent / "test_data"
_INPUT_DIR = _TEST_DATA / "sample_outputs"
_RULE_DIR = _PROJECT / "polyscript" / "polymerizer" / "rules"


# ---------------------------------------------------------------------------
# Temp runner script
# ---------------------------------------------------------------------------

_RUNNER = '''\
import os, sys
sys.path.insert(0, {project_dir!r})

from polyscript.utils.mpi_validator import MpiSeqValidator

validator = MpiSeqValidator(
    input_dir={input_dir!r},
    output_dir={out_dir!r},
    rule_dir={rule_dir!r},
    log_filepath="validator.log",
    use_rich=False,
)

stats = validator.run(p_classes={p_classes!r})
validator.close()

if validator._is_master:
    print(f"Stats: {{stats}}", file=sys.stderr)
    assert stats["completed"] > 0, "No tasks completed"
    assert stats["total_rows"] > 0, "No rows processed"

    # Verify expected output files exist
    output_dir = {out_dir!r}
    for fname in (
        os.path.join("linear_sets", "polyolefin_linear.parquet"),
        os.path.join("statistics", "validation_stats.json"),
    ):
        path = os.path.join(output_dir, fname)
        if os.path.exists(path):
            print(f"  OK: {{fname}}", file=sys.stderr)
        else:
            print(f"  MISSING: {{fname}}", file=sys.stderr)

    print("PASSED", file=sys.stderr)
'''


def _run_validator(
    tmp_root: str, p_classes: list[str]
) -> subprocess.CompletedProcess:
    """Launch the MPI validator inside *tmp_root*."""
    out_dir = os.path.join(tmp_root, "outputs")
    runner_script = os.path.join(tmp_root, "_runner.py")
    with open(runner_script, "w") as fh:
        fh.write(
            _RUNNER.format(
                project_dir=str(_PROJECT),
                input_dir=str(_INPUT_DIR),
                out_dir=out_dir,
                rule_dir=str(_RULE_DIR),
                p_classes=p_classes,
            )
        )

    return subprocess.run(
        ["mpirun", "-np", "2", sys.executable, runner_script],
        cwd=tmp_root,
        capture_output=True,
        text=True,
        timeout=120,
    )


# ---------------------------------------------------------------------------
# Skip markers
# ---------------------------------------------------------------------------

_MPIRUN = shutil.which("mpirun")
_MPI4PY = False
try:
    import mpi4py  # noqa: F401
    _MPI4PY = True
except ImportError:
    pass


# ---------------------------------------------------------------------------
# pytest tests
# ---------------------------------------------------------------------------


@pytest.mark.skipif(_MPIRUN is None, reason="mpirun not found on PATH")
@pytest.mark.skipif(not _MPI4PY, reason="mpi4py not installed")
class TestMpiSeqValidator:
    """Integration tests that launch the MPI validator via mpirun."""

    def test_validate_sample_outputs(self, tmp_path):
        """Validate polymerizer sample outputs and verify results."""
        result = _run_validator(str(tmp_path), p_classes=["polyolefin"])

        assert result.returncode == 0, (
            f"mpirun failed (rc={result.returncode})\n"
            f"STDERR:\n{result.stderr}\nSTDOUT:\n{result.stdout}"
        )
        assert "PASSED" in result.stderr

        # Verify key output files exist
        out_dir = os.path.join(str(tmp_path), "outputs")
        assert os.path.isfile(
            os.path.join(out_dir, "linear_sets", "polyolefin_linear.parquet")
        ), "Missing merged parquet output"
        assert os.path.isfile(
            os.path.join(out_dir, "statistics", "validation_stats.json")
        ), "Missing validation stats"


# ---------------------------------------------------------------------------
# Standalone entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    tmp_root = tempfile.mkdtemp(prefix="mpi_validator_test_")
    try:
        print("Launching mpirun -np 2 (polyolefin) …", file=sys.stderr)
        result = _run_validator(tmp_root, p_classes=["polyolefin"])
        if result.returncode != 0:
            print(result.stderr, file=sys.stderr)
            sys.exit(result.returncode)
        print(result.stderr, file=sys.stderr)
    finally:
        shutil.rmtree(tmp_root, ignore_errors=True)
