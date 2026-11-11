"""MPI polymerizer integration test.

Runs both as a standalone script and as a pytest test.  When run under
pytest, each test function launches ``mpirun`` as a subprocess.

Uses the pre-built ``cls_split_10/`` split-data directory (one small
batch per monomer type, 1–10 SMILES each).

Standalone:
    python distrib/tests/test_mpi_polymerize.py

pytest:
    pytest distrib/tests/test_mpi_polymerize.py -v
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
# Resolve project paths
# ---------------------------------------------------------------------------
_HERE = Path(__file__).resolve().parent
_PROJECT = _HERE.parent  # distrib/
_SPLIT_DATA = _HERE / "test_data" / "cls_split_10"


# ---------------------------------------------------------------------------
# Temp runner script (auto-generated for each test run)
# ---------------------------------------------------------------------------

_RUNNER_TEMPLATE = '''\
"""Auto-generated MPI polymerizer test runner — do not edit."""

import os, sys
sys.path.insert(0, {project_dir!r})

from polyscript.polymerizer.mpi import MpiPolymerizer

runner = MpiPolymerizer(
    exp_name="mpi_test",
    input_data_split_dir={split_dir!r},
    output_path={out_dir!r},
    ckpt_path="ckpt.db",
    rule_dir=None,
    log_filepath="polymerizer.log",
    col_name="psmiles",
    use_rich=False,
)

stats = runner.run(p_classes={p_classes!r})
runner.close()

if runner._is_master:
    print(f"Stats: {{stats}}", file=sys.stderr)
    assert stats["completed"] > 0, "No tasks completed"
    assert stats["failed"] == 0, f"{{stats['failed']}} tasks failed"

    for p_class in {p_classes!r}:
        result_dir = os.path.join({out_dir!r}, "mpi_test", p_class)
        if os.path.isdir(result_dir):
            count = len([f for f in os.listdir(result_dir) if f.endswith(".parquet")])
            print(f"  {{p_class}}: {{count}} parquet files", file=sys.stderr)
        else:
            print(f"  {{p_class}}: no output dir", file=sys.stderr)

    print("PASSED", file=sys.stderr)
'''


def _run_mpi_test(
    tmp_root: str, p_classes: list[str]
) -> subprocess.CompletedProcess:
    """Launch the MPI polymerizer test inside *tmp_root*.

    1. Writes a self-contained runner script to the temp directory.
    2. Launches it via ``mpirun -np 2``.
    3. Returns the ``CompletedProcess`` for inspection.
    """
    out_dir = os.path.join(tmp_root, "outputs")
    runner_script = os.path.join(tmp_root, "_runner.py")
    with open(runner_script, "w") as fh:
        fh.write(
            _RUNNER_TEMPLATE.format(
                project_dir=str(_PROJECT),
                split_dir=str(_SPLIT_DATA),
                out_dir=out_dir,
                p_classes=p_classes,
            )
        )

    return subprocess.run(
        ["mpirun", "-np", "2", sys.executable, runner_script],
        cwd=tmp_root,
        capture_output=True,
        text=True,
        timeout=300,
    )


# ---------------------------------------------------------------------------
# Skip markers
# ---------------------------------------------------------------------------

_MPIRUN_MISSING = shutil.which("mpirun") is None
_MPI4PY_MISSING = False
try:
    import mpi4py  # noqa: F401
except ImportError:
    _MPI4PY_MISSING = True


# ---------------------------------------------------------------------------
# pytest tests
# ---------------------------------------------------------------------------


@pytest.mark.skipif(_MPIRUN_MISSING, reason="mpirun not found on PATH")
@pytest.mark.skipif(_MPI4PY_MISSING, reason="mpi4py not installed")
class TestMpiPolymerizer:
    """Integration tests that launch the MPI polymerizer via mpirun.

    Uses the pre-built ``cls_split_10/`` dataset — one tiny batch per
    monomer type, so the total task count is small and memory stays low.
    """

    def test_polyolefin_all_tasks_complete(self, tmp_path):
        """Run only the 'polyolefin' class and verify all tasks succeed."""
        result = _run_mpi_test(str(tmp_path), p_classes=["polyolefin"])

        assert result.returncode == 0, (
            f"mpirun failed (rc={result.returncode})\n"
            f"STDERR:\n{result.stderr}\nSTDOUT:\n{result.stdout}"
        )
        assert "PASSED" in result.stderr

    def test_polyester_all_tasks_complete(self, tmp_path):
        """Run only the 'polyester' class and verify all tasks succeed."""
        result = _run_mpi_test(str(tmp_path), p_classes=["polyester"])

        assert result.returncode == 0, (
            f"mpirun failed (rc={result.returncode})\n"
            f"STDERR:\n{result.stderr}\nSTDOUT:\n{result.stdout}"
        )
        assert "PASSED" in result.stderr


# ---------------------------------------------------------------------------
# Standalone entry point (auto-launch with mpirun when called directly)
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    tmp_root = tempfile.mkdtemp(prefix="mpi_polymerizer_test_")
    try:
        print("Launching mpirun -np 2 (polyolefin + polyester) …", file=sys.stderr)
        result = _run_mpi_test(tmp_root, p_classes=["polyolefin", "polyester"])
        if result.returncode != 0:
            print(result.stderr, file=sys.stderr)
        else:
            print(result.stderr, file=sys.stderr)
        sys.exit(result.returncode)
    finally:
        shutil.rmtree(tmp_root, ignore_errors=True)
