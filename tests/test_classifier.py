"""Tests for PolyScriptClassifier and the classifier→MPI-polymerizer pipeline."""

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest
from polyscript.classifier import PolyScriptClassifier

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
_HERE = Path(__file__).resolve().parent
_PROJECT = _HERE.parent  # distrib/
_CLS_SPLIT = _HERE / "test_data" / "cls_split_10"


def _load_split_smiles(mon_type: str, n: int = 2) -> list[str]:
    """Load the first *n* SMILES from a cls_split_10 monomer batch."""
    df = pd.read_parquet(_CLS_SPLIT / f"{mon_type}_500" / f"{mon_type}_0.parquet")
    return df["psmiles"].head(n).tolist()


# ---------------------------------------------------------------------------
# MPI skip markers
# ---------------------------------------------------------------------------
_MPIRUN = shutil.which("mpirun")
_MPI4PY = False
try:
    import mpi4py  # noqa: F401
    _MPI4PY = True
except ImportError:
    pass

# ---------------------------------------------------------------------------
# MPI runner template
# ---------------------------------------------------------------------------

_MPI_DF_RUNNER = '''\
import os, sys
sys.path.insert(0, {project_dir!r})

import pandas as pd
from polyscript.polymerizer.mpi_df import MpiRunnerDf

df = pd.read_parquet({classified_path!r})

runner = MpiRunnerDf(
    exp_name="mpi_df_test",
    df=df,
    smiles_col="smiles",
    monomer_classes={monomer_classes!r},
    output_path={out_dir!r},
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
        result_dir = os.path.join({out_dir!r}, "mpi_df_test", p_class)
        if os.path.isdir(result_dir):
            count = len([
                f for f in os.listdir(result_dir) if f.endswith(".parquet")
            ])
            print(f"  {{p_class}}: {{count}} parquet files", file=sys.stderr)
    print("PASSED", file=sys.stderr)
'''


class TestPolyScriptClassifier:
    """Tests for the PolyScriptClassifier."""

    @pytest.fixture
    def sample_df(self):
        """A minimal DataFrame with two amino-acid-like SMILES."""
        datasets = ["CNC(C)(C(=O)O)C(C)(C)C", "CNC(CC(=O)O)C(=O)O"]
        return pd.DataFrame({"smiles": datasets})

    # ── carbonate inclusion ────────────────────────────────────────────

    def test_classify_with_carbonate_adds_co_and_hcho_rows(self, sample_df):
        """When include_carbonate=True, CO and HCHO rows are appended."""
        classifier = PolyScriptClassifier(include_carbonate=True)
        result = classifier.classify(sample_df, "smiles")
        assert result.shape[0] == 4, (
            f"Expected 4 rows (2 input + 2 carbonate), got {result.shape[0]}"
        )

    def test_classify_without_carbonate_keeps_original_row_count(self, sample_df):
        """When include_carbonate=False, row count should stay the same."""
        classifier = PolyScriptClassifier(include_carbonate=False)
        result = classifier.classify(sample_df, "smiles")
        assert result.shape[0] == 2, (
            f"Expected 2 rows, got {result.shape[0]}"
        )

    # ── column values ──────────────────────────────────────────────────

    @pytest.fixture
    def classified_df(self, sample_df):
        """Return a DataFrame already classified without carbonate rows."""
        return PolyScriptClassifier(include_carbonate=False).classify(
            sample_df, "smiles"
        )

    def test_aminCOOH_column(self, classified_df):
        """Both monomers should be classified as amino-carboxylic acids."""
        assert classified_df["aminCOOH"].to_dict() == {0: True, 1: True}

    def test_diCOOH_column(self, classified_df):
        """Only CNC(CC(=O)O)C(=O)O has two carboxylic acid groups."""
        assert classified_df["diCOOH"].to_dict() == {0: False, 1: True}

    def test_di_acid_chloride_column(self, classified_df):
        """Neither monomer contains acid chloride groups."""
        assert classified_df["di_acid_chloride"].to_dict() == {0: False, 1: False}


# ---------------------------------------------------------------------------
# Classifier → MPI polymerizer pipeline (runs via mpirun subprocess)
# ---------------------------------------------------------------------------


@pytest.mark.skipif(_MPIRUN is None, reason="mpirun not found on PATH")
@pytest.mark.skipif(not _MPI4PY, reason="mpi4py not installed")
class TestClassifierToMpiPolymerizer:
    """End-to-end: classify SMILES, then polymerize via MpiRunnerDf."""

    @pytest.fixture
    def multi_type_smiles_df(self):
        """Build a DataFrame with SMILES from multiple monomer types.

        Picks 1 SMILES each from vinyl, cOle, and diol batches so the
        classifier produces True columns for all three monomer types
        while keeping the test fast.
        """
        smi = []
        for mt in ("vinyl", "cOle", "diol"):
            smi.extend(_load_split_smiles(mt, n=1))
        return pd.DataFrame({"smiles": smi})

    @pytest.fixture
    def classified_multi(self, multi_type_smiles_df):
        """The classifier output ready for MpiRunnerDf."""
        return PolyScriptClassifier(include_carbonate=False).classify(
            multi_type_smiles_df, "smiles"
        )

    def test_pipeline_classify_then_mpi_polymerize(
        self, classified_multi, tmp_path
    ):
        """Full pipeline: classify → MpiRunnerDf → verify outputs."""
        tmp_root = str(tmp_path)
        out_dir = os.path.join(tmp_root, "outputs")

        # Save classified DataFrame so the MPI subprocess can load it
        classified_parquet = os.path.join(tmp_root, "classified.parquet")
        classified_multi.to_parquet(classified_parquet, index=False)

        # Write and launch the MPI runner
        runner_script = os.path.join(tmp_root, "_runner.py")
        with open(runner_script, "w") as fh:
            fh.write(
                _MPI_DF_RUNNER.format(
                    project_dir=str(_PROJECT),
                    classified_path=classified_parquet,
                    out_dir=out_dir,
                    p_classes=["polyolefin"],
                    monomer_classes=["vinyl", "cOle"],
                )
            )

        result = subprocess.run(
            ["mpirun", "-np", "2", sys.executable, runner_script],
            cwd=tmp_root,
            capture_output=True,
            text=True,
            timeout=60,
        )

        assert result.returncode == 0, (
            f"mpirun failed (rc={result.returncode})\n"
            f"STDERR:\n{result.stderr}\nSTDOUT:\n{result.stdout}"
        )
        assert "PASSED" in result.stderr

        # Verify output parquet files were produced
        res_dir = os.path.join(out_dir, "mpi_df_test", "polyolefin")
        assert os.path.isdir(res_dir), f"Missing: {res_dir}"
        pq_files = [f for f in os.listdir(res_dir) if f.endswith(".parquet")]
        assert len(pq_files) > 0, f"No parquet output in {res_dir}"
