"""Tests for polymerizer bipolymerization and post-processing."""

from pathlib import Path

import pandas as pd
import pytest
from polyscript.polymerizer import Polymerizer

# ---------------------------------------------------------------------------
# Paths — use the cls_split_10 dataset (mirrors production split layout)
# ---------------------------------------------------------------------------

_CLS_SPLIT = Path(__file__).resolve().parent / "test_data" / "cls_split_10"

# Map: short alias → (subdir, filename) inside cls_split_10/
_SPLIT_FILES = {
    "vinyl":    ("vinyl_500",    "vinyl_0.parquet"),
    "cOle":     ("cOle_500",     "cOle_0.parquet"),
    "diol":     ("diol_500",     "diol_0.parquet"),
    "CO":       ("CO_500",       "CO_0.parquet"),
    "hydCOOH":  ("hydCOOH_500",  "hydCOOH_0.parquet"),
}


def _load_monomer(alias: str) -> pd.DataFrame:
    """Load a monomer batch from the cls_split_10 dataset."""
    subdir, fname = _SPLIT_FILES[alias]
    return pd.read_parquet(_CLS_SPLIT / subdir / fname)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def polymerizer():
    """Module-scoped: loading rule files is expensive, so reuse across tests."""
    return Polymerizer()


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


class TestBipolymerization:
    """Bipolymerization with two monomer DataFrames."""

    # ── vinyl + vinyl ──────────────────────────────────────────────────

    @pytest.fixture
    def vinyl_vinyl_result(self, polymerizer):
        mon_df1 = _load_monomer("vinyl")
        mon_df2 = _load_monomer("vinyl")
        return polymerizer.bipolymerize(
            mon_df1, "vinyl", "polyolefin",
            mon_df2=mon_df2, mon_type2="vinyl",
            candidate_col_name="psmiles",
            output_path=None,
        )

    def test_vinyl_vinyl_shape(self, vinyl_vinyl_result):
        assert vinyl_vinyl_result.shape == (90, 4)

    def test_vinyl_vinyl_first_polym(self, vinyl_vinyl_result):
        assert vinyl_vinyl_result["polym"][0] == (
            "{C=Cc1ccc(-c2cc3ccccc3s2)c(F)c1}+"
            "{C=C(OC1C(=O)C(O)=C(c2ccccc2)S1(=O)=O)c1nccnc1Cl}=>"
            "[C&X3;H2,H1,H0;!R:1]=[C&X3;H2,H1,H0;!R:2]."
            "[C&X3;H2,H1,H0;!R:3]=[C&X3;H2,H1,H0;!R:4]>>"
            "*-[C&X4:1][C&X4:2][C&X4:3][C&X4:4]-*=>"
            "*CC(CC(*)(OC1C(=O)C(O)=C(c2ccccc2)S1(=O)=O)c1nccnc1Cl)"
            "c1ccc(-c2cc3ccccc3s2)c(F)c1"
        )

    # ── vinyl + cOle ───────────────────────────────────────────────────

    @pytest.fixture
    def vinyl_cole_result(self, polymerizer):
        mon_df1 = _load_monomer("vinyl")
        mon_df2 = _load_monomer("cOle")
        return polymerizer.bipolymerize(
            mon_df1, "vinyl", "polyolefin",
            mon_df2=mon_df2, mon_type2="cOle",
            candidate_col_name="psmiles",
            output_path=None,
        )

    def test_vinyl_cole_shape(self, vinyl_cole_result):
        assert vinyl_cole_result.shape == (324, 4)

    def test_vinyl_cole_first_polym(self, vinyl_cole_result):
        assert vinyl_cole_result["polym"][0] == (
            "{C=Cc1ccc(-c2cc3ccccc3s2)c(F)c1}+"
            "{O=C1CSC=CN1c1cc(N=Nc2ccccc2F)ccc1O}=>"
            "[C&X3;H2,H1,H0;!R:1]=[C&X3;H2,H1,H0;!R:2]."
            "[C&X3;H1,H0;R:3]=[C&X3;H1,H0;R:4]>>"
            "*-[C&X4:1][C&X4:2][C&X4:3][C&X4:4]-*=>"
            "*CC(c1ccc(-c2cc3ccccc3s2)c(F)c1)C1SCC(=O)"
            "N(c2cc(N=Nc3ccccc3F)ccc2O)C1*"
        )

    def test_vinyl_cole_polym_at_313(self, vinyl_cole_result):
        assert vinyl_cole_result["polym"][313] == (
            "{C=C(CBr)Cc1cc2ccc(O)cc2oc1=O}+"
            "{C=C(OC1(Cl)CC2C=CC1C2)c1cccc(N)c1Cl}=>"
            "[C&X3;H2,H1,H0;!R:1]=[C&X3;H2,H1,H0;!R:2]."
            "[C&X3;H1,H0;R:3]=[C&X3;H1,H0;R:4]>>"
            "*-[C&X4:1][C&X4:2][C&X4:3][C&X4:4]-*=>"
            "*C1C2CC(C1CC(*)(CBr)Cc1cc3ccc(O)cc3oc1=O)C(Cl)"
            "(OC(=C)c1cccc(N)c1Cl)C2|||"
            "{*C1C2CC(C1CC(*)(CBr)Cc1cc3ccc(O)cc3oc1=O)C(Cl)"
            "(OC(=C)c1cccc(N)c1Cl)C2}+{none}=>"
            "[C&X3:1]=[C&X3:2]>>"
            "*-[C&X4:1][C&X4:2]-*=>"
            "*CC(*)(OC1(Cl)CC2CC1C(CC(*)(CBr)Cc1cc3ccc(O)cc3oc1=O)C2*)"
            "c1cccc(N)c1Cl"
        )

    # ── diol + CO ──────────────────────────────────────────────────────

    @pytest.fixture
    def diol_co_result(self, polymerizer):
        mon_df1 = _load_monomer("diol")
        mon_df2 = _load_monomer("CO")
        return polymerizer.bipolymerize(
            mon_df1, "diol", "polyester",
            mon_df2=mon_df2, mon_type2="CO",
            candidate_col_name="psmiles",
            output_path=None,
        )

    def test_diol_co_shape(self, diol_co_result):
        assert diol_co_result.shape == (20, 4)

    def test_diol_co_first_polym(self, diol_co_result):
        assert diol_co_result["polym"][0] == (
            "{Cc1cc(-c2cc(Cl)ccc2O)ncc1O}+{[C-]#[O+]}=>"
            "([O,S;X2;H1;!$([O,S]C=*):1].[O,S;X2;H1;!$([O,S]C=*):2])."
            "[C&-]#[O&+]>>"
            "(*-[O,S;X2;!$([O,S]C=*):1]."
            "[O,S;X2;!$([O,S]C=*):2][C&X3](=O)-*)=>"
            "*Oc1ccc(Cl)cc1-c1cc(C)c(OC(*)=O)cn1"
        )


class TestHomopolymerization:
    """Homopolymerization with a single monomer DataFrame."""

    @pytest.fixture
    def hydcooh_none_result(self, polymerizer):
        mon_df1 = _load_monomer("hydCOOH")
        return polymerizer.bipolymerize(
            mon_df1, "hydCOOH", "polyester",
            mon_df2=None, mon_type2="none",
            candidate_col_name="psmiles",
            output_path=None,
        )

    def test_hydcooh_none_shape(self, hydcooh_none_result):
        assert hydcooh_none_result.shape == (10, 4)

    def test_hydcooh_none_first_polym(self, hydcooh_none_result):
        assert hydcooh_none_result["polym"][0] == (
            "{C=C(C(=O)O)C1CN(c2ccc(Cl)nc2O)C1}+{none}=>"
            "([O&X2&H1&!$(OC=*):1].[C&X3:2](=O)[O&X2&H1])>>"
            "(*-[O&X2:1].[C&X3:2](=O)-*)=>"
            "*Oc1nc(Cl)ccc1N1CC(C(=C)C(*)=O)C1"
        )


class TestPostProcess:
    """Post-processing of polymer results (Polars / pandas fallback)."""

    @pytest.fixture
    def raw_result(self, polymerizer):
        mon_df1 = _load_monomer("vinyl")
        mon_df2 = _load_monomer("cOle")
        return polymerizer.bipolymerize(
            mon_df1, "vinyl", "polyolefin",
            mon_df2=mon_df2, mon_type2="cOle",
            candidate_col_name="psmiles",
            output_path=None,
        )

    @pytest.fixture
    def processed(self, polymerizer, raw_result):
        return polymerizer.post_process(
            raw_result, drop_outro=True, reduce_column=["mon1", "mon2"]
        )

    def test_post_process_shape(self, processed):
        assert processed.shape == (324, 2)

    def test_post_process_first_rep_kept(self, processed):
        """Nested "|||" representations are reduced to the first one."""
        assert processed["polym"][313] == (
            "{C=C(CBr)Cc1cc2ccc(O)cc2oc1=O}+"
            "{C=C(OC1(Cl)CC2C=CC1C2)c1cccc(N)c1Cl}=>"
            "[C&X3;H2,H1,H0;!R:1]=[C&X3;H2,H1,H0;!R:2]."
            "[C&X3;H1,H0;R:3]=[C&X3;H1,H0;R:4]>>"
            "*-[C&X4:1][C&X4:2][C&X4:3][C&X4:4]-*=>"
            "*C1C2CC(C1CC(*)(CBr)Cc1cc3ccc(O)cc3oc1=O)C(Cl)"
            "(OC(=C)c1cccc(N)c1Cl)C2"
        )
