"""Tests for CanonicalizationAnalyzer."""

from pathlib import Path

import pytest
from polyscript.utils.canon_analyze import CanonicalizationAnalyzer

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

_TEST_DATA = Path(__file__).resolve().parent.parent / "test_data" / "sample_outputs"


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def analyzer():
    """Module-scoped: scanning parquet files is expensive, reuse across tests."""
    return CanonicalizationAnalyzer(base_dir=str(_TEST_DATA))


@pytest.fixture(scope="module")
def full_result(analyzer):
    """Collected + analysed result, computed once per module."""
    return analyzer.analyse(analyzer.collect())


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


class TestCanonicalizationAnalyzer:
    """Full analysis and max-rows-limited analysis."""

    def test_full_analysis(self, full_result):
        """Scan all parquet files and verify duplication stats."""
        assert full_result == {
            "unique": 2024,
            "total": 2052,
            "dup_groups": 14,
            "same_monomer": 0,
            "diff_monomer": 14,
        }

    def test_max_rows_limits_collection(self):
        """With max_rows=10, each class collects at most 10 rows."""
        analyzer = CanonicalizationAnalyzer(
            base_dir=str(_TEST_DATA), max_rows=10
        )
        result = analyzer.analyse(analyzer.collect())
        assert result == {
            "unique": 70,
            "total": 70,
            "dup_groups": 0,
            "same_monomer": 0,
            "diff_monomer": 0,
        }
