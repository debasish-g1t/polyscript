"""Shared pytest configuration — runs before any test module is imported."""

import pytest


@pytest.fixture(autouse=True)
def _suppress_rdkit_logs():
    """Suppress RDKit warnings / logs during all tests."""
    from rdkit import RDLogger
    RDLogger.logger().setLevel(RDLogger.ERROR)
