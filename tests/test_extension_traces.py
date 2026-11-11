"""Smoke tests for the conventional extension-trace catalogue.

The trace generator writes documentation and SVG files, so these tests only
build the selections and traces in memory.  They keep the known-polymer
examples aligned with the configured recipes without regenerating the page.
"""

import sys
from pathlib import Path

_PROJECT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_PROJECT / "benchmarks"))

import extract_extension_traces as traces  # noqa: E402


def test_conventional_extra_specs_cover_the_ring_and_romp_families():
    kind_by_prefix = {
        spec["prefix"]: spec["kind"]
        for spec in traces._CONVENTIONAL_EXTRA
    }
    assert kind_by_prefix["peo"] == "step_growth"
    assert kind_by_prefix["pcl"] == "ring_opening"
    assert kind_by_prefix["nylon6"] == "ring_opening"
    assert kind_by_prefix["pes"] == "ring_opening"
    assert kind_by_prefix["ppo"] == "step_growth"
    assert kind_by_prefix["romp"] == "romp"


def test_conventional_extra_traces_build_cleanly():
    for spec in traces._CONVENTIONAL_EXTRA:
        selection = traces._spec_selection(spec)
        result = traces._trace(selection["kind"], selection)

        assert result["intermediates"], spec["polymer_name"]
        assert result["final"], spec["polymer_name"]
        assert all(
            "*" not in smiles for smiles in result["final"]
        ), spec["polymer_name"]
        for step in result["intermediates"]:
            assert step
            assert all(isinstance(smiles, str) for smiles in step)


def test_known_polymer_names_are_present():
    names = {
        spec["polymer_name"] for spec in traces._CONVENTIONAL_EXTRA
    }
    assert "Poly(ethylene oxide) (PEO)" in names
    assert "Polycaprolactone (PCL)" in names
    assert "Nylon-6" in names
    assert "Poly(ethylene succinate) (PES)" in names
    assert (
        "Poly(2,6-dimethyl-1,4-phenylene oxide) (PPO)" in names
    )
    assert "Poly(norbornene) (PNB)" in names
