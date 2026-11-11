"""Unit tests for the benchmark harness helpers.

The benchmark is a development tool, but its provenance and baseline logic is
what makes the reported numbers trustworthy.  These tests keep that logic
honest without running the full 90-second evaluation.
"""

import json
import sys
from collections import Counter
from pathlib import Path

from rdkit import Chem

sys.path.insert(0, str(Path(__file__).resolve().parent))

import benchmark_extender as bench  # noqa: E402

from polyscript.extender.map import rxn_map  # noqa: E402


def _key_for_coupling(name: str) -> str:
    for key, config in rxn_map.items():
        if config and name in bench._coupling_names(key, config):
            return key
    raise AssertionError(f"no configured reaction uses {name!r}")


def test_end_group_catalogue_keeps_acid_chloride_terminal():
    labels = bench._end_group_labels(
        Chem.MolFromSmiles("O=C(Cl)CC(=O)OCCO")
    )
    assert "acid_halide" in labels
    assert "carboxylic_acid" not in labels
    assert "alcohol" in labels


def test_end_group_catalogue_recognises_imide_anhydride_end():
    labels = bench._end_group_labels(
        Chem.MolFromSmiles(
            "Nc1ccc(-n2c(=O)c3cc4c(=O)oc(=O)c4cc3c2=O)cc1"
        )
    )
    assert {"primary_amine", "anhydride"} <= labels


def test_documented_new_end_groups_are_explicit():
    epoxide_key = _key_for_coupling("epoxide_rop")
    assert bench._allowed_new_end_groups(
        rxn_map[epoxide_key]
    ) == {"alcohol"}

    lactone_key = _key_for_coupling("lactone_rop")
    assert bench._allowed_new_end_groups(
        rxn_map[lactone_key]
    ) == {"anhydride"}

    allowed_keys = {epoxide_key, lactone_key}
    for other_key, other_config in rxn_map.items():
        if other_config is None or other_key in allowed_keys:
            continue
        assert not bench._allowed_new_end_groups(other_config)


def test_atom_provenance_allowance_tracks_mapped_atoms():
    ester_key = _key_for_coupling("ester")
    ester_probe = bench._build_probe_driver(
        ester_key, rxn_map[ester_key], max_products=1
    )
    # The ester carbonyl oxygen is mapped through the reaction rather than
    # being recreated, so the template has no allowance for new heavy atoms.
    assert bench._allowed_new_atoms(ester_probe) == Counter()

    urethane_key = _key_for_coupling("urethane")
    urethane_probe = bench._build_probe_driver(
        urethane_key, rxn_map[urethane_key], max_products=1
    )
    # Every product atom is mapped back to a reactant atom.
    assert bench._allowed_new_atoms(urethane_probe) == Counter()


def test_baseline_diff_separates_regressions_and_fixes(tmp_path):
    baseline = {
        "record_index": {
            "f|0|2": {"ok": True, "failure": ""},
            "f|1|2": {"ok": False, "failure": "no_product"},
        }
    }
    path = tmp_path / "summary.json"
    path.write_text(json.dumps(baseline))

    records = [
        {
            "source_file": "f",
            "row_index": 0,
            "cycle": 2,
            "ok": False,
            "failure": "no_product",
        },
        {
            "source_file": "f",
            "row_index": 1,
            "cycle": 2,
            "ok": True,
            "failure": "",
        },
    ]
    diff = bench._compare_baseline(records, path)
    assert diff["regressions"] == ["f|0|2"]
    assert diff["fixes"] == ["f|1|2"]
    assert diff["regressions_by_cycle"] == {"2": 1}
    assert diff["fixes_by_cycle"] == {"2": 1}
    assert diff["unchanged_failures"] == []


def test_cycle_summary_separates_known_and_new_failures():
    base = {
        "warnings": [],
        "runtime_ms": 1.0,
        "heavy_atoms": 10,
        "delta_heavy": None,
        "n_products": 1,
    }
    records = [
        {**base, "ok": True, "failure": "", "known_gap": ""},
        {
            **base,
            "ok": False,
            "failure": "no_product",
            "known_gap": "",
            "n_products": 0,
        },
        {
            **base,
            "ok": False,
            "failure": "no_product",
            "known_gap": "known reason",
            "n_products": 0,
        },
    ]
    summary = bench._summarise_cycle(records)
    assert summary["attempted"] == 3
    assert summary["passed"] == 1
    assert summary["pass_rate"] == round(1 / 3, 4)
    assert summary["new_failures"] == 1
    assert summary["known_gap_failures"] == 1
    assert summary["failures"] == {"no_product": 2}
    assert summary["mean_products"] == round(1 / 3, 3)
    assert summary["max_products_observed"] == 1
    assert summary["rows_with_multiple_products"] == 0
