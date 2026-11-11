"""CSV-driven long-chain tests for the addition reaction driver.

Every polyolefin row in ``tests/test_data/extender_test_data_all.csv`` is matched
to a configured entry in ``rxn_map``.  The tests grow each row's monomer(s) to
increasing lengths with the frontier cap enabled, so they stay fast, and assert
structural invariants that must hold when a chain is extended one step at a
time:

* the heavy-atom count grows by exactly one monomer (or one monomer pair);
* the ring count grows by exactly the rings present in that unit;
* the terminated product contains no wildcard;
* the reference ``polym`` pSMILES in the CSV contains exactly the monomer
  heavy atoms and rings.

Initiation and termination are hydrogen caps (the wildcard is simply removed),
so they add no heavy atoms to the chain.
"""

from itertools import pairwise
from pathlib import Path

import pandas as pd
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from polyscript.extender import ChainGrowthReaction
from polyscript.extender.map import rxn_map

_CSV_PATH = (
    Path(__file__).resolve().parent
    / "test_data"
    / "extender_test_data_all.csv"
)
_DF = pd.read_csv(_CSV_PATH)
_POLYOLEFIN = _DF[_DF["polymer_class"] == "polyolefin"].reset_index(drop=True)
_MAX_LENGTH = 5


def _mol(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None, f"unparsable SMILES: {smiles}"
    return mol


def _free_slots(config) -> int:
    """Number of monomer slots in the propagation reaction (chain excluded)."""
    prop = AllChem.ReactionFromSmarts(config["params"]["prop_smarts"])
    return prop.GetNumReactantTemplates() - 1


def _possible_totals(first: int, second: int, pairs: int) -> set[int]:
    """Possible sums when each pair contributes 2*first, first+second, or 2*second."""
    totals = {0}
    for _ in range(pairs):
        totals = {
            total + contribution
            for total in totals
            for contribution in (2 * first, first + second, 2 * second)
        }
    return totals


def _row_data(row):
    config = rxn_map[row["reaction_smarts"]]
    mon1 = row["mon1"]
    mon2 = None if pd.isna(row["mon2"]) else row["mon2"]
    pool = [mon1] if mon2 is None else [mon1, mon2]

    mol1 = _mol(mon1)
    h1 = mol1.GetNumHeavyAtoms()
    r1 = mol1.GetRingInfo().NumRings()
    if mon2 is None:
        h2 = r2 = 0
    else:
        mol2 = _mol(mon2)
        h2 = mol2.GetNumHeavyAtoms()
        r2 = mol2.GetRingInfo().NumRings()

    return config, pool, mon1, mon2, h1, h2, r1, r2


def _ids() -> list[str]:
    return [
        f"{index}-{row['reaction_smarts'][:20]}"
        for index, row in _POLYOLEFIN.iterrows()
    ]


@pytest.mark.parametrize("row_index", _POLYOLEFIN.index, ids=_ids())
def test_polyolefin_grows_to_longer_lengths(row_index):
    row = _POLYOLEFIN.loc[row_index]
    config, pool, _mon1, mon2, h1, h2, r1, r2 = _row_data(row)
    assert config is not None
    free_slots = _free_slots(config)
    if free_slots == 1:
        assert mon2 is None
    else:
        assert mon2 is not None

    reaction = ChainGrowthReaction(**{**config["params"], "max_products": 1})
    heavy_counts = []
    ring_counts = []
    for length in range(_MAX_LENGTH + 1):
        products = reaction.extend(pool, length)
        assert len(products) == 1
        mol = _mol(products[0])
        assert not any(atom.GetAtomicNum() == 0 for atom in mol.GetAtoms())
        heavy_counts.append(mol.GetNumHeavyAtoms())
        ring_counts.append(mol.GetRingInfo().NumRings())

    for length, (heavy, rings) in enumerate(zip(heavy_counts, ring_counts)):
        if free_slots == 1:
            # Initiation and termination are hydrogen caps: no extra methyl
            # heavy atoms are added.
            assert heavy == (length + 1) * h1
            assert rings == (length + 1) * r1
        else:
            assert heavy in _possible_totals(h1, h2, length + 1)
            assert rings in _possible_totals(r1, r2, length + 1)

    # Consecutive lengths differ by exactly one added unit (or one added pair).
    for previous, current in pairwise(heavy_counts):
        if free_slots == 1:
            assert current - previous == h1
        else:
            assert current - previous in {2 * h1, h1 + h2, 2 * h2}
    for previous, current in pairwise(ring_counts):
        if free_slots == 1:
            assert current - previous == r1
        else:
            assert current - previous in {2 * r1, r1 + r2, 2 * r2}


@pytest.mark.parametrize("row_index", _POLYOLEFIN.index, ids=_ids())
def test_csv_reference_psmiles_matches_monomers(row_index):
    row = _POLYOLEFIN.loc[row_index]
    config, _pool, _mon1, _mon2, h1, h2, r1, r2 = _row_data(row)
    assert config is not None

    psmiles = row["polym"].rsplit("=>", 1)[-1]
    reference = _mol(psmiles)
    assert reference.GetNumHeavyAtoms() == h1 + h2
    assert reference.GetRingInfo().NumRings() == r1 + r2


def test_all_csv_reactions_are_registered():
    assert set(_DF["reaction_smarts"]) <= set(rxn_map)


def test_rxn_map_holds_lazy_smarts_data():
    """The registry must stay pure data; SMARTS compile only when loaded."""
    from rdkit.Chem.rdChemReactions import ChemicalReaction

    names_by_type = {
        "addition": {"init_smarts", "prop_smarts", "term_smarts"},
        "step_growth": {"coupling_smarts"},
        "ring_opening": {"prop_smarts"},
        "romp": {"initiator", "terminator"},
    }
    for config in rxn_map.values():
        assert not isinstance(config, ChemicalReaction)
        if config is None:
            continue

        expected = names_by_type[config["type"]]
        params = config["params"]
        assert expected <= set(params)

        if config["type"] == "addition":
            assert set(params) == expected
            assert all(isinstance(params[name], str) for name in expected)
        elif config["type"] == "romp":
            assert set(params) == expected
            assert all(isinstance(params[name], str) for name in expected)
        else:
            key = (
                "prop_smarts"
                if config["type"] == "ring_opening"
                else "coupling_smarts"
            )
            assert isinstance(params[key], list)
            assert params[key]
            assert all(
                isinstance(smarts, str)
                for smarts in params[key]
            )
            if "reagents" in params:
                assert isinstance(params["reagents"], list)
                assert all(
                    isinstance(smiles, str) for smiles in params["reagents"]
                )
            if "opening_reagents" in params:
                assert isinstance(params["opening_reagents"], list)
                assert all(
                    isinstance(smiles, str)
                    for smiles in params["opening_reagents"]
                )


def test_polyolefin_rows_are_configured():
    for smarts in _POLYOLEFIN["reaction_smarts"]:
        config = rxn_map[smarts]
        assert config is not None
        assert config["polymer_class"] == "polyolefin"
        if config["type"] == "romp":
            assert set(config["params"]) == {"initiator", "terminator"}
            continue
        assert config["type"] == "addition"
        assert set(config["params"]) == {
            "init_smarts",
            "prop_smarts",
            "term_smarts",
        }


_PAIR_ROWS = [
    index
    for index, row in _POLYOLEFIN.iterrows()
    if _free_slots(rxn_map[row["reaction_smarts"]]) == 2
]


@pytest.mark.parametrize(
    "row_index", _PAIR_ROWS, ids=[f"pair-{index}" for index in _PAIR_ROWS]
)
def test_pair_rows_respect_frontier_cap(row_index):
    row = _POLYOLEFIN.loc[row_index]
    config, pool, *_ = _row_data(row)
    reaction = ChainGrowthReaction(**{**config["params"], "max_products": 3})
    for length in range(_MAX_LENGTH + 1):
        products = reaction.extend(pool, length)
        assert 1 <= len(products) <= 3


def test_return_intermediate_traces_a_csv_row():
    row_index = next(
        index
        for index, row in _POLYOLEFIN.iterrows()
        if _free_slots(rxn_map[row["reaction_smarts"]]) == 1
    )
    row = _POLYOLEFIN.loc[row_index]
    config, pool, _mon1, _mon2, h1, _h2, _r1, _r2 = _row_data(row)
    reaction = ChainGrowthReaction(**{**config["params"], "max_products": 1})
    # Toggle the loaded reaction rather than rebuilding it or touching rxn_map.
    reaction.return_intermediate = True
    result = reaction.extend(pool, 3)

    assert len(result["intermediates"]) == 4  # n = 0..3
    assert all(len(step) == 1 for step in result["intermediates"])
    assert len(result["final"]) == 1

    heavy_counts = []
    for step in result["intermediates"]:
        mol = _mol(step[0])
        assert any(atom.GetAtomicNum() == 0 for atom in mol.GetAtoms())
        heavy_counts.append(mol.GetNumHeavyAtoms())

    # Initiation is a hydrogen cap; each propagation adds one monomer.
    assert heavy_counts[0] == h1
    for previous, current in pairwise(heavy_counts):
        assert current - previous == h1

    final = _mol(result["final"][0])
    assert not any(atom.GetAtomicNum() == 0 for atom in final.GetAtoms())
    # Termination removes the dummy and leaves a hydrogen end; it adds no
    # heavy atom.
    assert final.GetNumHeavyAtoms() == heavy_counts[-1]
