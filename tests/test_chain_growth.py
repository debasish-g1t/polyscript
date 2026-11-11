"""Tests for the chain-growth reaction driver.

The SMARTS here are deliberately simple hydrocarbon stand-ins so the chain
mechanics can be checked independently of the production rules in
:mod:`polyscript.extender.map`.
"""

import pytest

from polyscript.extender import ChainGrowthReaction

# Ethylene + a methyl initiator, head-to-tail propagation, methyl termination.
# The wildcard marks the growing end: init introduces it, propagation moves it,
# termination consumes it.
INIT = "[C:1]=[C:2]>>C[C:1][C:2]-*"
PROP = "[C:1][C:2]-*.[C:3]=[C:4]>>[C:1][C:2][C:3][C:4]-*"
TERM = "[C:1][C:2]-*>>[C:1][C:2]-C"


@pytest.fixture
def reaction():
    return ChainGrowthReaction(INIT, PROP, TERM)


def test_extend_counts_propagation_steps(reaction):
    # init: 3 C, each propagation adds 2 C, termination adds 1 C.
    assert reaction.extend(["C=C"], 0) == ["CCCC"]
    assert reaction.extend(["C=C"], 1) == ["CCCCCC"]
    assert reaction.extend(["C=C"], 3) == ["CCCCCCCCCC"]


def test_wildcard_reactant_matches_only_growing_end(reaction):
    """A bare ``*`` on the reactant side must not match backbone atoms."""
    chains = reaction.chain_init(["C=C"])
    assert [len(chain.GetAtoms()) for chain in chains] == [4]  # 3 C + dummy

    # If the reactant wildcard matched ordinary bonds, propagation would
    # branch and retain the old dummy; neither happens here.
    extended = reaction.chain_extend(["C=C"])
    assert [len(chain.GetAtoms()) for chain in extended] == [6]  # 5 C + dummy


def test_stateful_chain_methods(reaction):
    assert [len(chain.GetAtoms()) for chain in reaction.chain_init(["C=C"])] == [4]
    assert [len(chain.GetAtoms()) for chain in reaction.chain_extend(["C=C"])] == [6]
    terminated = reaction.chain_term(["C=C"])
    assert [len(chain.GetAtoms()) for chain in terminated] == [6]  # 6 C, no dummy


def test_propagation_can_add_two_units():
    prop_dimer = (
        "[C:1][C:2]-*.[C:3]=[C:4].[C:5]=[C:6]"
        ">>[C:1][C:2][C:3][C:4][C:5][C:6]-*"
    )
    reaction = ChainGrowthReaction(INIT, prop_dimer, TERM)
    # 3 C (init) + 2 steps * 4 C + 1 C (term) = 12 C.
    assert reaction.extend(["C=C"], 2) == ["CCCCCCCCCCCC"]


def test_monomers_may_be_mols(reaction):
    from rdkit import Chem

    assert reaction.extend([Chem.MolFromSmiles("C=C")], 1) == ["CCCCCC"]


def test_invalid_length(reaction):
    with pytest.raises(ValueError):
        reaction.extend(["C=C"], -1)


def test_initiation_rejects_wildcard_reactant():
    reaction = ChainGrowthReaction("*[C:1]>>[C:1]", PROP, TERM)
    with pytest.raises(ValueError, match="Initiation SMARTS"):
        reaction.chain_init(["C=C"])


def test_max_products_defaults_to_one():
    capped = ChainGrowthReaction(INIT, PROP, TERM)
    exhaustive = ChainGrowthReaction(INIT, PROP, TERM, max_products=None)

    # Propene has two propagation orientations, so exhaustive enumeration branches.
    assert len(exhaustive.extend(["CC=C"], 3)) > 1
    assert len(capped.extend(["CC=C"], 3)) == 1


def test_max_products_caps_every_frontier():
    reaction = ChainGrowthReaction(INIT, PROP, TERM, max_products=2)

    assert len(reaction.chain_init(["CC=C"])) <= 2
    assert len(reaction.chain_extend(["CC=C"])) <= 2
    assert len(reaction.chain_term([])) <= 2
    assert len(reaction.extend(["CC=C"], 5)) <= 2


def test_max_products_must_be_positive():
    with pytest.raises(ValueError, match="max_products"):
        ChainGrowthReaction(INIT, PROP, TERM, max_products=0)


def test_capped_result_is_deterministic():
    first = ChainGrowthReaction(INIT, PROP, TERM).extend(["CC=C"], 4)
    second = ChainGrowthReaction(INIT, PROP, TERM).extend(["CC=C"], 4)
    assert first == second


def test_return_intermediate_defaults_to_false(reaction):
    result = reaction.extend(["C=C"], 2)
    assert isinstance(result, list)
    assert result == ["CCCCCCCC"]


def test_return_intermediate_traces_every_step(reaction):
    from rdkit import Chem

    result = reaction.extend(["C=C"], 2, return_intermediate=True)
    assert set(result) == {"intermediates", "final"}
    assert len(result["intermediates"]) == 3  # init + two propagation steps
    assert all(len(step) == 1 for step in result["intermediates"])

    heavy_counts = []
    for step in result["intermediates"]:
        mol = Chem.MolFromSmiles(step[0])
        assert mol is not None
        assert any(atom.GetAtomicNum() == 0 for atom in mol.GetAtoms())
        heavy_counts.append(mol.GetNumHeavyAtoms())
    assert heavy_counts == [3, 5, 7]

    final = Chem.MolFromSmiles(result["final"][0])
    assert final is not None
    assert not any(atom.GetAtomicNum() == 0 for atom in final.GetAtoms())
    assert final.GetNumHeavyAtoms() == 8


def test_return_intermediate_at_length_zero(reaction):
    result = reaction.extend(["C=C"], 0, return_intermediate=True)
    assert len(result["intermediates"]) == 1
    assert len(result["intermediates"][0]) == 1
    assert result["final"] == ["CCCC"]


def test_return_intermediate_respects_frontier_cap():
    reaction = ChainGrowthReaction(INIT, PROP, TERM, max_products=2)
    result = reaction.extend(["CC=C"], 2, return_intermediate=True)
    assert all(len(step) <= 2 for step in result["intermediates"])
    assert len(result["final"]) <= 2


def test_return_intermediate_setter_toggles_trace(reaction):
    assert reaction.return_intermediate is False
    assert isinstance(reaction.extend(["C=C"], 2), list)

    reaction.return_intermediate = True
    assert reaction.return_intermediate is True
    assert isinstance(reaction.extend(["C=C"], 2), dict)


def test_return_intermediate_per_call_override(reaction):
    reaction.return_intermediate = True
    assert isinstance(reaction.extend(["C=C"], 2), dict)
    assert isinstance(
        reaction.extend(["C=C"], 2, return_intermediate=False), list
    )

    reaction.return_intermediate = False
    assert isinstance(
        reaction.extend(["C=C"], 2, return_intermediate=True), dict
    )


def test_return_intermediate_setter_rejects_non_bool(reaction):
    with pytest.raises(TypeError, match="return_intermediate"):
        reaction.return_intermediate = "yes"


def test_a1_scope_excludes_cyclic_alkenes_and_cumulenes():
    from polyscript.extender.map import rxn_map

    key = "[C&X3;H2,H1,H0:1]=[C;H2,H1,H0:2]>>*-[C:1][C:2]-*"
    reaction = ChainGrowthReaction(**rxn_map[key]["params"])

    # A terminal vinyl and an exocyclic vinyl on a ring are both acyclic
    # double bonds and must propagate.
    assert reaction.extend(["C=CC"], 1)
    assert reaction.extend(["C=C1CC(C)(c2cnc[nH]2)C1"], 1)

    # Endocyclic alkenes (norbornene) and cumulenes (allene) must not.
    assert reaction.extend(["C1CC2C=CC1C2"], 1) == []
    assert reaction.extend(["C=C=C"], 1) == []
