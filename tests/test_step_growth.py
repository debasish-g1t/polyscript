"""Tests for the stepwise step-growth driver.

The driver consumes one functional group from the chain and one from a
monomer.  The monomer's other groups are carried over untouched, so these
tests focus on the real terminal chemistry produced by the monomers rather
than on wildcard termination.
"""

from pathlib import Path

import pandas as pd
import pytest
from rdkit import Chem, RDLogger

from polyscript.extender import (
    RingOpeningMetathesisReaction,
    RingOpeningReaction,
    StepGrowthReaction,
)
from polyscript.extender.map import (
    RING_OPENING_PROPAGATIONS,
    STEP_GROWTH_COUPLINGS,
    reaction_config,
    rxn_map,
)

RDLogger.DisableLog("rdApp.*")

_CSV_PATH = (
    Path(__file__).resolve().parent
    / "test_data"
    / "extender_test_data_all.csv"
)
_DF = pd.read_csv(_CSV_PATH)


def _mol(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None, f"unparsable SMILES: {smiles}"
    return mol


def _count(mol, smarts: str) -> int:
    return len(mol.GetSubstructMatches(Chem.MolFromSmarts(smarts)))


def _find_type_key(cls: str, rxn_type: str, *needles: str, exclude=()) -> str:
    matches = [
        key
        for key, config in rxn_map.items()
        if config
        and config.get("type") == rxn_type
        and config["polymer_class"] == cls
        and all(needle in key for needle in needles)
        and not any(item in key for item in exclude)
    ]
    assert len(matches) == 1, (cls, rxn_type, needles, matches)
    return matches[0]


def _find_step_key(cls: str, *needles: str, exclude=()) -> str:
    return _find_type_key(cls, "step_growth", *needles, exclude=exclude)


def _find_ring_opening_key(
    cls: str, *needles: str, exclude=()
) -> str:
    return _find_type_key(
        cls, "ring_opening", *needles, exclude=exclude
    )


def _reaction(key: str, **kwargs):
    config = rxn_map[key]
    params = reaction_config(key)
    if config["type"] == "ring_opening":
        return RingOpeningReaction(**params, **kwargs)
    if config["type"] == "romp":
        return RingOpeningMetathesisReaction(**params, **kwargs)
    return StepGrowthReaction(**params, **kwargs)


POLYESTER_AABB = _find_step_key(
    "polyester", "[C&X3:1](=O)[O&X2&H1,Cl,Br]"
)
POLYESTER_AB = _find_step_key(
    "polyester", ">>(*-[O&X2:1].[C&X3:2](=O)-*)"
)
POLYESTER_CARBONATE = _find_step_key("polyester", "[C&-]#[O&+]")
POLYAMIDE_AB = _find_step_key("polyamide", "!$(OC=*)")
POLYURETHANE = _find_step_key(
    "polyurethane", "[N&X2:1]=[C&X2]=[O&X1,S&X1:2]"
)
PEEK = _find_step_key("polyether", "(F)", exclude=("[S&X4]",))
POLYETHER_EPOXIDE = _find_step_key(
    "polyether", "[C&X4;H2,H1,H0;R:1]1[O&R]"
)
POLYESTER_LACTONE_ROP = _find_type_key(
    "polyester", "ring_opening", "[C&X3&R:1](=[O&X1])[O&X2&R:2]"
)
POLYAMIDE_LACTAM_ROP = _find_type_key(
    "polyamide", "ring_opening", "[C&X3&R:1](=[O&X1])[N&X3&R:2]"
)
ROCOP = _find_type_key("polyester", "ring_opening", "[C,c;R:1]")
PHENOL_OXIDATIVE = _find_step_key("polyether", "c1([O&H1:1])")
ROMP = next(
    key
    for key, config in rxn_map.items()
    if config and config.get("type") == "romp"
)
POLYIMIDE = _find_step_key("polyimide", "[C&X3,c;R:1]")
POLYOXAZOLIDONE = _find_step_key(
    "polyoxazolidone", "[C&X4;H2,H1,H0;R:1]1[O&X2&R:2]"
)

ESTER_AABB_MONOMERS = ["OC(=O)CC(=O)O", "OCCO"]
ESTER_AB_MONOMERS = ["OCC(=O)O"]
AMIDE_MONOMERS = ["NCC(=O)O"]
URETHANE_MONOMERS = ["O=C=Nc1ccc(N=C=O)cc1", "OCCO"]
PEEK_MONOMERS = [
    "O=C(c1ccc(F)cc1)c1ccc(F)cc1",
    "Oc1ccc(O)cc1",
]
EPOXIDE_MONOMERS = ["CC1CO1"]
IMIDE_MONOMERS = [
    "O=C1OC(=O)c2cc3c(=O)oc(=O)c3cc12",
    "Nc1ccc(N)cc1",
]
OXAZOLIDONE_MONOMERS = [
    "O=C1CCCCCCC12OC2C12OC1C1CCCC12",
    "O=C=Nc1cccc(-c2ccc(N=C=O)c(Cl)c2)c1",
]
CARBONATE_MONOMERS = ["Oc1ccc(O)cc1"]


def test_length_one_returns_the_monomer():
    reaction = _reaction(POLYESTER_AB)
    assert reaction.extend(ESTER_AB_MONOMERS, 1) == ["O=C(O)CO"]


def test_polyester_preserves_acid_and_alcohol_ends():
    product = _mol(
        _reaction(POLYESTER_AABB).extend(ESTER_AABB_MONOMERS, 2)[0]
    )
    assert not any(atom.GetAtomicNum() == 0 for atom in product.GetAtoms())
    assert _count(product, "[CX3](=O)[OX2H1]") == 1
    assert _count(product, "[OX2H1;!$(OC=*)]") == 1
    assert _count(product, "[OX2;!H1][CX3](=O)") == 1


def test_polyester_preserves_acid_chloride_end():
    product = _mol(
        _reaction(POLYESTER_AABB).extend(
            ["O=C(Cl)CC(=O)Cl", "OCCO"], 2
        )[0]
    )
    # One acid chloride reacted; the other must remain -COCl, not become
    # -COOH through a generic cap.
    assert _count(product, "[CX3](=O)[Cl]") == 1
    assert _count(product, "[CX3](=O)[OX2H1]") == 0
    assert _count(product, "[OX2H1]") == 1


def test_polyamide_preserves_amine_and_acid_ends():
    product = _mol(_reaction(POLYAMIDE_AB).extend(AMIDE_MONOMERS, 3)[0])
    assert _count(product, "[NX3;H2]") == 1
    assert _count(product, "[CX3](=O)[OX2H1]") == 1


def test_polyurethane_preserves_isocyanate_and_alcohol_ends():
    product = _mol(
        _reaction(POLYURETHANE).extend(URETHANE_MONOMERS, 2)[0]
    )
    assert _count(product, "[NX2]=[CX2]=[OX1]") == 1
    assert _count(product, "[OX2H1]") == 1


def test_polyurethane_preserves_isothiocyanate_end():
    product = _mol(
        _reaction(POLYURETHANE).extend(
            ["S=C=Nc1ccc(N=C=S)cc1", "OCCO"], 2
        )[0]
    )
    assert _count(product, "[NX2]=[CX2]=[SX1]") == 1
    assert _count(product, "[NX2]=[CX2]=[OX1]") == 0
    assert _count(product, "[OX2H1]") == 1


def test_peek_preserves_aryl_fluoride_and_phenol_ends():
    product = _mol(_reaction(PEEK).extend(PEEK_MONOMERS, 2)[0])
    assert _count(product, "[c][F]") == 1
    assert _count(product, "[OX2H1]") == 1


def test_polyimide_preserves_amine_and_anhydride_ends():
    product = _mol(_reaction(POLYIMIDE).extend(IMIDE_MONOMERS, 2)[0])
    assert _count(product, "[NX3;H2]") == 1
    assert _count(
        product, "[#6&X3](=[#8&X1])[#8&X2][#6&X3](=[#8&X1])"
    ) == 1


def test_polyoxazolidone_preserves_isocyanate_and_epoxide_ends():
    product = _mol(
        _reaction(POLYOXAZOLIDONE).extend(OXAZOLIDONE_MONOMERS, 2)[0]
    )
    assert _count(product, "[NX2]=[CX2]=[OX1]") == 1
    assert _count(product, "[C&X4;R]1[O&X2&R][C&X4&R]1") == 1


def test_carbonate_preserves_phenol_ends():
    product = _mol(
        _reaction(POLYESTER_CARBONATE).extend(CARBONATE_MONOMERS, 2)[0]
    )
    assert _count(product, "[OX2][CX3](=O)[OX2]") == 1
    assert _count(product, "[OX2H1]") == 2


def test_epoxide_rop_preserves_hydroxyl_ends():
    product = _mol(
        _reaction(POLYETHER_EPOXIDE).extend(EPOXIDE_MONOMERS, 3)[0]
    )
    assert _count(product, "[OX2H1]") == 2
    assert _count(product, "[OX2]([CX4])[CX4]") == 1


def test_lactone_opens_to_hydroxyacid_and_esterifies():
    product = _mol(
        _reaction(POLYESTER_LACTONE_ROP).extend(
            ["O=C1CCCCCO1"], 2
        )[0]
    )
    assert _count(product, "[CX3](=O)[OX2H1]") == 1
    assert _count(product, "[OX2H1;!$(OC=*)]") == 1


def test_lactam_opens_to_aminoacid_and_amidates():
    product = _mol(
        _reaction(POLYAMIDE_LACTAM_ROP).extend(
            ["O=C1CCCCCN1"], 2
        )[0]
    )
    assert _count(product, "[NX3;H2]") == 1
    assert _count(product, "[CX3](=O)[OX2H1]") == 1


def test_thioimidate_lactam_propagates_past_the_generic_amide_guard():
    monomer = "O=C1CC(c2cccc3c2NC(=O)CC3=O)N=C(S)N1"
    products = _reaction(POLYAMIDE_LACTAM_ROP).extend([monomer], 3)
    assert len(products) == 1
    product = _mol(products[0])
    assert not any(atom.GetAtomicNum() == 0 for atom in product.GetAtoms())
    # The amidine/thioimidate amine is still a primary amine; the fused
    # imide N-H is not consumed as a nucleophile.
    assert _count(product, "[NX3;H2]") == 1
    assert _count(product, "[CX3](=O)[OX2H1]") == 1


def test_cyclic_carbamate_lactone_propagates():
    monomer = "CC1CN(C(=O)CCB2OC(C)(C)C(C)(C)O2)C(=O)OC1=O"
    products = _reaction(POLYESTER_LACTONE_ROP).extend([monomer], 3)
    assert len(products) == 1
    product = _mol(products[0])
    assert not any(atom.GetAtomicNum() == 0 for atom in product.GetAtoms())
    # The opened unit retains one real carboxylic acid end; the acyloxy O-H
    # (carbamic acid) is what propagates, and it is still a C(=O)-O-H group.
    assert _count(product, "[CX3](=O)[OX2H1]") == 2
    assert _count(product, "[NX3][CX3](=O)[OX2H1]") == 1


def test_anhydride_epoxide_rocop_inserts_alternating_pairs():
    rows = _DF[_DF["reaction_smarts"] == ROCOP]
    row = rows.iloc[0]
    reaction = _reaction(ROCOP)
    monomer_oxygen = _count(
        _mol(row["mon1"]), "[#8]"
    ) + _count(_mol(row["mon2"]), "[#8]")
    counts = []
    oxygen = []
    for length in (1, 2, 3):
        products = reaction.extend([row["mon1"], row["mon2"]], length)
        assert len(products) == 1
        product = _mol(products[0])
        assert not any(
            atom.GetAtomicNum() == 0 for atom in product.GetAtoms()
        )
        counts.append(product.GetNumHeavyAtoms())
        oxygen.append(_count(product, "[#8]"))

    # The balanced opening consumes one water molecule.  Its O is the extra
    # oxygen in residue 1; each later residue adds only the monomers' own
    # oxygens.
    assert oxygen[0] == monomer_oxygen + 1
    assert oxygen[1] == monomer_oxygen * 2 + 1
    assert oxygen[2] == monomer_oxygen * 3 + 1
    assert counts[1] > counts[0] > 0
    assert counts[2] > counts[1]


def _label_heavy_atoms(mol, base: int):
    rw = Chem.RWMol(mol)
    for index, atom in enumerate(rw.GetAtoms()):
        if atom.GetAtomicNum() > 1:
            atom.SetIsotope(base + index)
    return rw.GetMol()


def test_rocop_opening_consumes_water_without_inventing_atoms():
    """The first residue contains the water O, not an invented oxygen.

    Label every heavy atom in the monomers and the opening water.  The
    balanced opening leaves no unlabelled heavy atom in the opened unit; the
    previous unbalanced key left several.  This is the provenance check the
    benchmark cannot perform because it labels the pool after ring opening.
    """
    rows = _DF[_DF["reaction_smarts"] == ROCOP]
    row = rows.iloc[0]
    reaction = _reaction(ROCOP)
    reaction.opening_reagents = [
        _label_heavy_atoms(_mol("O"), 9000)
    ]
    pool = reaction.prepare_pool(
        [
            _label_heavy_atoms(_mol(row["mon1"]), 1000),
            _label_heavy_atoms(_mol(row["mon2"]), 2000),
        ]
    )
    assert pool
    monomer_oxygen = _count(
        _mol(row["mon1"]), "[#8]"
    ) + _count(_mol(row["mon2"]), "[#8]")
    for mol in pool:
        unlabelled = [
            atom.GetSymbol()
            for atom in mol.GetAtoms()
            if atom.GetAtomicNum() > 1 and atom.GetIsotope() == 0
        ]
        assert not unlabelled
        assert _count(mol, "[#8]") == monomer_oxygen + 1


def test_ring_opening_trace_has_complete_intermediates():
    reaction = _reaction(POLYESTER_LACTONE_ROP)
    result = reaction.extend(
        ["O=C1CCCCCO1"], 3, return_intermediate=True
    )
    assert len(result["intermediates"]) == 3
    for step in result["intermediates"]:
        assert "*" not in step[0]
    assert len(result["final"]) == 1
    assert "*" not in result["final"][0]


def test_phenol_oxidative_couples_phenols_through_the_para_position():
    rows = _DF[_DF["reaction_smarts"] == PHENOL_OXIDATIVE]
    row = rows.iloc[0]
    product = _mol(
        _reaction(PHENOL_OXIDATIVE).extend([row["mon1"]], 2)[0]
    )
    assert _count(product, "[c][OX2][c]") == 1
    assert _count(product, "[c][OX2H1]") == 1


def test_romp_inserts_cycloalkenes_without_wildcards():
    reaction = _reaction(ROMP)
    for monomer in ("C1=CCCCC1", "C1CC2C=CC1C2"):
        heavy = []
        double_bonds = []
        rings = []
        original = _mol(monomer)
        unit = original.GetNumHeavyAtoms()
        unit_rings = original.GetRingInfo().NumRings()
        for length in (0, 1, 2, 3):
            products = reaction.extend([monomer], length)
            assert len(products) == 1
            product = _mol(products[0])
            assert not any(
                atom.GetAtomicNum() == 0 for atom in product.GetAtoms()
            )
            heavy.append(product.GetNumHeavyAtoms())
            double_bonds.append(_count(product, "C=C"))
            rings.append(product.GetRingInfo().NumRings())
        assert heavy[0] == unit + 2  # initiator and terminator fragments
        for previous, current in zip(heavy, heavy[1:]):
            assert current - previous == unit
        # Metathesis cleaves the ring C=C and joins units through new C=C
        # bonds: n inserted rings give n + 1 backbone/end double bonds.
        assert double_bonds == [length + 2 for length in (0, 1, 2, 3)]
        # Each opened ring loses its alkene ring, so a norbornene unit keeps
        # one ring and a cyclohexene unit keeps none.
        assert rings == [
            (length + 1) * (unit_rings - 1)
            for length in (0, 1, 2, 3)
        ]


def test_ab_length_counts_monomer_residues():
    reaction = _reaction(POLYESTER_AB)
    for length in (1, 2, 3):
        product = _mol(reaction.extend(ESTER_AB_MONOMERS, length)[0])
        assert len(
            product.GetSubstructMatches(Chem.MolFromSmarts("[CX4;H2]"))
        ) == length


def test_aabb_length_counts_both_residues():
    reaction = _reaction(POLYESTER_AABB)
    product = _mol(reaction.extend(ESTER_AABB_MONOMERS, 2)[0])

    # Malonic acid contributes two carbonyls; ethylene glycol contributes
    # one CH2-CH2 pair.
    assert len(
        product.GetSubstructMatches(Chem.MolFromSmarts("[#6&X3]=[#8&X1]"))
    ) == 2
    assert len(
        product.GetSubstructMatches(Chem.MolFromSmarts("[CH2][CH2]"))
    ) == 1


def test_trace_returns_intermediates_without_wildcards():
    reaction = _reaction(POLYESTER_AB)
    result = reaction.extend(
        ESTER_AB_MONOMERS, 3, return_intermediate=True
    )

    assert set(result) == {"intermediates", "final"}
    assert len(result["intermediates"]) == 3
    for step in result["intermediates"]:
        assert len(step) == 1
        assert "*" not in step[0]

    assert len(result["final"]) == 1
    assert "*" not in result["final"][0]


def test_return_intermediate_setter_toggles_trace():
    reaction = _reaction(POLYESTER_AB)
    assert isinstance(reaction.extend(ESTER_AB_MONOMERS, 2), list)

    reaction.return_intermediate = True
    assert isinstance(reaction.extend(ESTER_AB_MONOMERS, 2), dict)

    reaction.return_intermediate = False
    assert isinstance(reaction.extend(ESTER_AB_MONOMERS, 2), list)


def test_length_zero_returns_empty_chain():
    reaction = _reaction(POLYESTER_AB)
    assert reaction.extend(ESTER_AB_MONOMERS, 0) == []
    assert reaction.extend(
        ESTER_AB_MONOMERS, 0, return_intermediate=True
    ) == {"intermediates": [], "final": []}


def test_invalid_length():
    reaction = _reaction(POLYESTER_AB)
    with pytest.raises(ValueError, match="length"):
        reaction.extend(ESTER_AB_MONOMERS, -1)
    with pytest.raises(ValueError, match="length"):
        reaction.extend(ESTER_AB_MONOMERS, 1.5)


def test_max_products_must_be_positive():
    with pytest.raises(ValueError, match="max_products"):
        StepGrowthReaction(
            coupling_smarts=STEP_GROWTH_COUPLINGS["ester"],
            max_products=0,
        )


def test_max_products_caps_the_frontier():
    reaction = _reaction(POLYESTER_AB, max_products=1)
    assert len(reaction.extend(ESTER_AB_MONOMERS, 3)) == 1


def test_coupling_smarts_compile():
    from rdkit.Chem import AllChem

    for name, smarts in {
        **STEP_GROWTH_COUPLINGS,
        **RING_OPENING_PROPAGATIONS,
    }.items():
        alternatives = smarts if isinstance(smarts, list) else [smarts]
        for alternative in alternatives:
            reaction = AllChem.ReactionFromSmarts(alternative)
            assert reaction is not None, (name, alternative)


def test_amide_guard_distinguishes_amidine_from_amide_and_sulfonamide():
    from rdkit.Chem import AllChem

    reaction = AllChem.ReactionFromSmarts(STEP_GROWTH_COUPLINGS["amide"])
    acid = _mol("CC(=O)O")
    amidine = _mol("NC(=N)N")
    amide = _mol("CC(=O)N")
    sulfonamide = _mol("CS(=O)(=O)N")
    phosphoramide = _mol("CP(=O)(N)N")

    # The amidine N-H is a free amine; the amide, sulfonamide, and
    # phosphoramide N-H are not.
    assert reaction.RunReactants((acid, amidine))
    assert not reaction.RunReactants((acid, amide))
    assert not reaction.RunReactants((acid, sulfonamide))
    assert not reaction.RunReactants((acid, phosphoramide))


def test_ester_guard_rejects_heteroatom_acids_and_peroxides():
    from rdkit.Chem import AllChem

    reaction = AllChem.ReactionFromSmarts(
        STEP_GROWTH_COUPLINGS["ester"][0]
    )
    acid = _mol("CC(=O)O")

    # Real alcohol and phenol nucleophiles react.
    assert reaction.RunReactants((acid, _mol("CCO")))
    assert reaction.RunReactants((acid, _mol("Oc1ccccc1")))

    # Heteroatom acids, oximes, and peroxides are not alcohols.
    for smiles in (
        "OB(O)c1ccccc1",
        "OS(=O)(=O)c1ccccc1",
        "CC=NO",
        "COO",
    ):
        assert not reaction.RunReactants((acid, _mol(smiles))), smiles


def test_ester_transesterification_accepts_methyl_esters_only():
    from rdkit.Chem import AllChem

    reaction = AllChem.ReactionFromSmarts(
        STEP_GROWTH_COUPLINGS["ester"][1]
    )
    diol = _mol("OCCO")

    # DMT (dimethyl terephthalate) is the industrial transesterification
    # monomer; its methyl esters react with the diol.
    dmt = _mol("COC(=O)c1ccc(C(=O)OC)cc1")
    assert reaction.RunReactants((dmt, diol))

    # A backbone-like ethyl/aryl ester is deliberately not matched: the
    # methyl-specific guard prevents the chain from cleaving its own esters.
    backbone = _mol("CCOC(=O)c1ccccc1")
    assert not reaction.RunReactants((backbone, diol))


def test_lactone_rop_allows_a_carbamic_acid_oh_only():
    from rdkit.Chem import AllChem

    standard, carbamic = [
        AllChem.ReactionFromSmarts(smarts)
        for smarts in RING_OPENING_PROPAGATIONS["lactone_rop"]
    ]
    acid = _mol("CC(=O)O")
    carbamic_acid = _mol("CC(=O)NC(=O)O")

    # The generic ester coupling rejects the carbamic acid O-H; the dedicated
    # ROP alternative accepts exactly that end group.
    assert not standard.RunReactants((acid, carbamic_acid))
    assert carbamic.RunReactants((acid, carbamic_acid))


def test_imide_requires_a_cyclic_anhydride():
    from rdkit.Chem import AllChem

    reaction = AllChem.ReactionFromSmarts(STEP_GROWTH_COUPLINGS["imide"])
    amine = _mol("Nc1ccc(N)cc1")

    # Acyclic acetic anhydride must not be scrambled into an imide.
    assert not reaction.RunReactants((_mol("CC(=O)OC(=O)C"), amine))
    # Succinic (five-membered) and glutaric (six-membered) rings react.
    assert reaction.RunReactants((_mol("O=C1OC(=O)CC1"), amine))
    assert reaction.RunReactants((_mol("O=C1OC(=O)CCC1"), amine))


def test_imide_accepts_aromatic_fused_and_biphenyl_anhydrides():
    from rdkit.Chem import AllChem

    reaction = AllChem.ReactionFromSmarts(STEP_GROWTH_COUPLINGS["imide"])
    amine = _mol("Nc1ccc(N)cc1")
    pmda = _mol("O=C1OC(=O)c2cc3c(=O)oc(=O)c3cc12")
    bpda = _mol("O=C1OC(=O)c2cc(-c3ccc4c(c3)C(=O)OC4=O)ccc21")

    # PMDA is the fully aromatic fused anhydride; BPDA is the biphenyl
    # dianhydride.  Both must form the imide without a kekulisation error.
    assert reaction.RunReactants((pmda, amine))
    assert reaction.RunReactants((bpda, amine))


def test_oxazolidinone_attacks_the_less_hindered_epoxide_carbon():
    from rdkit.Chem import AllChem

    reaction = AllChem.ReactionFromSmarts(
        STEP_GROWTH_COUPLINGS["oxazolidinone"][0]
    )
    products = {
        Chem.MolToSmiles(result[0])
        for result in reaction.RunReactants(
            (_mol("CC1CO1"), _mol("O=C=Nc1ccccc1"))
        )
    }
    assert products == {"CC1CN(c2ccccc2)C(=O)O1"}

    product = _mol(next(iter(products)))
    # The methyl stays on the O-bearing (more substituted) carbon: a
    # 5-substituted oxazolidin-2-one, not the 4-substituted regioisomer.
    assert _count(
        product, "[CH3][C&X4]1[O&X2][C&X3](=[O])[N&X3][C&X4]1"
    ) == 1
    assert _count(
        product, "[CH3][C&X4]1[N&X3][C&X3](=[O])[O&X2][C&X4]1"
    ) == 0


def test_epoxide_rop_attacks_the_less_hindered_carbon():
    from rdkit.Chem import AllChem

    reaction = AllChem.ReactionFromSmarts(
        STEP_GROWTH_COUPLINGS["epoxide_rop"][0]
    )
    products = {
        Chem.MolToSmiles(result[0])
        for result in reaction.RunReactants(
            (_mol("CCO"), _mol("CC1CO1"))
        )
    }
    assert products == {"CCOCC(C)O"}


def test_aryl_ether_nas_requires_para_activation():
    from rdkit.Chem import AllChem

    ketone = AllChem.ReactionFromSmarts(
        STEP_GROWTH_COUPLINGS["aryl_ether_nas"][0]
    )
    sulfone = AllChem.ReactionFromSmarts(
        STEP_GROWTH_COUPLINGS["aryl_ether_nas"][11]
    )
    nitrile = AllChem.ReactionFromSmarts(
        STEP_GROWTH_COUPLINGS["aryl_ether_nas"][2]
    )
    ortho_nitrile = AllChem.ReactionFromSmarts(
        STEP_GROWTH_COUPLINGS["aryl_ether_nas"][3]
    )
    biphenyl_ketone = AllChem.ReactionFromSmarts(
        STEP_GROWTH_COUPLINGS["aryl_ether_nas"][4]
    )
    phenol = _mol("Oc1ccc(O)cc1")

    # Fluoride templates are tried first, matching SNAr reactivity.
    assert "([F])" in STEP_GROWTH_COUPLINGS["aryl_ether_nas"][0]
    assert "([Cl,Br,I])" in STEP_GROWTH_COUPLINGS["aryl_ether_nas"][10]

    # Unactivated aryl halides do not undergo uncatalysed SNAr.
    assert not ketone.RunReactants((_mol("Fc1ccccc1"), phenol))
    assert not ketone.RunReactants((_mol("Clc1ccc(Cl)cc1"), phenol))

    # The para ketone (PEEK) and para sulfone (PES) monomers do.
    assert ketone.RunReactants(
        (_mol("O=C(c1ccc(F)cc1)c1ccc(F)cc1"), phenol)
    )
    assert sulfone.RunReactants(
        (_mol("O=S(=O)(c1ccc(Cl)cc1)c1ccc(Cl)cc1"), phenol)
    )

    # Nitrile-activated (PAEN) and phenyl-spaced (biphenyl) aryl halides are
    # also para-activated NAS substrates.
    assert nitrile.RunReactants((_mol("N#Cc1ccc(F)cc1"), phenol))
    assert ortho_nitrile.RunReactants(
        (_mol("Fc1cccc(F)c1C#N"), phenol)
    )
    assert biphenyl_ketone.RunReactants(
        (
            _mol(
                "O=C(c1ccc(-c2ccc(F)cc2)cc1)c1ccc(-c2ccc(F)cc2)cc1"
            ),
            phenol,
        )
    )


def test_phenol_oxidative_requires_blocked_ortho_positions():
    product = _reaction(PHENOL_OXIDATIVE).extend(
        ["Cc1cccc(C)c1O"], 3
    )
    assert product
    with pytest.raises(ValueError):
        _reaction(PHENOL_OXIDATIVE).extend(["Oc1ccccc1"], 3)


def test_rocop_opens_heteroatom_anhydrides():
    reaction = _reaction(ROCOP)
    product = _mol(
        reaction.extend(["O=C1COC(=O)O1", "C1CO1"], 2)[0]
    )
    assert _count(product, "[CX3](=O)[OX2H1]") == 1
    assert _count(product, "[OX2][CX4][CX4][OX2H1]") == 1


def test_carbonate_uses_phosgene_not_direct_co_insertion():
    config = rxn_map[POLYESTER_CARBONATE]
    params = config["params"]
    assert params["reagents"] == ["ClC(=O)Cl"]
    assert all("[C-]#[O+]" not in smarts for smarts in params["coupling_smarts"])


def test_reaction_config_injects_the_map_key():
    assert reaction_config(POLYESTER_LACTONE_ROP)[
        "ring_open_smarts"
    ] == POLYESTER_LACTONE_ROP
    assert reaction_config(POLYAMIDE_LACTAM_ROP)[
        "ring_open_smarts"
    ] == POLYAMIDE_LACTAM_ROP
    rocop = reaction_config(ROCOP)
    assert "[O&X2;H2:9]" in rocop["ring_open_smarts"]
    assert rocop["opening_reagents"] == ["O"]
    assert reaction_config(ROMP)["ring_bond_smarts"] == ROMP


def test_configured_entries():
    by_type = {}
    for config in rxn_map.values():
        if config is None:
            continue
        by_type.setdefault(config["type"], []).append(config)
    # A1 has a second compatibility key for the recorded dataset SMARTS.
    assert len(by_type["addition"]) == 5
    assert len(by_type["step_growth"]) == 14
    assert len(by_type["ring_opening"]) == 3
    assert len(by_type["romp"]) == 1
    assert [key for key, config in rxn_map.items() if config is None] == []


def test_sample_data_couplings():
    step_like = {"step_growth", "ring_opening"}
    rows = _DF[_DF["reaction_smarts"].isin(
        key for key, config in rxn_map.items()
        if config and config.get("type") in step_like
    )]
    failures = []
    for reaction_smarts, group in rows.groupby("reaction_smarts"):
        reaction = _reaction(reaction_smarts)
        row = group.iloc[0]
        monomers = [row["mon1"]]
        if isinstance(row["mon2"], str):
            monomers.append(row["mon2"])
        try:
            products = reaction.extend(monomers, 2)
            assert products
            assert all("*" not in smiles for smiles in products)
        except Exception as exc:  # noqa: BLE001 - report the first row
            failures.append((reaction_smarts[:50], type(exc).__name__, str(exc)))
    assert not failures, failures
