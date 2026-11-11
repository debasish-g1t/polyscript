"""Tests for representation adapters (PSMILES, SMILES, SELFIES, BigSMILES)."""

import pytest
from polyscript.utils import (
    BigSMILESAdapter,
    PSMILESAdapter,
    SELFIESAdapter,
    SMILESAdapter,
)


# ---------------------------------------------------------------------------
# Test polymers
# ---------------------------------------------------------------------------

PS1 = (
    "{O=C(O)CCc1ccc(CCC(=O)O)nc1}"
    "+{Nc1ccc(C(c2cccc(B(O)O)c2)(C(F)(F)F)C(F)(F)F)cc1N}"
    "=>([C&X3:1](=O)[O&X2&H1,Cl,Br].[C&X3:2](=O)[O&X2&H1,Cl,Br])."
    "([N&X3;H2,H1;!$(NC=*):3].[N&X3;H2,H1;!$(NC=*):4])"
    ">>(*-[C&X3:1]=O.[C&X3:2](=O)-[N&X3&!$(NC=*):3].[N&X3&!$(NC=*):4]-*)"
    "=>*C(=O)C1C(C(=O)Nc2nc(N)nc(N(*)COc3cc(Cl)ccc3O)n2)C(C(=O)O)C1C(=O)O"
)

PS2 = (
    "{O=C1CN(C2C(Cl)C3CC2C2C(=O)OC(=O)C32)CC(=O)O1}"
    "+{O=C1CN(C2C(Cl)C3CC2C2C(=O)OC(=O)C32)CC(=O)O1}"
    "=>([C&X3,c;R:1](=[O&X1])[O&X2,o;R][C&X3,c;R:2]=[O&X1]."
    "[C&X3,c;R:3](=[O&X1])[O&X2,o;R][C&X3,c;R:4]=[O&X1])."
    "([C,c:5][N&X3&H2&!$(N[C,S]=*)].[C,c:6][N&X3&H2&!$(N[C,S]=*)])"
    ">>([C&X3,c;R:1](=[O&X1])[N&X3&R]([C,c:5])[C&X3,c;R:2]=[O&X1]."
    "[C,c:6]-*.[C&X3,c;R:3](=[O&X1])[N&X3&R](-*)[C&X3&R:4]=[O&X1])"
    "=>*c1c(CC#CC2CC2)ccc(N2C(=O)CN(C3C(Cl)C4CC3C3C(=O)N(*)C(=O)C43)CC2=O)c1"
)

PS3 = (
    "{CC(=O)O}+{NC}"
    "=>([C&X3:1](=O)[O&X2]).[N&X3&H2:2]"
    ">>(*-[C&X3:1]=O.[C&X3:1](=O)-[N&X3:2]-*)"
    "=>*C(=O)NC*"
)

POLYMERS = [
    pytest.param(PS1, id="complex-polyamide"),
    pytest.param(PS2, id="benchmark-polymer"),
    pytest.param(PS3, id="minimal-sanity"),
]


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def psmiles_adapter():
    return PSMILESAdapter()


@pytest.fixture
def smiles_adapter():
    return SMILESAdapter()


@pytest.fixture
def selfies_adapter():
    return SELFIESAdapter()


@pytest.fixture
def bigsmiles_adapter():
    return BigSMILESAdapter()


# ---------------------------------------------------------------------------
# polyscript → each format
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("polym", POLYMERS)
def test_polyscript_to_psmiles(psmiles_adapter, polym):
    """Extract PSMILES from a valid PolyScript string."""
    result = psmiles_adapter.convert(polym, _from="polyscript")
    assert isinstance(result, str)
    assert len(result) > 0


@pytest.mark.parametrize("polym", POLYMERS)
def test_polyscript_to_smiles(smiles_adapter, polym):
    """Extract dot-separated monomer SMILES from PolyScript."""
    result = smiles_adapter.convert(polym, _from="polyscript")
    assert isinstance(result, str)
    assert "." in result  # at least two monomers


@pytest.mark.parametrize("polym", POLYMERS)
def test_polyscript_to_selfies(selfies_adapter, polym):
    """Encode the polymer as SELFIES."""
    result = selfies_adapter.convert(polym, _from="polyscript")
    assert isinstance(result, str)
    assert result.startswith("[")  # SELFIES always starts with bracket


@pytest.mark.parametrize("polym", POLYMERS)
def test_polyscript_to_bigsmiles(bigsmiles_adapter, polym):
    """Convert polymer to BigSMILES representation."""
    result = bigsmiles_adapter.convert(polym, _from="polyscript")
    assert isinstance(result, str)
    assert len(result) > 0


# ---------------------------------------------------------------------------
# Cross-format conversions
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("polym", POLYMERS)
def test_psmiles_to_selfies(psmiles_adapter, selfies_adapter, polym):
    """PSMILES → SELFIES round-trip through encoder."""
    ps = psmiles_adapter.convert(polym, _from="polyscript")
    result = selfies_adapter.convert(ps, _from="psmiles")
    assert isinstance(result, str)
    assert result.startswith("[")


@pytest.mark.parametrize("polym", POLYMERS)
def test_psmiles_to_bigsmiles(psmiles_adapter, bigsmiles_adapter, polym):
    """PSMILES → BigSMILES works."""
    ps = psmiles_adapter.convert(polym, _from="polyscript")
    result = bigsmiles_adapter.convert(ps, _from="psmiles")
    assert isinstance(result, str)
    assert len(result) > 0


@pytest.mark.parametrize("polym", POLYMERS)
def test_selfies_to_psmiles(psmiles_adapter, selfies_adapter, polym):
    """Decode SELFIES back to PSMILES."""
    sf = selfies_adapter.convert(polym, _from="polyscript")
    result = psmiles_adapter.convert(sf, _from="selfies")
    assert isinstance(result, str)
    assert len(result) > 0


@pytest.mark.parametrize("polym", POLYMERS)
def test_selfies_to_smiles(smiles_adapter, selfies_adapter, polym):
    """Decode SELFIES back to SMILES."""
    sf = selfies_adapter.convert(polym, _from="polyscript")
    result = smiles_adapter.convert(sf, _from="selfies")
    assert isinstance(result, str)
    assert len(result) > 0


@pytest.mark.parametrize("polym", POLYMERS)
def test_smiles_to_selfies(smiles_adapter, selfies_adapter, polym):
    """Encode monomer SMILES as SELFIES."""
    smi = smiles_adapter.convert(polym, _from="polyscript")
    result = selfies_adapter.convert(smi, _from="smiles")
    assert isinstance(result, str)
    assert result.startswith("[")


# ---------------------------------------------------------------------------
# Edge cases: error handling
# ---------------------------------------------------------------------------


def test_invalid_polyscript_raises_valueerror(psmiles_adapter):
    """Malformed PolyScript input raises ValueError."""
    with pytest.raises(ValueError):
        psmiles_adapter.convert("not-a-valid-polyscript", _from="polyscript")


def test_unknown_from_raises_valueerror(smiles_adapter):
    """Unsupported _from keyword raises ValueError."""
    with pytest.raises(ValueError):
        smiles_adapter.convert("CCO", _from="xyz")


def test_bigsmiles_from_unsupported_source_raises(bigsmiles_adapter):
    """BigSMILES from an unsupported _from raises ValueError."""
    with pytest.raises(ValueError):
        bigsmiles_adapter.convert("CCO", _from="smiles")
