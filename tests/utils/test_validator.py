"""Tests for SeqValidator."""

import pytest
from polyscript.utils import SeqValidator


# ---------------------------------------------------------------------------
# Test candidates
# ---------------------------------------------------------------------------

VALID = (
    "{CNC(C)(C(=O)O)C(C)(C)C}+"
    "{CNC(CC(=O)O)C(=O)O}=>"
    "([N&X3;H2,H1;!$(NC=*):1].[C&X3:2](=O)[O&X2&H1])."
    "([N&X3;H2,H1;!$(NC=*):3].[C&X3:4](=O)[O&X2&H1])>>"
    "(*-[N&X3&!$(NC=*):1].[C&X3:2](=O)[N&X3&!$(NC=*):3].[C&X3:4](=O)-*)=>"
    "*C(=O)CC(C(=O)O)N(C)C(=O)C(C)(N(*)C)C(C)(C)C|||"
    "{*C(=O)CC(C(=O)O)N(C)C(=O)C(C)(N(*)C)C(C)(C)C}+{none}=>"
    "[C&X3:1](=O)[O&X2&H1,F,Cl,Br,I]>>[C&X3:1](=O)-*=>"
    "*C(=O)CC(C(*)=O)N(C)C(=O)C(C)(N(*)C)C(C)(C)C"
)

ONE_BAD_MONOMER = VALID.replace(
    "CNC(C)(C(=O)O)C(C)(C)C", "CNC_invalid_(C)(C(=O)O)C(C)(C)C"
)

TWO_BAD_MONOMERS = (
    ONE_BAD_MONOMER.replace(
        "CNC(CC(=O)O)C(=O)O", "CNC__invalid__(CC(=O)O)C(=O)O"
    )
)

BAD_REACTION = VALID.replace(
    "[N&X3;H2,H1;!$(NC=*):1]", "[N&X3;H2,H1;!$(NC=*):1]<invalid_rxn>"
)

BAD_POLYMER = VALID.replace(
    "*C(=O)CC(C(*)=O)N(C)C(=O)C(C)(N(*)C)C(C)(C)C",
    "*C(=O)CC(C(*)=O)N(C)C(=O)C(C)(N(*)__invalid__C)C(C)(C)C",
)

WRONG_POLYMER = VALID.replace(
    "*C(=O)CC(C(*)=O)N(C)C(=O)C(C)(N(*)C)C(C)(C)C",
    "*C(=O)CC(C(*)=O)N(C)C(=O)C(C)(N(*)C)C(C)(C)CC",
)


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


class TestSeqValidator:
    """Tests for sequence validation of PolyScript candidates."""

    @pytest.fixture
    def validator(self):
        return SeqValidator()

    def test_valid_sequence_no_errors(self, validator):
        """A fully valid nested candidate produces zero errors."""
        errors = validator.validate(VALID)
        assert errors == []

    def test_one_bad_monomer_reports_invalid_smiles(self, validator):
        """One invalid monomer SMILES triggers monomer + reaction errors."""
        errors = validator.validate(ONE_BAD_MONOMER)
        assert errors == [
            "<E-chem-|invalid-monomer-smiles|>",
            "<E-chem-|reaction-execution-error|>: reaction called with None "
            "reactants",
        ]

    def test_two_bad_monomers_reports_both(self, validator):
        """Two invalid monomer SMILES produce two invalid-monomer errors."""
        errors = validator.validate(TWO_BAD_MONOMERS)
        assert errors == [
            "<E-chem-|invalid-monomer-smiles|>",
            "<E-chem-|invalid-monomer-smiles|>",
            "<E-chem-|reaction-execution-error|>: reaction called with None "
            "reactants",
        ]

    def test_bad_reaction_reports_loading_error(self, validator):
        """A malformed reaction SMARTS triggers a reaction-loading error."""
        errors = validator.validate(BAD_REACTION)
        assert errors == [
            "<E-chem-|reaction-loading-error|>: "
            "ChemicalReactionParserException: multi-step reactions not supported"
        ]

    def test_bad_polymer_reports_validation_error_and_no_match(
        self, validator
    ):
        """An invalid polymer SMILES yields validation + no-match errors."""
        errors = validator.validate(BAD_POLYMER)
        assert errors == [
            "<E-chem-|polymer-validation-error|>: invalid polymer SMILES",
            "<E-chem-|no-polymer-match|>",
        ]

    def test_wrong_polymer_reports_no_match(self, validator):
        """A valid-but-wrong polymer SMILES yields only a no-match error."""
        errors = validator.validate(WRONG_POLYMER)
        assert errors == ["<E-chem-|no-polymer-match|>"]
