"""Tests for PolyScriptParser and NestedPolyScriptParser."""

import pytest
from polyscript.utils import NestedPolyScriptParser, PolyScriptParser


# ---------------------------------------------------------------------------
# Test data
# ---------------------------------------------------------------------------

VALID_SINGLE = (
    "{CNC(C)(C(=O)O)C(C)(C)C}+"
    "{CNC(C(=O)O)C(C)C}=>"
    "([N&X3;H2,H1;!$(NC=*):1].[C&X3:2](=O)[O&X2&H1])."
    "([N&X3;H2,H1;!$(NC=*):3].[C&X3:4](=O)[O&X2&H1])>>"
    "(*-[N&X3&!$(NC=*):1].[C&X3:2](=O)[N&X3&!$(NC=*):3].[C&X3:4](=O)-*)=>"
    "*C(=O)C(C(C)C)N(C)C(=O)C(C)(N(*)C)C(C)(C)C"
)

INVALID_SINGLE = (
    "{CNC(C)(C(=O)O)C(C)(C)C}+"
    "{CNC(C(=O)O)C(C)C}"
    "<_invalid_sign_>"
    "([N&X3;H2,H1;!$(NC=*):1].[C&X3:2](=O)[O&X2&H1])."
    "([N&X3;H2,H1;!$(NC=*):3].[C&X3:4](=O)[O&X2&H1])>>"
    "(*-[N&X3&!$(NC=*):1].[C&X3:2](=O)[N&X3&!$(NC=*):3].[C&X3:4](=O)-*)=>"
    "*C(=O)C(C(C)C)N(C)C(=O)C(C)(N(*)C)C(C)(C)C"
)

VALID_NESTED = (
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

INVALID_NESTED = (
    "{CNC(C)(C(=O)O)C(C)(C)C}+"
    "{CNC(CC(=O)O)C(=O)O}=>"
    "([N&X3;H2,H1;!$(NC=*):1].[C&X3:2](=O)[O&X2&H1])."
    "([N&X3;H2,H1;!$(NC=*):3].[C&X3:4](=O)[O&X2&H1])>>"
    "(*-[N&X3&!$(NC=*):1].[C&X3:2](=O)[N&X3&!$(NC=*):3].[C&X3:4](=O)-*)=>"
    "*C(=O)CC(C(=O)O)N(C)C(=O)C(C)(N(*)C)C(C)(C)C"
    "|<_invalid_sign_>|"
    "{*C(=O)CC(C(=O)O)N(C)C(=O)C(C)(N(*)C)C(C)(C)C}+{none}=>"
    "[C&X3:1](=O)[O&X2&H1,F,Cl,Br,I]>>[C&X3:1](=O)-*=>"
    "*C(=O)CC(C(*)=O)N(C)C(=O)C(C)(N(*)C)C(C)(C)C"
)


# ---------------------------------------------------------------------------
# PolyScriptParser (single sequence)
# ---------------------------------------------------------------------------


class TestPolyScriptParser:
    """Tests for single-sequence PolyScriptParser."""

    @pytest.fixture
    def parser(self):
        return PolyScriptParser()

    def test_parse_valid_single_sequence(self, parser):
        """Valid single sequence yields monomers, reaction, polymer, no errors."""
        errors = parser.parse(VALID_SINGLE)
        assert len(errors) == 0
        assert parser.monomers == [
            "CNC(C)(C(=O)O)C(C)(C)C",
            "CNC(C(=O)O)C(C)C",
        ]
        assert parser.reaction == (
            "([N&X3;H2,H1;!$(NC=*):1].[C&X3:2](=O)[O&X2&H1])."
            "([N&X3;H2,H1;!$(NC=*):3].[C&X3:4](=O)[O&X2&H1])>>"
            "(*-[N&X3&!$(NC=*):1].[C&X3:2](=O)[N&X3&!$(NC=*):3]."
            "[C&X3:4](=O)-*)"
        )
        assert parser.polymer == (
            "*C(=O)C(C(C)C)N(C)C(=O)C(C)(N(*)C)C(C)(C)C"
        )

    def test_parse_invalid_single_sequence(self, parser):
        """Invalid single sequence returns parse error and null attributes."""
        errors = parser.parse(INVALID_SINGLE)
        assert parser.monomers is None
        assert parser.reaction is None
        assert parser.polymer is None
        assert errors == ["<E-parse-|invalid-number-of-parts|>"]


# ---------------------------------------------------------------------------
# NestedPolyScriptParser (multi-sequence, |||-separated)
# ---------------------------------------------------------------------------


class TestNestedPolyScriptParser:
    """Tests for nested (|||-separated) PolyScriptParser."""

    @pytest.fixture
    def parser(self):
        return NestedPolyScriptParser()

    def test_parse_valid_nested_sequence(self, parser):
        """Valid nested sequence yields two sub-sequences, no errors."""
        errors = parser.parse(VALID_NESTED)
        assert errors == [[], []]

        assert len(parser.sequences) == 2

        # First sub-sequence
        m1, r1, p1 = parser.sequences[0]
        assert m1 == [
            "CNC(C)(C(=O)O)C(C)(C)C",
            "CNC(CC(=O)O)C(=O)O",
        ]
        assert r1 == (
            "([N&X3;H2,H1;!$(NC=*):1].[C&X3:2](=O)[O&X2&H1])."
            "([N&X3;H2,H1;!$(NC=*):3].[C&X3:4](=O)[O&X2&H1])>>"
            "(*-[N&X3&!$(NC=*):1].[C&X3:2](=O)[N&X3&!$(NC=*):3]."
            "[C&X3:4](=O)-*)"
        )
        assert p1 == (
            "*C(=O)CC(C(=O)O)N(C)C(=O)C(C)(N(*)C)C(C)(C)C"
        )

        # Second sub-sequence
        m2, r2, p2 = parser.sequences[1]
        assert m2 == [
            "*C(=O)CC(C(=O)O)N(C)C(=O)C(C)(N(*)C)C(C)(C)C",
            "none",
        ]
        assert r2 == "[C&X3:1](=O)[O&X2&H1,F,Cl,Br,I]>>[C&X3:1](=O)-*"
        assert p2 == (
            "*C(=O)CC(C(*)=O)N(C)C(=O)C(C)(N(*)C)C(C)(C)C"
        )

    def test_parse_invalid_nested_sequence(self, parser):
        """Invalid nested sequence reports errors for the broken segment."""
        errors = parser.parse(INVALID_NESTED)
        assert errors == [["<E-parse-|too-many-parts|>"]]

        # The first (broken) sub-sequence has null attributes
        m0, r0, p0 = parser.sequences[0]
        assert m0 is None
        assert r0 is None
        assert p0 is None
