"""
Full chemical validator for PolyScript sequences.

Validates monomers (SMILES), reaction (SMARTS), reaction execution, and
polymer product matching via RDKit.

Handles multi-sequences separated by ``|||`` in the ``polym`` column.
A datapoint is considered **entirely valid** only when **every** sub-sequence
passes all validation stages.
"""

from __future__ import annotations

from types import SimpleNamespace

from rdkit import Chem
from rdkit.Chem import AllChem
from polyscript.utils.parsers import NestedPolyScriptParser

class SeqValidator:
    """Validation for generated sequences in the PolyScript format.

    Format::

        ${monomer1}+{monomer2}=>{reaction_smarts}=>{polymer_smiles}$

    Multi-sequences are separated by ``|||``::

        ${m1}+{m2}=>{r1}=>{p1}$|||${m1}+{m2}=>{r2}=>{p2}$
    """

    def __init__(self):
        """Initialise a SeqValidator with a NestedPolyScriptParser."""
        self.parser = NestedPolyScriptParser()
        # Holds the current sub-sequence being validated during
        # multi-sequence iteration.  Set by validate().
        self._seq: SimpleNamespace | None = None
        self.errors: list[str] = []

    # ------------------------------------------------------------------
    # Per-sub-sequence validation helpers
    # ------------------------------------------------------------------

    def _validate_monomers(self) -> None:
        """Validate that each monomer is a parseable SMILES string.

        Appends ``<E-chem-|invalid-monomer-smiles|>`` to ``self.errors``
        for any monomer that RDKit cannot parse.  Monomers whose value is
        ``"none"`` are skipped.
        """
        for monomer in self._seq.monomers:
            if monomer == "none":
                continue
            mol = Chem.MolFromSmiles(monomer)
            if mol is None:
                self.errors.append("<E-chem-|invalid-monomer-smiles|>")

    def _validate_reaction(self):
        """Parse and validate the reaction SMARTS.

        Returns:
            An RDKit ``ChemicalReaction`` object, or ``None`` if the
            SMARTS is invalid (an error token is appended to
            ``self.errors``).
        """
        try:
            reaction = AllChem.ReactionFromSmarts(self._seq.reaction)
            if reaction is None:
                self.errors.append("<E-chem-|invalid-reaction-smarts|>")
            return reaction
        except Exception as e:
            self.errors.append(f"<E-chem-|reaction-loading-error|>: {str(e)}")
            return None

    def _validate_reaction_execution(self, reaction):
        """Run the reaction on the validated monomers and sanitize products.

        Args:
            reaction: An RDKit ``ChemicalReaction`` object.

        Returns:
            A list of product tuples (each tuple contains sanitized
            ``Mol`` objects), or ``None`` if execution fails (an error
            token is appended to ``self.errors``).
        """
        try:
            monomers = [
                Chem.MolFromSmiles(s) for s in self._seq.monomers if s != "none"
            ]
            products = reaction.RunReactants(monomers)
            if len(products) == 0:
                self.errors.append("<E-chem-|no-products|>")
                return None
            # Sanitize each product mol
            sanitized = []
            for tup in products:
                sanitized_tup = []
                for mol in tup:
                    try:
                        Chem.SanitizeMol(mol)
                    except Exception:
                        pass
                    sanitized_tup.append(mol)
                sanitized.append(tuple(sanitized_tup))
            return sanitized
        except Exception as e:
            self.errors.append(f"<E-chem-|reaction-execution-error|>: {str(e)}")
            return None

    def _has_polymer_match(self, products) -> bool:
        """Check whether any reaction product matches the expected polymer.

        Args:
            products: List of product tuples from reaction execution.

        Returns:
            ``True`` if at least one product's canonical SMILES matches
            ``self._seq.polymer``; ``False`` otherwise (and an error token
            is appended to ``self.errors``).
        """
        # Validate polymer SMILES
        try:
            mol = Chem.MolFromSmiles(self._seq.polymer)
            if mol is None:
                self.errors.append(
                    "<E-chem-|polymer-validation-error|>: invalid polymer SMILES"
                )
                return False
        except Exception as e:
            self.errors.append(f"<E-chem-|polymer-validation-error|>: {str(e)}")
            return False

        # Flatten products (list of tuples) into a flat list of mols
        flat_products = [p for t in products for p in t]
        for p in flat_products:
            if not isinstance(p, str):
                p = Chem.MolToSmiles(p, canonical=True)
            if self._seq.polymer == p:
                return True
        return False

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def _validate(
        self,
        monomers: list[str],
        reaction_smarts: str,
        polymer_smiles: str,
    ) -> list[str]:
        """Run the full validation pipeline on already-parsed values.

        Useful when the caller has pre-parsed (monomers, reaction, polymer)
        and wants to skip the parse step.

        Returns a flat list of error tokens for **this** sub-sequence;
        an empty list means this sub-sequence is valid.
        """
        self._seq = SimpleNamespace(
            monomers=monomers,
            reaction=reaction_smarts,
            polymer=polymer_smiles,
        )
        self._validate_monomers()
        reaction = self._validate_reaction()
        if reaction is not None:
            products = self._validate_reaction_execution(reaction)
            if products is not None:
                if not self._has_polymer_match(products):
                    self.errors.append("<E-chem-|no-polymer-match|>")
        return self.errors

    def validate(self, sequence: str) -> list[str]:
        """Validate **every** sub-sequence (monomers, reaction, execution,
        polymer match).  Returns a flat list of error tokens; empty means
        all sub-sequences are valid.
        """
        self.errors = []

        errs_per_seq = self.parser.parse(sequence)
        # Flatten per-sub-sequence parse errors
        flat_parse_errs = [e for sub in errs_per_seq for e in sub]
        if flat_parse_errs:
            self.errors = flat_parse_errs
            return self.errors

        for monomers, reaction_smarts, polymer_smiles in self.parser.sequences:
            self._validate(monomers, reaction_smarts, polymer_smiles)
        return self.errors

    def is_valid(self, sequence: str) -> bool:
        """Return True when **every** sub-sequence passes all validation."""
        return len(self.validate(sequence)) == 0


if __name__ == "__main__":
    validator = SeqValidator()
    cand = "{CNC_invalid_smiles(C)(C(=O)O)C(C)(C)C}+{CNC(CC(=O)O)C(=O)O}=>([N&X3;H2,H1;!$(NC=*):1].[C&X3:2](=O)[O&X2&H1]).([N&X3;H2,H1;!$(NC=*):3].[C&X3:4](=O)[O&X2&H1])>>(*-[N&X3&!$(NC=*):1].[C&X3:2](=O)[N&X3&!$(NC=*):3].[C&X3:4](=O)-*)=>*C(=O)CC(C(=O)O)N(C)C(=O)C(C)(N(*)C)C(C)(C)C|||{*C(=O)CC(C(=O)O)N(C)C(=O)C(C)(N(*)C)C(C)(C)C}+{none}=>[C&X3:1](=O)[O&X2&H1,F,Cl,Br,I]>>[C&X3:1](=O)-*=>*C(=O)CC(C(*)=O)N(C)C(=O)C(C)(N(*)C)C(C)(C)C"
    errors = validator.validate(cand)
    print("errors >>", errors)
