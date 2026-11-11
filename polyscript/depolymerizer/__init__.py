"""
PolyScript Depolymerizer: parse PSMILES sequences, run reverse
(depolymerization) and forward (repolymerization) reactions to verify
monomer recovery.

Usage:
    from depolymerizer.poly_script import PolyScriptDepolymerizer
    dp = PolyScriptDepolymerizer()
    results = dp.validate("*Oc1ccccc1")
    for r in results:
        print(r)
"""

from dataclasses import dataclass, field
from typing import Optional

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

from polyscript.depolymerizer.reaction_reference import polym_to_forward_rxn, polym_to_reverse_rxn
from polyscript.utils.executor import BaseExecutor
from polyscript.utils.logger import INFO
from polyscript.utils.parsers import PolyScriptParser

# this should suppress rdkit logging (removed once the rdkit log flooding is fixed)
RDLogger.logger().setLevel(RDLogger.CRITICAL)

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------


def _strip_outer_parens(s: str) -> str:
    """Strip a single matching pair of outer parentheses from a string.

    Only removes parentheses when they enclose the entire string
    (i.e. the outermost pair balances and ends at the final character).
    Nested parentheses are handled correctly.

    Args:
        s: Input string that may have outer parentheses.

    Returns:
        str: The string with outer parens removed, or the original
        string if no matching outer pair exists.
    """
    if s.startswith("(") and s.endswith(")"):
        depth = 0
        for i, ch in enumerate(s):
            if ch == "(":
                depth += 1
            elif ch == ")":
                depth -= 1
            if depth == 0 and i == len(s) - 1:
                return s[1:-1]
            if depth == 0 and i < len(s) - 1:
                break
    return s


def _match_polymer_substructure(cand_smiles: str, poly_smarts: str) -> bool:
    """Check if *cand_smiles* contains the polymer SMARTS substructure."""
    cand_mol = Chem.MolFromSmiles(cand_smiles)
    if cand_mol is None:
        return False
    try:
        Chem.SanitizeMol(cand_mol)
    except Exception:
        return False
    poly_mol = Chem.MolFromSmarts(poly_smarts)
    if poly_mol is None:
        return False
    return cand_mol.HasSubstructMatch(poly_mol)


def _split_reactant_groups(rxn_smarts: str) -> list[list[str]]:
    """Split the reactant side of a reaction SMARTS into template groups.

    Handles:
      - Simple:    'A.B'        -> [['A'], ['B']]
      - Grouped:   '(A.B).(C.D)' -> [['A', 'B'], ['C', 'D']]
      - Single:    'A'          -> [['A']]

    Each group is a list of SMARTS templates (OR logic within a group).
    """
    reactant = rxn_smarts.split(">>")[0]
    groups: list[list[str]] = []
    i = 0
    while i < len(reactant):
        if reactant[i] == "(":
            depth = 1
            j = i + 1
            while j < len(reactant) and depth > 0:
                if reactant[j] == "(":
                    depth += 1
                elif reactant[j] == ")":
                    depth -= 1
                j += 1
            inner = reactant[i + 1 : j - 1]
            groups.append(inner.split("."))
            i = j
            if i < len(reactant) and reactant[i] == ".":
                i += 1
        else:
            j = i
            while j < len(reactant) and reactant[j] != ".":
                j += 1
            groups.append([reactant[i:j]])
            i = j
            if i < len(reactant) and reactant[i] == ".":
                i += 1
    return groups


def _clean_reactant(mol) -> Optional[Chem.Mol]:
    """Strip query properties from a mol via SMILES round-trip, then sanitize."""
    try:
        smi = Chem.MolToSmiles(mol)
        clean = Chem.MolFromSmiles(smi)
        if clean is not None:
            Chem.SanitizeMol(clean)
        return clean
    except Exception:
        return None


def _init_rings(mol) -> bool:
    """Initialize RingInfo on a mol (lightweight — no valence changes)."""
    try:
        Chem.FastFindRings(mol)
        return True
    except Exception:
        return False


def _mol_from_smiles(smi: str) -> Optional[Chem.Mol]:
    """Parse and sanitize a SMILES string. Returns None on failure."""
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    try:
        Chem.SanitizeMol(mol)
        return mol
    except Exception:
        return None


# ---------------------------------------------------------------------------
# Result dataclass
# ---------------------------------------------------------------------------


@dataclass
class DepolymerizeResult:
    """Result of a successful reverse -> forward round-trip."""

    polymer: str
    reaction: str
    monomers: list[str]
    input_monomers: list[str] = field(default_factory=list)
    polymer_type: str = ""
    polymer_pattern: str = ""
    seq_index: int = -1  # index into input sequence list

    @property
    def is_exact_match(self) -> bool:
        """True if recovered monomers match input (ignoring 'none' entries)."""
        out = set(m for m in self.monomers if m != "none")
        inp = set(m for m in self.input_monomers if m != "none")
        return out == inp

    def __repr__(self) -> str:
        """Return a compact string representation for debugging."""
        return (
            f"DepolymerizeResult(idx={self.seq_index}, "
            f"polymer={self.polymer!r}, "
            f"reaction={self.reaction!r}, monomers={self.monomers!r})"
        )


# ---------------------------------------------------------------------------
# Main depolymerizer class
# ---------------------------------------------------------------------------


def _normalize_psmiles(smi: str) -> Optional[str]:
    """Normalize a polymer SMILES to canonical form for comparison."""
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    try:
        return Chem.MolToSmiles(mol, canonical=True)
    except Exception:
        return None


class PolyScriptDepolymerizer(BaseExecutor):
    """Validate monomer recovery through reverse -> forward reaction round-trip.

    Uses ``polym_to_forward_rxn`` and ``polym_to_reverse_rxn`` from
    ``reaction_reference`` for multi-step depolymerization.
    """

    _log_level = INFO  # class-level default, override via set_log_level()

    def __init__(self, **kwargs):
        """Initialize the depolymerizer with reaction reference tables.

        Loads polymer backbone patterns and their associated forward
        reaction SMARTS from ``reaction_reference``, stripping outer
        parentheses from pattern keys for matching.

        Args:
            **kwargs: Passed to ``BaseExecutor.__init__``.
        """
        super().__init__(**kwargs)
        # Map: stripped_pattern -> list of forward reactions
        self._pattern_to_rxns: dict[str, list[str]] = {}
        self._poly_patterns: list[str] = []
        for key, rxns in polym_to_forward_rxn.items():
            stripped = _strip_outer_parens(key)
            self._poly_patterns.append(stripped)
            self._pattern_to_rxns[stripped] = rxns

        self.logger.debug(
            "PolyScriptDepolymerizer initialised (%d polymer patterns)",
            len(self._poly_patterns),
        )

    # -- reverse reaction -------------------------------------------------

    def _try_reverse_reaction(
        self, rxn: str, polymer_smiles: str
    ) -> Optional[list[str]]:
        """Run multi-step reverse reaction to recover monomer candidates.

        Returns a list of monomer SMILES (without wildcard '*') or None.
        """
        polym_pattern = rxn.split(">>")[-1]
        mol = _mol_from_smiles(polymer_smiles)
        if mol is None:
            self.logger.debug("mol is None for SMILES: %s", polymer_smiles)
            return None

        if (
            polym_pattern not in polym_to_reverse_rxn
            and polym_pattern not in polym_to_forward_rxn
        ):
            self.logger.debug(
                "polym_pattern not in reverse or forward rxn: %s", polym_pattern
            )
            return None

        rxn_fwd_list = polym_to_forward_rxn.get(polym_pattern, [])
        rxn_fwd_idx = -1
        for i, smarts in enumerate(rxn_fwd_list):
            if smarts == rxn:
                rxn_fwd_idx = i
                break
        if rxn_fwd_idx == -1:
            self.logger.debug("rxn_fwd_idx == -1 for rxn: %s", rxn)
            return None

        rxn_smarts_list = polym_to_reverse_rxn[polym_pattern][rxn_fwd_idx]

        reactants = [mol]
        master_products: list = []
        for smarts_list in rxn_smarts_list:
            all_products = []
            for smarts in smarts_list:
                if smarts is None:
                    self.logger.debug("smarts is None in reverse rxn")
                    return None
                reaction = AllChem.ReactionFromSmarts(smarts)
                for r in reactants:
                    try:
                        r_clean = _clean_reactant(r)
                        if r_clean is None:
                            continue
                        products = reaction.RunReactants([r_clean])
                        flat_products = [
                            p for d in products for p in d if p is not None
                        ]
                        sanitized = []
                        for p in flat_products:
                            if _init_rings(p):
                                sanitized.append(p)
                        if sanitized:
                            all_products.extend(sanitized)
                    except Exception:
                        continue
            if all_products:
                reactants = all_products
                master_products.extend(all_products)

        master_smiles = [Chem.MolToSmiles(m) for m in master_products if m is not None]
        # Remove wildcard-containing entries
        return list(set(r for r in master_smiles if "*" not in r))

    # -- forward reaction -------------------------------------------------

    def _create_forward_combinations(
        self, monomer_list: list[str], rxn: str
    ) -> list[tuple[list[str], tuple]]:
        """Try forward reaction with classified monomer combinations.

        Returns a list of (product_smiles_list, monomer_tuple) pairs.
        """
        reaction = AllChem.ReactionFromSmarts(rxn)
        rxn_products: list[tuple[list[str], tuple]] = []

        groups = _split_reactant_groups(rxn)

        if len(groups) == 1:
            # Single-monomer reaction (homopolymerization)
            for monomer in monomer_list:
                mol = _mol_from_smiles(monomer)
                if mol is None:
                    continue
                products = reaction.RunReactants([mol])
                flat_products: list[str] = []
                for sublist in products:
                    for p in sublist:
                        if p is None:
                            continue
                        if _init_rings(p):
                            try:
                                flat_products.append(
                                    Chem.MolToSmiles(p, canonical=True)
                                )
                            except Exception:
                                continue
                if flat_products:
                    rxn_products.append((flat_products, (monomer,)))

        elif len(groups) >= 2:
            # Multi-monomer: classify by template match, cross-react groups
            tmpl_mols_0 = [Chem.MolFromSmarts(t) for t in groups[0] if t]
            tmpl_mols_1 = [Chem.MolFromSmarts(t) for t in groups[1] if t]

            set_0: list[tuple[str, Chem.Mol]] = []
            set_1: list[tuple[str, Chem.Mol]] = []

            for monomer in monomer_list:
                mol = _mol_from_smiles(monomer)
                if mol is None:
                    continue
                in_0 = any(
                    mol.HasSubstructMatch(t) for t in tmpl_mols_0 if t is not None
                )
                in_1 = any(
                    mol.HasSubstructMatch(t) for t in tmpl_mols_1 if t is not None
                )
                if in_0:
                    set_0.append((monomer, mol))
                if in_1:
                    set_1.append((monomer, mol))

            for smi0, mol0 in set_0:
                for smi1, mol1 in set_1:
                    if smi0 == smi1:
                        continue
                    products = reaction.RunReactants([Chem.Mol(mol0), Chem.Mol(mol1)])
                    flat_products = []
                    for sublist in products:
                        for p in sublist:
                            if p is None:
                                continue
                            if _init_rings(p):
                                try:
                                    flat_products.append(
                                        Chem.MolToSmiles(p, canonical=True)
                                    )
                                except Exception:
                                    continue
                    if flat_products:
                        rxn_products.append((flat_products, (smi0, smi1)))

        return rxn_products

    # -- validation -------------------------------------------------------

    def _find_match(
        self,
        rxn_products: list[tuple[list[str], tuple]],
        polymer_smiles: str,
    ) -> Optional[tuple[list[str], tuple]]:
        """Find the (products, monomers) pair that contains polymer_smiles.
        Normalizes SMILES to handle [*] vs * notation differences."""
        target = _normalize_psmiles(polymer_smiles)
        if target is None:
            return None
        for products, monomer in rxn_products:
            for p in products:
                if _normalize_psmiles(p) == target:
                    return (products, monomer)
        return None

    # -- main pipeline ----------------------------------------------------

    def convert_psmiles(self, psmiles: str) -> list[DepolymerizeResult]:
        """Convert a raw PSMILES by finding matching reactions from the
        pattern library and reverse-depolymerizing.

        Returns all successful reaction pathways.
        """
        mol = _mol_from_smiles(psmiles)
        if mol is None:
            return []
        results: list[DepolymerizeResult] = []
        for poly_pattern in self._poly_patterns:
            if not _match_polymer_substructure(psmiles, poly_pattern):
                self.logger.debug("pattern mismatch: %s", poly_pattern)
                continue
            self.logger.debug("pattern matched: %s", poly_pattern)
            for rxn in self._pattern_to_rxns.get(poly_pattern, []):
                rev = self._try_reverse_reaction(rxn, psmiles)
                self.logger.debug("reverse reaction result: %s", rev)
                if rev is None:
                    continue
                rxn_products = self._create_forward_combinations(rev, rxn)
                match = self._find_match(rxn_products, psmiles)
                if match is not None:
                    _p, monomers = match
                    results.append(
                        DepolymerizeResult(
                            polymer=psmiles,
                            reaction=rxn,
                            monomers=list(monomers),
                            polymer_pattern=poly_pattern,
                        )
                    )
        return results

    def convert_all(self, psmiles_list: list[str]):
        """Batch-convert a list of raw PSMILES strings.

        Runs ``convert_psmiles`` on each input and collects all
        successful pathways.

        Args:
            psmiles_list: List of PSMILES strings to convert.

        Returns:
            tuple: ``(all_results, stats)`` where *all_results* is a
            list of ``DepolymerizeResult`` objects (with ``seq_index``
            set to the position in the input list) and *stats* is a
            dict with keys ``total``, ``success``, ``failed``.
        """
        n_total = len(psmiles_list)
        self.logger.info("Starting batch convert (%d PSMILES)", n_total)

        all_results = []
        seen: set[int] = set()
        for idx, ps in enumerate(psmiles_list):
            for r in self.convert_psmiles(ps):
                r.seq_index = idx
                all_results.append(r)
                seen.add(idx)

            if (idx + 1) % max(n_total // 10, 1) == 0 or idx == n_total - 1:
                self.logger.info(
                    "convert progress: %d/%d (found %d)",
                    idx + 1, n_total, len(seen),
                )
                self.logger.flush()

        stats = {
            "total": n_total,
            "success": len(seen),
            "failed": n_total - len(seen),
        }
        self.logger.info(
            "Batch convert complete: %d/%d matched", stats["success"], n_total
        )
        self.logger.flush()
        return all_results, stats

    def validate(self, polym_sequence: str) -> list[DepolymerizeResult]:
        """Parse a PolyScript sequence and validate monomer recovery.

        Returns a list of ``DepolymerizeResult`` for each successful
        reverse -> forward round-trip.
        """
        parser = PolyScriptParser()
        errors = parser.parse(polym_sequence)
        if errors:
            self.logger.warning("Parse errors: %s", errors)
            return []

        if not parser.polymer or not parser.reaction:
            self.logger.debug("Missing polymer or reaction after parse")
            return []

        results: list[DepolymerizeResult] = []

        for poly_pattern in self._poly_patterns:
            if not _match_polymer_substructure(parser.polymer, poly_pattern):
                continue

            rev_monomers = self._try_reverse_reaction(parser.reaction, parser.polymer)
            if rev_monomers is None:
                continue

            rxn_products = self._create_forward_combinations(
                rev_monomers, parser.reaction
            )
            match = self._find_match(rxn_products, parser.polymer)
            if match is not None:
                _products, matched_monomers = match
                results.append(
                    DepolymerizeResult(
                        polymer=parser.polymer,
                        reaction=parser.reaction,
                        monomers=list(matched_monomers),
                        input_monomers=list(parser.monomers),
                        polymer_pattern=poly_pattern,
                    )
                )

        return results

    def validate_all(
        self, polym_sequences: list[str]
    ) -> tuple[list[DepolymerizeResult], dict]:
        """Run ``validate`` on a batch. Returns ALL pathways per sequence
        with ``seq_index`` set so results can be mapped back to input.

        Returns (all_results, stats).  Stats keys: total, success, failed.
        """
        n_total = len(polym_sequences)
        self.logger.info("Starting batch validation (%d sequences)", n_total)

        all_results: list[DepolymerizeResult] = []
        seen_success: set[int] = set()

        for idx, seq in enumerate(polym_sequences):
            seq_results = self.validate(seq)
            if seq_results:
                seen_success.add(idx)
                for r in seq_results:
                    r.seq_index = idx
                    all_results.append(r)

            if (idx + 1) % max(n_total // 10, 1) == 0 or idx == n_total - 1:
                self.logger.info(
                    "validate progress: %d/%d (found %d)",
                    idx + 1, n_total, len(seen_success),
                )
                self.logger.flush()

        stats = {
            "total": n_total,
            "success": len(seen_success),
            "failed": n_total - len(seen_success),
        }
        self.logger.info(
            "Batch validation complete: %d/%d succeeded",
            stats["success"], n_total,
        )
        self.logger.flush()
        return all_results, stats

    # -- display -----------------------------------------------------------

    def print_results(self, results: list[DepolymerizeResult]) -> None:
        """Log a list of ``DepolymerizeResult`` at INFO level."""
        if not results:
            self.logger.info("No results.")
            return
        for i, r in enumerate(results):
            self.logger.info(
                "Case %d: monomers=%s reaction=%s polymer=%s",
                i + 1, r.monomers, r.reaction, r.polymer,
            )
        self.logger.flush()
