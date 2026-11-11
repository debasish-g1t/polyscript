"""Step-growth polymerisation by stepwise coupling of real monomer ends.

The map key describes the polymerisation chemistry.  A step-growth chain is
built by repeatedly running a *coupling* reaction that consumes **one**
functional group from the growing chain and **one** functional group from a
monomer.  The monomer's other functional groups are not part of the reaction
template, so RDKit carries them over unchanged.  The chain therefore keeps the
monomers' real terminal groups throughout growth: an acid chloride monomer
leaves ``-COCl`` ends, a difluoro monomer leaves ``Ar-F`` ends, a diisocyanate
leaves ``N=C=O`` ends, and so on.

No wildcard or end-cap synthesis is involved.  The wildcard convention in
:mod:`polyscript.extender.map` is only the repeat-unit/cut-point notation used
by the depolymerizer; it is not a place to invent terminal chemistry.

``length`` is the number of monomer residues in the final chain.  The driver
starts from the supplied monomers (one residue each) and performs
``length - 1`` coupling steps.  For an ``AA`` + ``BB`` pair, ``length = 2 * n``
therefore gives ``n`` units of each monomer if the sequence alternates.
"""

from itertools import product

from rdkit import Chem

from polyscript.extender.base import BaseReaction, get_rxn, prepare_monomers


class StepGrowthReaction(BaseReaction):
    """Build step-growth polymers by repeatedly coupling real chain ends."""

    def __init__(
        self,
        coupling_smarts,
        ring_open_smarts=None,
        reagents=None,
        rxn_type="step_growth",
        max_products: int | None = 1,
    ):
        super().__init__(rxn_type)
        if max_products is not None and (
            not isinstance(max_products, int) or max_products < 1
        ):
            raise ValueError("max_products must be a positive integer or None")

        if isinstance(coupling_smarts, str):
            coupling_smarts = [coupling_smarts]
        if not coupling_smarts:
            raise ValueError("At least one coupling reaction SMARTS is required")

        self.max_products = max_products
        self._return_intermediate = False
        self.coupling_rxns = get_rxn(coupling_smarts)
        self.ring_open_rxn = (
            get_rxn(ring_open_smarts)[0] if ring_open_smarts else None
        )
        self.reagents = prepare_monomers(reagents or [])
        self._chains = []

    # ------------------------------------------------------------------
    # public API
    # ------------------------------------------------------------------
    @property
    def return_intermediate(self) -> bool:
        """Whether :meth:`extend` returns the growth trace by default."""
        return self._return_intermediate

    @return_intermediate.setter
    def return_intermediate(self, value: bool) -> None:
        if not isinstance(value, bool):
            raise TypeError("return_intermediate must be a bool")
        self._return_intermediate = value

    def extend(
        self,
        monomers,
        length: int,
        return_intermediate: bool | None = None,
    ):
        """Build chains containing ``length`` monomer residues.

        Args:
            monomers: Monomer pool.  Each item may be a SMILES string or an
                ``rdkit.Chem.Mol``.  Both components of an ``AA`` + ``BB``
                pair are supplied in the same pool; the coupling reactions
                determine which combinations can react.
            length: Number of monomer residues in the final chain.  ``1``
                returns the supplied monomers, ``0`` returns an empty chain.
            return_intermediate: Per-call override.  When enabled, returns the
                chain after initialization and after every coupling step.

        Returns:
            list[str]: Canonical SMILES of the final chains.

            dict: When tracing is enabled::

                {
                    "intermediates": [
                        ["..."],   # one monomer residue (n=1)
                        ["..."],   # after one coupling (n=2)
                        ...
                    ],
                    "final": ["..."],
                }
        """
        if not isinstance(length, int) or length < 0:
            raise ValueError("length must be a non-negative integer")
        if return_intermediate is not None and not isinstance(
            return_intermediate, bool
        ):
            raise TypeError("return_intermediate must be a bool or None")

        trace = (
            self.return_intermediate
            if return_intermediate is None
            else return_intermediate
        )
        intermediates: list = []

        if length == 0:
            self._chains = []
            if trace:
                return {"intermediates": intermediates, "final": []}
            return []

        pool = self._prepare_pool(monomers)
        if not pool:
            raise ValueError("No monomers provided for step-growth extension")

        chains = self._dedupe(pool)
        if trace:
            intermediates.append([Chem.MolToSmiles(mol) for mol in chains])

        for _ in range(length - 1):
            if chains:
                chains = self._couple(chains, pool)
            if trace:
                intermediates.append([Chem.MolToSmiles(mol) for mol in chains])
            if not chains:
                continue

        self._chains = self._trim(chains)
        final = [Chem.MolToSmiles(mol) for mol in self._chains]

        if trace:
            return {"intermediates": intermediates, "final": final}
        return final

    def couple(self, chain, monomers):
        """Run every coupling reaction on *chain* and a monomer/reagent pool.

        Both reactant orders and all chain/template assignments are tried.
        The returned molecules are single connected products with no
        wildcards; the monomers' unreacted functional groups are carried over.
        """
        pool = prepare_monomers(monomers)
        products = []
        for rxn in self.coupling_rxns:
            products.extend(self._run_coupling(rxn, chain, pool))
        return self._dedupe(products)

    # ------------------------------------------------------------------
    # internals
    # ------------------------------------------------------------------
    def _prepare_pool(self, monomers) -> list:
        """Return the monomer pool used for growth.

        Step growth uses the supplied monomers unchanged (plus reagents), or
        their ring-opened forms when ``ring_open_smarts`` is configured.  The
        ring-opening hook is overridden by
        :class:`~polyscript.extender.ring_opening.RingOpeningReaction`.
        """
        pool = prepare_monomers(monomers)
        if self.ring_open_rxn is not None:
            pool = self._dedupe(
                [self._open_monomer(mol) for mol in pool]
            )
        return self._merge_pool(pool, self.reagents)

    def prepare_pool(self, monomers) -> list:
        """Public wrapper around :meth:`_prepare_pool`.

        The benchmark uses this to obtain exactly the same monomer pool as
        :meth:`extend`, including ring opening and reagents.
        """
        return self._prepare_pool(monomers)

    @staticmethod
    def _merge_pool(monomers: list, reagents: list) -> list:
        """Merge monomers and reagents, removing duplicate structures."""
        merged = []
        seen = set()
        for mol in list(monomers) + list(reagents):
            smiles = Chem.MolToSmiles(mol)
            if smiles in seen:
                continue
            seen.add(smiles)
            merged.append(Chem.Mol(mol))
        return merged

    def _open_monomer(self, mol):
        """Convert a ring monomer into its real open-form monomer.

        A lactone becomes a hydroxyacid and a lactam becomes an aminoacid.
        The wildcard at an ``O``/``N``/``S`` end becomes an ``H`` (the free
        heteroatom end); the wildcard on a carbonyl carbon becomes an ``O``
        (the carboxylic acid end).  This is the ring-opened monomer, not a
        synthesized polymer end cap.
        """
        ring_open = self.ring_open_rxn
        if ring_open is None:
            return mol
        try:
            results = ring_open.RunReactants([mol])
        except Exception:
            return mol

        for result in results:
            if len(result) != 1:
                continue
            candidate = result[0]
            if any(
                atom.GetAtomicNum() == 0
                for atom in candidate.GetAtoms()
            ):
                opened = self._replace_wildcards(candidate)
            else:
                try:
                    Chem.SanitizeMol(candidate)
                    opened = candidate
                except Exception:
                    opened = None
            if opened is not None:
                return opened
        return mol

    @staticmethod
    def _replace_wildcards(mol):
        """Replace terminal wildcards with the open-form end atoms."""
        from rdkit.Chem import RWMol

        if mol is None:
            return None
        try:
            Chem.SanitizeMol(mol)
        except Exception:
            return None

        rw = RWMol(mol)
        dummies = [
            atom.GetIdx()
            for atom in rw.GetAtoms()
            if atom.GetAtomicNum() == 0
        ]
        if len(dummies) != 2:
            return None

        replacements = []
        for index in dummies:
            atom = rw.GetAtomWithIdx(index)
            if atom.GetDegree() != 1:
                return None
            neighbour = atom.GetNeighbors()[0]
            carbonyl = any(
                bond.GetBondType() == Chem.BondType.DOUBLE
                and bond.GetOtherAtom(neighbour).GetAtomicNum() == 8
                for bond in neighbour.GetBonds()
            )
            replacements.append(("O" if carbonyl else "H", neighbour.GetIdx()))

        for offset, (_kind, neighbour) in enumerate(replacements):
            rw.GetAtomWithIdx(neighbour).SetAtomMapNum(100 + offset)
        for index in sorted(dummies, reverse=True):
            rw.RemoveAtom(index)

        opened = rw.GetMol()
        rw = RWMol(opened)
        for offset, (kind, _neighbour) in enumerate(replacements):
            target = next(
                atom.GetIdx()
                for atom in rw.GetAtoms()
                if atom.GetAtomMapNum() == 100 + offset
            )
            rw.GetAtomWithIdx(target).SetAtomMapNum(0)
            if kind == "O":
                oxygen = Chem.Atom(8)
                oxygen_index = rw.AddAtom(oxygen)
                rw.AddBond(target, oxygen_index, Chem.BondType.SINGLE)

        opened = rw.GetMol()
        try:
            Chem.SanitizeMol(opened)
        except Exception:
            return None
        if any(atom.GetAtomicNum() == 0 for atom in opened.GetAtoms()):
            return None
        return opened

    def _couple(self, chains: list, pool: list) -> list:
        """Add one monomer residue to each chain, respecting the frontier cap."""
        extended = []
        seen = set()
        errors = []

        for chain in chains:
            for rxn in self.coupling_rxns:
                for product_mol in self._run_coupling(rxn, chain, pool):
                    smiles = Chem.MolToSmiles(product_mol)
                    if smiles in seen:
                        continue
                    seen.add(smiles)
                    extended.append(product_mol)
                    if (
                        self.max_products is not None
                        and len(extended) >= self.max_products
                    ):
                        return extended
            if not extended:
                errors.append(
                    "no coupling product for chain "
                    f"{Chem.MolToSmiles(chain)!r}"
                )

        if not extended and errors:
            raise ValueError(
                "All step-growth couplings failed; first error: "
                f"{errors[0]}"
            )
        return extended

    def _run_coupling(self, rxn, chain, pool: list) -> list:
        """Run one coupling reaction with *chain* in each possible slot."""
        num_templates = rxn.GetNumReactantTemplates()
        if num_templates == 0:
            raise ValueError("Coupling SMARTS must declare at least one reactant")

        products = []
        seen = set()
        for chain_slot in range(num_templates):
            free_slots = [
                index
                for index in range(num_templates)
                if index != chain_slot
            ]
            for combo in product(pool, repeat=len(free_slots)):
                reactants = [None] * num_templates
                reactants[chain_slot] = chain
                for slot, mol in zip(free_slots, combo):
                    reactants[slot] = mol
                try:
                    results = rxn.RunReactants(reactants)
                except Exception:
                    continue

                for result in results:
                    if len(result) != 1:
                        continue
                    candidate = result[0]
                    if candidate is None or candidate.GetNumAtoms() == 0:
                        continue
                    try:
                        Chem.SanitizeMol(candidate)
                    except Exception:
                        continue
                    if len(Chem.GetMolFrags(candidate)) != 1:
                        continue
                    smiles = Chem.MolToSmiles(candidate)
                    if smiles in seen:
                        continue
                    seen.add(smiles)
                    products.append(candidate)
                    if (
                        self.max_products is not None
                        and len(products) >= self.max_products
                    ):
                        return products
        return products

    @staticmethod
    def _dedupe(mols: list) -> list:
        """Remove duplicate molecules by canonical SMILES, preserving order."""
        unique = []
        seen = set()
        for mol in mols:
            smiles = Chem.MolToSmiles(mol)
            if smiles in seen:
                continue
            seen.add(smiles)
            unique.append(mol)
        return unique

    def _trim(self, mols: list) -> list:
        """Deduplicate and apply the ``max_products`` frontier cap."""
        unique = self._dedupe(mols)
        if self.max_products is None:
            return unique
        return unique[: self.max_products]
