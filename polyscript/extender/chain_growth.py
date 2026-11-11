"""Chain-growth polymerisation from init/propagation/termination SMARTS.

The three reaction SMARTS describe one growing active chain:

* ``init_smarts`` consumes one or more monomer units and produces a chain that
  carries a single wildcard (``*``) at its active end.
* ``prop_smarts`` reacts that wildcard with the next monomer unit (or two
  monomer units, if the template declares two free reactants) and places a new
  wildcard at the end of the extended chain.
* ``term_smarts`` consumes the final wildcard and caps the chain with the
  terminating group.

Termination is chemistry-specific and is passed in per reaction entry; it is
never a global cap.  This driver implements chain-growth only.  Polyolefin
(addition) polymerisation is the primary configured use, but the same
initiation/propagation/termination lifecycle applies to other chain-growth
families such as ROMP.  A condensation/step-growth chain has reactive groups
at both ends and must use a driver whose termination preserves those end
groups.

``max_products`` caps the frontier after initiation, after every propagation
step, and after termination.  The default is ``1``: the driver follows one
deterministic chain instead of enumerating every possible sequence.  Pass
``None`` to restore exhaustive enumeration (useful for small exploratory
cases, but it grows exponentially with ``length``).

The growing-chain reactant is located automatically: in propagation and
termination the unique reactant template containing a wildcard atom is filled
with the current chain, while every other reactant slot is filled from the
monomer pool.  This keeps the class independent of any particular monomer
chemistry, so no code change is needed when the SMARTS in
:mod:`polyscript.extender.map` are updated.
"""

from itertools import product

from rdkit import Chem

from polyscript.extender.base import (
    BaseReaction,
    get_rxn,
    normalize_reactant_wildcards,
    prepare_monomers,
)


class ChainGrowthReaction(BaseReaction):
    """Build chain-growth polymers by iteratively applying reaction SMARTS."""

    def __init__(
        self,
        init_smarts,
        prop_smarts,
        term_smarts,
        rxn_type="chain_growth",
        max_products: int | None = 1,
    ):
        super().__init__(rxn_type)
        if max_products is not None and (
            not isinstance(max_products, int) or max_products < 1
        ):
            raise ValueError("max_products must be a positive integer or None")
        self.max_products = max_products
        self._return_intermediate = False
        self.init_rxn = self._compile(init_smarts)
        self.prop_rxn = self._compile(prop_smarts)
        self.term_rxn = self._compile(term_smarts)
        # Most recently built chains, used when chain_extend/chain_term are
        # called without an explicit ``chains`` argument.
        self._chains = []

    # ------------------------------------------------------------------
    # public API
    # ------------------------------------------------------------------
    @property
    def return_intermediate(self) -> bool:
        """Whether :meth:`extend` returns the growth trace by default.

        Set this after loading a reaction from ``rxn_map`` to toggle tracing
        without re-creating the reaction or recompiling its SMARTS::

            reaction = ChainGrowthReaction(**config["params"])
            reaction.return_intermediate = True
            result = reaction.extend(monomers, length)
        """
        return self._return_intermediate

    @return_intermediate.setter
    def return_intermediate(self, value: bool) -> None:
        if not isinstance(value, bool):
            raise TypeError("return_intermediate must be a bool")
        self._return_intermediate = value

    def extend(
        self,
        monomers: list,
        length: int,
        return_intermediate: bool | None = None,
    ):
        """Build terminated chains containing ``length`` propagation steps.

        Args:
            monomers: Monomer pool.  Each item may be a SMILES string or an
                ``rdkit.Chem.Mol``.  Every propagation step may pick any
                member of the pool, including the same one repeatedly, so a
                multi-monomer pool enumerates all sequences of that length.
            length: Number of propagation steps performed after initiation.
                ``0`` returns the terminated initiation product.
            return_intermediate: Per-call override.  ``None`` (default)
                uses the instance setting
                :attr:`~ChainGrowthReaction.return_intermediate`.  ``True`` also
                returns the growing product after initiation and after every
                propagation step; ``False`` returns only the final molecules.
                Intermediates retain the wildcard at the growing end; the
                ``final`` entry is the terminated molecule.  The trace is not
                built when tracing is disabled.

        Returns:
            list[str]: Canonical SMILES of the terminated chains when tracing
            is disabled.  At most ``max_products`` chains are retained after
            initiation and after every propagation step.

            dict: When tracing is enabled::

                {
                    "intermediates": [
                        ["...*"],   # after initiation (n=0)
                        ["...*"],   # after propagation step 1
                        ...
                    ],
                    "final": ["..."],  # after termination
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

        chains = self.chain_init(monomers)
        if trace:
            intermediates.append([Chem.MolToSmiles(mol) for mol in chains])

        for _ in range(length):
            if chains:
                chains = self.chain_extend(monomers, chains)
            else:
                # Keep one trace entry per requested step even if growth
                # stopped early.
                if trace:
                    intermediates.append([])
                continue
            if trace:
                intermediates.append([Chem.MolToSmiles(mol) for mol in chains])

        terminated = self.chain_term(monomers, chains)
        final = [Chem.MolToSmiles(mol) for mol in terminated]

        if trace:
            return {"intermediates": intermediates, "final": final}
        return final

    def chain_init(self, monomers: list):
        """Apply the initiation reaction to the monomer pool.

        The init reaction must not contain a wildcard on its reactant side;
        the wildcard is introduced in its product to mark the growing end.

        Returns:
            list[rdkit.Chem.Mol]: Initiated chains ending in a wildcard,
            capped to ``max_products``.
        """
        pool = self._prepare_monomers(monomers)
        self._chains = self._run(
            self.init_rxn, pool, limit=self.max_products
        )
        return self._chains

    def chain_extend(self, monomers: list, chains: list | None = None):
        """Extend every chain by one propagation step.

        Args:
            monomers: Monomer pool available for the propagation reaction.
            chains: Chains to extend.  Defaults to the result of the most
                recent :meth:`chain_init`/:meth:`chain_extend` call.

        Returns:
            list[rdkit.Chem.Mol]: Extended chains ending in a wildcard,
            capped to ``max_products``.  Once the cap is reached, remaining
            parent chains and reaction assignments are not evaluated.
        """
        if chains is None:
            chains = self._chains
        if not chains:
            return []

        pool = self._prepare_monomers(monomers)
        extended = []
        for chain in chains:
            remaining = (
                None
                if self.max_products is None
                else self.max_products - len(extended)
            )
            if remaining is not None and remaining <= 0:
                break
            extended.extend(
                self._run(self.prop_rxn, pool, chain=chain, limit=remaining)
            )

        self._chains = self._trim(extended)
        return self._chains

    def chain_term(self, monomers: list, chains: list | None = None):
        """Apply the termination reaction to close the wildcard ends.

        Args:
            monomers: Monomer pool available to terminating reactants (empty
                or ignored when the termination SMARTS needs only the chain).
            chains: Chains to terminate.  Defaults to the result of the most
                recent :meth:`chain_init`/:meth:`chain_extend` call.

        Returns:
            list[rdkit.Chem.Mol]: Terminated chains, capped to
            ``max_products``.
        """
        if chains is None:
            chains = self._chains
        if not chains:
            return []

        pool = self._prepare_monomers(monomers)
        terminated = []
        for chain in chains:
            remaining = (
                None
                if self.max_products is None
                else self.max_products - len(terminated)
            )
            if remaining is not None and remaining <= 0:
                break
            terminated.extend(
                self._run(self.term_rxn, pool, chain=chain, limit=remaining)
            )

        self._chains = self._trim(terminated)
        return self._chains

    # ------------------------------------------------------------------
    # internals
    # ------------------------------------------------------------------
    @classmethod
    def _compile(cls, smarts):
        """Compile reaction SMARTS after normalising growing-end wildcards."""
        return get_rxn(cls._normalize_reactant_wildcards(smarts))[0]

    @staticmethod
    def _normalize_reactant_wildcards(smarts: str) -> str:
        """Rewrite growing-end wildcards on the reactant side as ``[#0]``."""
        return normalize_reactant_wildcards(smarts)

    @staticmethod
    def _prepare_monomers(monomers) -> list:
        """Normalise monomers to Mol objects and remove duplicates."""
        return prepare_monomers(monomers)

    @staticmethod
    def _has_wildcard(mol) -> bool:
        """True if *mol* contains a ``*`` (dummy) atom."""
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 0:
                return True
        return False

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

    def _run(self, rxn, monomers: list, chain=None, limit=None) -> list:
        """Run one reaction over all monomer combinations.

        When *chain* is given, the unique wildcard-bearing reactant template
        is filled with the chain and all remaining templates are filled from
        *monomers*.  Without a chain (initiation), every template is filled
        from *monomers*.

        *limit* stops collection as soon as that many unique products have
        been produced, so a small frontier cap also avoids materialising the
        combinatorial explosion.
        """
        num_templates = rxn.GetNumReactantTemplates()
        if num_templates == 0:
            raise ValueError("Reaction SMARTS must declare at least one reactant")

        wildcard_slots = [
            index
            for index in range(num_templates)
            if self._has_wildcard(rxn.GetReactantTemplate(index))
        ]

        if chain is None:
            if wildcard_slots:
                raise ValueError(
                    "Initiation SMARTS must not contain a wildcard reactant; "
                    "the growing-end wildcard belongs in the product."
                )
            chain_slot = None
        else:
            if len(wildcard_slots) != 1:
                raise ValueError(
                    "Propagation/termination SMARTS must contain exactly one "
                    "wildcard-bearing reactant template (the growing chain); "
                    f"found {len(wildcard_slots)}."
                )
            if not self._has_wildcard(chain):
                raise ValueError(
                    "The chain passed to a propagation/termination step has "
                    "no wildcard attachment point."
                )
            chain_slot = wildcard_slots[0]

        free_slots = [i for i in range(num_templates) if i != chain_slot]
        if free_slots and not monomers:
            raise ValueError(
                "The monomer pool is empty but the reaction needs "
                f"{len(free_slots)} free reactant(s)."
            )

        products = []
        seen = set()
        for combo in product(monomers, repeat=len(free_slots)):
            reactants = [None] * num_templates
            if chain_slot is not None:
                reactants[chain_slot] = chain
            for slot, mol in zip(free_slots, combo):
                reactants[slot] = mol

            try:
                results = rxn.RunReactants(reactants)
            except Exception:
                continue

            for result in results:
                if not result:
                    continue
                mol = result[0]
                if mol is None or mol.GetNumAtoms() == 0:
                    continue
                try:
                    Chem.SanitizeMol(mol)
                except Exception:
                    continue
                smiles = Chem.MolToSmiles(mol)
                if smiles in seen:
                    continue
                seen.add(smiles)
                products.append(mol)
                if limit is not None and len(products) >= limit:
                    return products

        return products
