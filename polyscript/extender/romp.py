"""Ring-opening metathesis polymerisation (ROMP) by insertion.

ROMP is a chain-growth process: a metal carbene opens a cycloalkene and the
ring becomes one repeat unit of the growing chain.  The repeat-unit map key

``[C&X3&H1&R:1]=[C&X3&H1&R:2][C&X4&R:3]>>(*=[C&X3&H1:1].[C&X3&H1:2]=*)``

cuts the metathesis ``C1=C2`` double bond and leaves both cut carbons as
double-bonded ends.  The allylic ``C2-C3`` single bond is *not* cut.  That cut
is topology-only: it cannot be turned into a generic reaction SMARTS for every
ring size.  A product template that breaks the ring must describe the path
between ``C1`` and ``C2``; for a variable-size ring that path is made of
unmapped atoms, and RDKit partitions it across dot-separated product
components when a second reactant is present, returning the opened monomer as
two fragments.

This driver therefore performs the cut programmatically:

* **initiation** opens one monomer at ``C1=C2``, caps ``C1`` with the
  initiator through a double bond, and leaves ``C2`` double-bonded to a dummy
  atom as the growing carbene end;
* **propagation** opens the next monomer at ``C1=C2``, forms a double bond
  between the chain's ``C2`` end and the new ``C1``, and leaves the new ``C2``
  as the growing end;
* **termination** replaces the dummy with the terminator through a double
  bond.

No heavy atoms are synthesised beyond the configured initiator/terminator
fragments.  Each repeat unit keeps one double bond, but the double bond is at
the junction between repeat units (the real ROMP topology), not inside the
ring path.
:class:`~polyscript.extender.chain_growth.ChainGrowthReaction`
(``chain_init`` / ``chain_extend`` / ``chain_term`` / ``extend``) so the same
benchmark and tests can drive it.
"""

from rdkit import Chem

from polyscript.extender.base import BaseReaction, prepare_monomers


class RingOpeningMetathesisReaction(BaseReaction):
    """Build ROMP chains by opening and inserting cycloalkene rings."""

    def __init__(
        self,
        ring_bond_smarts,
        initiator="C",
        terminator="C",
        rxn_type="romp",
        max_products: int | None = 1,
    ):
        super().__init__(rxn_type)
        if max_products is not None and (
            not isinstance(max_products, int) or max_products < 1
        ):
            raise ValueError("max_products must be a positive integer or None")

        query_smarts = str(ring_bond_smarts).split(">")[0]
        self.query = Chem.MolFromSmarts(query_smarts)
        if self.query is None:
            raise ValueError(
                f"Invalid ROMP ring-bond SMARTS: {query_smarts!r}"
            )
        if self.query.GetNumAtoms() < 3:
            raise ValueError(
                "ROMP ring-bond SMARTS must map the double-bond atom pair "
                "and the adjacent allylic ring atom"
            )
        self.initiator = self._fragment(initiator)
        self.terminator = self._fragment(terminator)
        self.max_products = max_products
        self._return_intermediate = False
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
        """Build ROMP chains with ``length`` propagation steps.

        ``length == 0`` returns the terminated initiation product.  Each
        propagation step inserts one cycloalkene ring.  The returned products
        contain no dummy atoms.
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
            if trace:
                intermediates.append([Chem.MolToSmiles(mol) for mol in chains])

        terminated = self.chain_term(monomers, chains)
        final = [Chem.MolToSmiles(mol) for mol in terminated]

        if trace:
            return {"intermediates": intermediates, "final": final}
        return final

    def chain_init(self, monomers):
        """Open each matching monomer and cap one side of the ring."""
        pool = prepare_monomers(monomers)
        products = []
        for mol in pool:
            for match in self._matches(mol):
                try:
                    products.append(self._init_monomer(Chem.Mol(mol), match))
                except Exception:
                    continue
                if (
                    self.max_products is not None
                    and len(products) >= self.max_products
                ):
                    self._chains = self._trim(products)
                    return self._chains
        self._chains = self._trim(products)
        return self._chains

    def chain_extend(self, monomers, chains=None):
        """Insert one additional opened monomer into every chain."""
        if chains is None:
            chains = self._chains
        if not chains:
            return []

        pool = prepare_monomers(monomers)
        extended = []
        for chain in chains:
            if not self._has_wildcard(chain):
                raise ValueError(
                    "The chain passed to ROMP propagation has no growing end"
                )
            for mol in pool:
                for match in self._matches(mol):
                    try:
                        extended.append(
                            self._insert(
                                Chem.Mol(chain),
                                Chem.Mol(mol),
                                match,
                            )
                        )
                    except Exception:
                        continue
                    if (
                        self.max_products is not None
                        and len(extended) >= self.max_products
                    ):
                        self._chains = self._trim(extended)
                        return self._chains
        self._chains = self._trim(extended)
        return self._chains

    def chain_term(self, monomers, chains=None):
        """Replace the growing dummy atom with the terminator fragment."""
        if chains is None:
            chains = self._chains
        terminated = []
        for chain in chains:
            try:
                terminated.append(self._terminate(Chem.Mol(chain)))
            except Exception:
                continue
            if (
                self.max_products is not None
                and len(terminated) >= self.max_products
            ):
                break
        self._chains = self._trim(terminated)
        return self._chains

    # ------------------------------------------------------------------
    # internals
    # ------------------------------------------------------------------
    @staticmethod
    def _fragment(smiles: str):
        mol = Chem.MolFromSmiles(str(smiles))
        if mol is None or mol.GetNumAtoms() == 0:
            raise ValueError(f"Invalid ROMP fragment: {smiles!r}")
        return mol

    def _matches(self, mol) -> list:
        return list(mol.GetSubstructMatches(self.query, uniquify=True))

    @staticmethod
    def _has_wildcard(mol) -> bool:
        return any(
            atom.GetAtomicNum() == 0 for atom in mol.GetAtoms()
        )

    @staticmethod
    def _dummy_index(mol) -> int | None:
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 0:
                return atom.GetIdx()
        return None

    @staticmethod
    def _attach_fragment(
        rw, atom_index: int, fragment, bond_type=Chem.BondType.SINGLE
    ) -> None:
        """Append *fragment* and bond its first atom to *atom_index*."""
        offset = rw.GetNumAtoms()
        for atom in fragment.GetAtoms():
            new_atom = Chem.Atom(atom.GetAtomicNum())
            new_atom.SetFormalCharge(atom.GetFormalCharge())
            new_atom.SetIsAromatic(atom.GetIsAromatic())
            rw.AddAtom(new_atom)
        for bond in fragment.GetBonds():
            rw.AddBond(
                offset + bond.GetBeginAtomIdx(),
                offset + bond.GetEndAtomIdx(),
                bond.GetBondType(),
            )
        rw.AddBond(atom_index, offset, bond_type)

    def _init_monomer(self, mol, match):
        c1, c2 = match[:2]
        rw = Chem.RWMol(mol)
        # Metathesis cleaves the C=C double bond, not the allylic single bond.
        rw.RemoveBond(c1, c2)
        self._attach_fragment(
            rw, c1, self.initiator, Chem.BondType.DOUBLE
        )
        dummy = rw.AddAtom(Chem.Atom(0))
        rw.AddBond(c2, dummy, Chem.BondType.DOUBLE)
        product = rw.GetMol()
        Chem.SanitizeMol(product)
        return product

    def _insert(self, chain, monomer, match):
        dummy = self._dummy_index(chain)
        if dummy is None:
            raise ValueError("growing chain has no dummy atom")
        neighbour = chain.GetAtomWithIdx(dummy).GetNeighbors()[0].GetIdx()

        rw = Chem.RWMol(chain)
        rw.RemoveAtom(dummy)
        if neighbour > dummy:
            neighbour -= 1
        clean = rw.GetMol()

        offset = clean.GetNumAtoms()
        combined = Chem.RWMol(Chem.CombineMols(clean, monomer))
        c1, c2 = (index + offset for index in match[:2])
        combined.RemoveBond(c1, c2)
        # Join the previous repeat unit and the new one through a double bond.
        combined.AddBond(neighbour, c1, Chem.BondType.DOUBLE)
        new_dummy = combined.AddAtom(Chem.Atom(0))
        combined.AddBond(c2, new_dummy, Chem.BondType.DOUBLE)

        product = combined.GetMol()
        Chem.SanitizeMol(product)
        return product

    def _terminate(self, chain):
        dummy = self._dummy_index(chain)
        if dummy is None:
            raise ValueError("chain has no dummy atom to terminate")
        neighbour = chain.GetAtomWithIdx(dummy).GetNeighbors()[0].GetIdx()
        rw = Chem.RWMol(chain)
        rw.RemoveAtom(dummy)
        if neighbour > dummy:
            neighbour -= 1
        self._attach_fragment(
            rw, neighbour, self.terminator, Chem.BondType.DOUBLE
        )
        product = rw.GetMol()
        Chem.SanitizeMol(product)
        return product

    @staticmethod
    def _dedupe(mols: list) -> list:
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
        unique = self._dedupe(mols)
        if self.max_products is None:
            return unique
        return unique[: self.max_products]
