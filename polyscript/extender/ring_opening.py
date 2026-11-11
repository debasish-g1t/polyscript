"""Ring-opening polymerisation (ROP/ROCOP) by open-then-insert propagation.

Ring-opening polymerisation is an insertion/chain-growth process: the growing
chain's active end opens an intact ring and inserts the whole ring as one
repeat unit.  It is not a generic step-growth condensation, even though the
two ends of the ring-opened monomer can look like an ``AB`` step-growth pair.

This driver therefore separates the two chemical events:

1. **Ring opening (cycle 1).**  The map key opens the cyclic monomer into its
   real open form: a lactone becomes a hydroxyacid, a lactam an aminoacid, and
   a cyclic anhydride + epoxide pair becomes the alternating anhydride/epoxide
   hydroxyacid unit.  Ring opening is a single-reactant reaction for ROP and a
   two-reactant reaction for ROCOP.  RDKit handles ring opening reliably in
   this form because the ring topology is described by one reaction template.
2. **Propagation (cycle 2+).**  A ROP-specific coupling template reacts the
   growing chain's terminal group with the opened monomer's complementary
   terminal group.  The monomer's other end is carried over unchanged, so the
   chain keeps real terminal groups throughout growth.

Why not attack the intact ring directly in the propagation SMARTS?  A generic
direct-attack template would have to break the ring in its *product* template.
For a variable-size ring the path between the scissile atoms is made of
unmapped atoms; when the template is used with two reactants, RDKit partitions
that path across the dot-separated product components and returns the opened
monomer as two fragments.  The open-then-insert form avoids that limitation and
keeps the propagation template specific to the opened end groups instead of
falling back to the generic step-growth exclusions (``!$(NC=*)`` and
``!$([O,S]C=*)``) that reject amidine/thioimidate amines and acyloxy O-H ends.

``length`` has the same meaning as in :class:`StepGrowthReaction`: the number
of opened residues (for ROCOP, one anhydride + one epoxide pair counts as one
residue).  ``length == 1`` returns the opened monomer; every further cycle
inserts one more opened monomer.

The map key is supplied as ``ring_open_smarts`` by
:func:`polyscript.extender.map.reaction_config`.  Propagation templates live in
:data:`polyscript.extender.map.RING_OPENING_PROPAGATIONS`.
"""

from itertools import product

from rdkit import Chem

from polyscript.extender.base import prepare_monomers
from polyscript.extender.step_growth import StepGrowthReaction


class RingOpeningReaction(StepGrowthReaction):
    """Grow a chain by ring-opening and inserting opened monomers."""

    def __init__(
        self,
        prop_smarts,
        ring_open_smarts,
        reagents=None,
        opening_reagents=None,
        rxn_type="ring_opening",
        max_products: int | None = 1,
    ):
        super().__init__(
            coupling_smarts=prop_smarts,
            ring_open_smarts=ring_open_smarts,
            reagents=reagents,
            rxn_type=rxn_type,
            max_products=max_products,
        )
        # Reagents consumed by the ring-opening step (for example water in a
        # hydrolytic ROCOP initiation).  They are not added to the propagation
        # pool, because they are consumed by the opening reaction.
        self.opening_reagents = prepare_monomers(opening_reagents or [])
        # Propagation reactions are exposed under their chemistry name; the
        # inherited ``coupling_rxns`` remain the same compiled templates so
        # the provenance helpers and the benchmark can use one attribute.
        self.prop_rxns = self.coupling_rxns

    # ------------------------------------------------------------------
    # pool preparation
    # ------------------------------------------------------------------
    def _prepare_pool(self, monomers) -> list:
        """Open the monomer pool and add reagents.

        A one-reactant ring-opening template is applied to each monomer.  A
        multi-reactant template (ROCOP) is applied to every ordered
        combination of monomers; the first successful opened unit per row is a
        real alternating residue.  If no combination opens, the raw pool is
        retained so the failure is reported as a propagation failure rather
        than silently dropping the row.
        """
        pool = prepare_monomers(monomers)
        ring_open = self.ring_open_rxn
        if ring_open is not None:
            if ring_open.GetNumReactantTemplates() == 1:
                pool = self._dedupe(
                    [self._open_monomer(mol) for mol in pool]
                )
            else:
                opened = self._open_combinations(
                    self._dedupe(pool + self.opening_reagents)
                )
                if opened:
                    pool = opened
        return self._merge_pool(pool, self.reagents)

    def _open_combinations(self, pool: list) -> list:
        """Run a multi-reactant ring-opening template over monomer pairs."""
        ring_open = self.ring_open_rxn
        if ring_open is None:
            return []
        templates = ring_open.GetNumReactantTemplates()
        opened = []
        for combo in product(pool, repeat=templates):
            try:
                results = ring_open.RunReactants(list(combo))
            except Exception:
                continue
            for result in results:
                if len(result) != 1:
                    continue
                candidate = result[0]
                if any(
                    atom.GetAtomicNum() == 0
                    for atom in candidate.GetAtoms()
                ):
                    candidate = self._replace_wildcards(candidate)
                if candidate is None:
                    continue
                try:
                    Chem.SanitizeMol(candidate)
                except Exception:
                    continue
                opened.append(candidate)
        return self._dedupe(opened)
