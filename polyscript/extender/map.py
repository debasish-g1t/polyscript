"""Reaction SMARTS registry and mechanism recipes for chain extension.

``rxn_map`` is keyed by the reaction SMARTS recorded for a polymer family.
Every value is pure data; the reaction objects are compiled lazily by the
driver selected for the entry's ``type``.  Each entry also carries a
``mechanism`` code (``A1``–``R5``) that ties it to
``benchmarks/mechanism_table.tex`` and the persistent example manifest
``benchmarks/mechanism_examples.json``.  The mechanisms are:

``addition``
    Chain-growth from ``init_smarts`` / ``prop_smarts`` / ``term_smarts``.
    The growing end is a wildcard; initiation and termination supply the
    end groups.  The configured polyolefin entries use hydrogen caps: the
    termination product removes the wildcard, so no extra methyl carbon is
    invented.  This is the chain-transfer/β-hydride-elimination end group of
    coordination polymerisation and the H-ended product of radical
    disproportionation, rather than a methyl-capped chain.

``step_growth``
    Condensation by stepwise coupling of real monomer functional groups.
    ``coupling_smarts`` consumes one functional group from the chain and one
    from a monomer; the monomer's other groups are carried over unchanged.
    This is not wildcard propagation and does not synthesize end groups.  It
    covers ester, amide, urethane, imide, oxazolidinone, carbonate, NAS aryl
    ether, and oxidative phenol coupling.

``ring_open`` / ``ring_opening``
    Ring-opening polymerisation (ROP) and copolymerisation (ROCOP) in the
    open-then-condense form.  The map key opens the cyclic monomer into its
    real open form; the opened monomer is then inserted with a ROP-specific
    coupling.  The final chain carries real terminal groups.

    This is a **stoichiometric/topological** model, not an active-end insertion
    simulation.  It does not track a growing dummy or a separate termination
    event; it reproduces the terminated chain directly.  The stoichiometry
    includes the initiator explicitly:

    * a water-initiated lactone/lactam opens to a hydroxyacid/aminoacid, so
      the first residue has one more oxygen than the cyclic monomer;
    * the ROCOP opening is the balanced anhydride + epoxide + water reaction,
      so the first residue has one more oxygen than the two monomers.

    Later residues add only the monomer's own atoms.  The extra first-residue
    oxygen is the water initiator, not an invented atom.  ``ring_open`` is the
    legacy spelling of the same driver with a generic coupling; the registry
    uses ``ring_opening`` with the ROP-specific couplings.

``romp``
    Ring-opening metathesis polymerisation.  The reactant side of the map key
    names the ring alkene and an allylic atom so a variable-size ring can be
    found.  The driver cuts the ``C=C`` double bond (the metathesis event),
    caps one cut end with the initiator, leaves the other as a double-bonded
    growing end, and joins each new repeat unit through a new double bond.
    Termination caps the growing end with the configured terminator.  It does
    not keep the ring alkene and cleave the allylic single bond.

Why ring-opening propagation has its own couplings
--------------------------------------------------
A cyclic monomer is not an ``AB`` step-growth monomer.  Opening a lactone gives
a hydroxyacid, opening a lactam an aminoacid, and opening an anhydride + epoxide
pair gives a hydroxyacid with an ester linkage.  Those open forms can look like
step-growth monomers, but the propagation event is an insertion into an intact
ring.  The generic ester/amide guards used to reject two real opened ends:

* a broad amide guard rejected an amine attached to an amidine/thioimidate
  carbon.  A thioimidate lactam (``R11``) opens to exactly that amine.
* a broad ester guard rejected an acyloxy O-H whose carbonyl is attached to
  nitrogen.  A cyclic carbamate/carbonic anhydride such as the
  oxazine-2,6-dione in ``R12`` opens to that carbamic-acid end.

The ROP-specific couplings in :data:`RING_OPENING_PROPAGATIONS` keep the
chemically correct guard and add only the missing carbamic-acid case:

* ``lactone_rop`` excludes carboxylic acid, ester, carbonate, sulfonic,
  boronic, phosphoric, oxime, and peroxide O-H from the nucleophile slot, but
  allows an alcohol and the O-H of a carbamic acid (``N-C(=O)-O-H``).
* ``lactam_rop`` excludes amide/imide/thioamide, sulfonamide, and
  phosphoramide nitrogen, but allows an amine attached to an
  amidine/thioimidate carbon (``C=N``).

A literal direct attack on the intact ring is the chemically cleanest form,
but a single reaction SMARTS cannot describe the variable-length ring path
between the scissile atoms.  When such a template is used with two reactants,
RDKit partitions the unmapped path across the dot-separated product components
and returns the opened monomer as two fragments.  The open-then-condense form
avoids that RDKit limitation and still produces the same terminated chain as
an active-end insertion model for the known polymer examples.  Small rings
whose atoms can all be mapped are the exception: the three-membered epoxide in
``epoxide_rop`` is attacked directly and remains a ``step_growth`` entry.

ROCOP oxygen accounting
-----------------------
The anhydride + epoxide dataset key is a repeat-unit topology, not a balanced
elementary step.  The registry uses :data:`ROCOP_OPENING_SMARTS`, a balanced
three-reactant opening (anhydride + epoxide + water):

* the anhydride ring O becomes the ester O;
* the water O becomes the remaining acid O-H;
* the epoxide O becomes the terminal alcohol.

``opening_reagents`` supplies the water to the opening step only; it is consumed
and not added to the propagation pool.  For PES:

::

    succinic anhydride (3 O) + ethylene oxide (1 O) + water (1 O)
        -> HOOC-CH2-CH2-COO-CH2-CH2-OH (5 O)

The extra oxygen is the initiator.  One cycle consumes one anhydride + one
epoxide pair, which matches the alternating ROCOP repeat unit.

Chemoselectivity guards
-----------------------
The step-growth couplings use explicit functional-group guards instead of the
old negative lookaheads:

* an alcohol/phenol/thiol nucleophile must not be an acyloxy O-H (carboxylic
  acid, ester, carbonate), a sulfonic/boronic/phosphoric acid, an oxime or
  hydroxylamine O-H, a peroxide, or an enol ether;
* a free amine must not be an amide/imide/thioamide (``N-C(=O/S)``),
  sulfonamide, or phosphoramide, but amidine/thioimidate ``N-C(=N)`` is
  allowed;
* NAS aryl ether coupling requires the aryl halide to be para-activated by a
  ketone, sulfone, or nitrile; unactivated aryl halides do not match.  The
  activating group may be on the halide ring or on a directly attached
  phenyl spacer (PEEK/PES/PAEN analogs);
* imide formation requires a cyclic (five- or six-membered) anhydride, so an
  acyclic anhydride cannot be scrambled into an acyclic imide;
* ROCOP opening accepts five- and six-membered anhydrides whose ring atoms
  flanking the carbonyls may be carbon, nitrogen, oxygen, or sulfur, so
  diglycolic and other heteroatom anhydrides open as well;
* oxidative phenol coupling requires both ortho positions of the coupling
  ring to be blocked, so unsubstituted phenol is not accepted as a linear
  PPO monomer;
* polyester coupling supports both acid/acyl-halide esterification and
  methyl-ester transesterification (the DMT route).  The transesterification
  template is methyl-specific so it cannot cleave backbone esters;
* isocyanate + epoxide and epoxide ROP attack the less hindered epoxide
  carbon; regioisomeric mixtures are only produced where the two carbons are
  equivalent by substitution;
* carbonate linkages use phosgene as the carbonyl source
  (``2 R-OH + Cl2C=O -> R-O-C(=O)-O-R + 2 HCl``), not direct CO insertion,
  which would remove H2 without an oxidant.

Mechanistic scope
-----------------
The step-growth recipes reproduce the **net** repeat-unit chemistry.  They do
not track transient intermediates whose lifetime is shorter than one coupling
cycle:

* a carbonate coupling is written as the balanced two-alcohol + phosgene net
  reaction; a real process can also stop at a chloroformate or half-carbonate
  ester;
* imide formation is written as the net cyclodehydration; the open
  poly(amic acid) is not represented as a separate residue;
* the R15/R16 polyolefin keys name two alkene comonomers.  The template
  inserts one unit of each in the same propagation step, so the model gives a
  strictly alternating ABAB copolymer that matches the dataset's comonomer
  pair.  It is not a conjugated-diene 1,4-addition (which would retain one
  backbone C=C);
* epoxide ROP is represented by the water-initiated coupling of the intact
  epoxide.  The product stoichiometry is correct, but the driver does not
  simulate a metal alkoxide active end.

Oxidative phenol coupling
-------------------------
Poly(phenylene oxide) synthesis is a step-growth C-O coupling, not a
ring-opening reaction.  It is configured with ``phenol_oxidative`` in
:data:`STEP_GROWTH_COUPLINGS`; the phenolic O-H of one ring couples to the
para C-H of another ring, and the second ring's O-H is carried over as the
real terminal group.  Both ortho positions of the coupling ring must be
blocked, so unsubstituted phenol is not accepted as a linear PPO monomer.

Terminal-group and atom provenance
----------------------------------
* ``addition`` entries: end groups come from ``init_smarts``/``term_smarts``.
  Wildcards in the repeat-unit key are main-chain cut points, not end groups.
  The configured polyolefin termination removes the wildcard, giving a
  hydrogen end rather than an extra methyl carbon.
* ``step_growth`` entries: a coupling consumes one functional group from the
  chain and one from a monomer.  RDKit carries the rest of the monomer through
  unchanged, so the polymer keeps the monomer's exact unreacted end groups
  (``-COOH``, ``-COCl``, ``-OH``, ``-SH``, ``-NH2``, ``-N=C=O``, ``-N=C=S``,
  ``Ar-F``/``Ar-Cl``, epoxide, anhydride, phenol, ...).  No end group is
  synthesized from a wildcard.  The carbonate recipe consumes an explicit
  phosgene reagent; its carbonyl carbon and oxygen are reagent atoms, not
  invented atoms.
* ``ring_open``/``ring_opening`` entries: the opened monomer is the source of
  end groups.  Ring opening is the only inference; the balanced ROCOP opening
  additionally consumes an explicit water initiator.  The first residue's
  extra oxygen is that initiator.  No heavy atom is synthesized from a
  wildcard.
* ``romp`` entries: the initiator and terminator fragments are explicit recipe
  parameters.  Ring opening cuts the ``C=C`` bond and does not synthesize a
  heavy atom.
* Entries mapped to ``None`` would be unconfigured.  The registry currently has
  no ``None`` entries; a mechanism that needs new chemistry should be added
  here with an explicit recipe rather than guessed at run time.

Use :func:`reaction_config` to build the driver kwargs.  For ``ring_open`` and
``ring_opening`` entries it injects the map key as ``ring_open_smarts``; for
``romp`` entries it injects the map key as ``ring_bond_smarts``.
"""

# ---------------------------------------------------------------------------
# Chemoselectivity guards
# ---------------------------------------------------------------------------
# The generic couplings below are deliberately narrow.  In particular they do
# not treat every O-H or N-H as a nucleophile:
#
# * an acyloxy O-H (carboxylic acid, ester, carbonate) is not an alcohol;
# * a sulfonic, boronic, phosphoric, hydroxylamine/oxime, or peroxide O-H is
#   not an alcohol;
# * an amide/imide/thioamide N-H is not a free amine, while an amidine or
#   thioimidate N-H is (its carbon is C=N, not C=O/S);
# * a sulfonamide/phosphoramide N-H is not a free amine.
#
# The same guards are used for the ROP propagation templates so an opened
# monomer is inserted by its real end group rather than by a broad "any O-H"
# or "any N-H" rule.
_ACYL = r'[C&X3](=[O,S])'


def _nucleophile(map_num: int, *, water: bool = False) -> str:
    """SMARTS for an alcohol/thiol/phenol nucleophile.

    ``water=True`` also accepts water (``O`` with two H and one connection),
    which is the explicit initiator for the epoxide step-growth entry.
    """
    connectivity = 'H1,H2' if water else 'X2;H1'
    return (
        f'[#8,#16;{connectivity};'
        f'!$([#8,#16]{_ACYL});'
        '!$([#8,#16][S,P,B,N]);'
        '!$([#8,#16][#8]);'
        '!$([#8,#16][C&X3]=[C,N])'
        f':{map_num}]'
    )


def _amine(map_num: int, *, primary_only: bool = False) -> str:
    """SMARTS for a free primary/secondary amine.

    Amidine and thioimidate N-H are allowed (C=N is not an acyl group).
    Amide, imide, thioamide, sulfonamide, and phosphoramide N-H are excluded.
    """
    hydrogens = 'H2' if primary_only else 'H2,H1'
    return (
        f'[N&X3;{hydrogens};'
        '!$([N][C&X3](=[O,S]));'
        '!$([N][S,P](=[O]))'
        f':{map_num}]'
    )


def _carbamic_oh(map_num: int) -> str:
    """SMARTS for the O-H of a carbamic acid (N-C(=O)-O-H)."""
    return (
        '[#8;X2;H1;'
        '$([#8][C&X3](=[O,S])[#7])'
        f':{map_num}]'
    )


# An aryl halide only undergoes uncatalysed SNAr when the ring is activated by
# a para carbonyl/sulfonyl/nitrile group.  The halogen is deliberately left
# unmapped: it is the leaving group and is removed by the reaction.  The ring
# atoms and the EWG are mapped so the aromatic ring is carried through intact.
#
# The dataset uses para-activated PEEK (diaryl ketone) and PES (diaryl
# sulfone) monomers.  Ketone, sulfone, and nitrile activation are covered, as
# are monomers where the EWG sits on a directly attached phenyl spacer.
# Fluoride templates come first so an activated Ar-F is preferred over an
# activated Ar-Cl/Br/I, matching SNAr reactivity (Ar-F >> Ar-Cl > Ar-Br >
# Ar-I).  Unactivated aryl halides do not match.
def _ar_nas_ketone(halogen: str) -> str:
    return (
        f'[c:1]1([{halogen}])[c:3][c:4][c:5]'
        '([C&X3:6](=[O:7])[#6:8])[c:9][c:10]1.'
        + _nucleophile(11) +
        '>>[c:1]1([#8,#16;X2:11])[c:3][c:4][c:5]'
        '([C&X3:6](=[O:7])[#6:8])[c:9][c:10]1'
    )


def _ar_nas_sulfone(halogen: str) -> str:
    return (
        f'[c:1]1([{halogen}])[c:3][c:4][c:5]'
        '([S&X4:6](=[O:7])(=[O:8])[#6:9])[c:12][c:13]1.'
        + _nucleophile(11) +
        '>>[c:1]1([#8,#16;X2:11])[c:3][c:4][c:5]'
        '([S&X4:6](=[O:7])(=[O:8])[#6:9])[c:12][c:13]1'
    )


def _ar_nas_nitrile(halogen: str) -> str:
    return (
        f'[c:1]1([{halogen}])[c:3][c:4][c:5]([C&X2:6]#[N&X1:7])'
        '[c:9][c:10]1.'
        + _nucleophile(11) +
        '>>[c:1]1([#8,#16;X2:11])[c:3][c:4][c:5]'
        '([C&X2:6]#[N&X1:7])[c:9][c:10]1'
    )


def _ar_nas_nitrile_ortho(halogen: str) -> str:
    """ArX ortho to a nitrile (2,6-difluorobenzonitrile PAEN topology)."""
    return (
        f'[c:1]1([{halogen}])[c:2]([C&X2:6]#[N&X1:7])[c:3][c:4][c:5]'
        '[c:10]1.'
        + _nucleophile(11) +
        '>>[c:1]1([#8,#16;X2:11])[c:2]([C&X2:6]#[N&X1:7])'
        '[c:3][c:4][c:5][c:10]1'
    )


def _ar_nas_ketone_biphenyl(halogen: str) -> str:
    """ArX para to a phenyl spacer that carries the ketone."""
    return (
        f'[c:1]1([{halogen}])[c:3][c:4][c:5](-[c:9]2[c:10][c:11]'
        '[c:12]([C&X3:6](=[O:7])[#6:8])[c:13][c:14]2)[c:15][c:16]1.'
        + _nucleophile(17) +
        '>>[c:1]1([#8,#16;X2:17])[c:3][c:4][c:5](-[c:9]2[c:10][c:11]'
        '[c:12]([C&X3:6](=[O:7])[#6:8])[c:13][c:14]2)[c:15][c:16]1'
    )


def _ar_nas_sulfone_biphenyl(halogen: str) -> str:
    """ArX para to a phenyl spacer that carries the sulfone."""
    return (
        f'[c:1]1([{halogen}])[c:3][c:4][c:5](-[c:9]2[c:10][c:11]'
        '[c:12]([S&X4:6](=[O:7])(=[O:8])[#6:18])[c:13][c:14]2)'
        '[c:15][c:16]1.'
        + _nucleophile(17) +
        '>>[c:1]1([#8,#16;X2:17])[c:3][c:4][c:5](-[c:9]2[c:10][c:11]'
        '[c:12]([S&X4:6](=[O:7])(=[O:8])[#6:18])[c:13][c:14]2)'
        '[c:15][c:16]1'
    )


def _ar_nas_nitrile_biphenyl(halogen: str) -> str:
    """ArX para to a phenyl spacer that carries the nitrile."""
    return (
        f'[c:1]1([{halogen}])[c:3][c:4][c:5](-[c:9]2[c:10][c:11]'
        '[c:12]([C&X2:6]#[N&X1:7])[c:13][c:14]2)[c:15][c:16]1.'
        + _nucleophile(17) +
        '>>[c:1]1([#8,#16;X2:17])[c:3][c:4][c:5](-[c:9]2[c:10][c:11]'
        '[c:12]([C&X2:6]#[N&X1:7])[c:13][c:14]2)[c:15][c:16]1'
    )


_KETONE_EWG = '[C&X3:6](=[O:7])[#6:8]'
_SULFONE_EWG = '[S&X4:6](=[O:7])(=[O:8])[#6:19]'
_NITRILE_EWG = '[C&X2:6]#[N&X1:7]'


def _ar_nas_naphthyl(halogen: str, ewg: str) -> str:
    """ArX on one naphthalene ring, EWG on the other ring."""
    return (
        f'[c:1]1([{halogen}])[c:10][c:11]2[c:12][c:13][c:14]'
        f'({ewg})[c:15][c:16]2[c:17][c:18]1.'
        + _nucleophile(9) +
        '>>[c:1]1([#8,#16;X2:9])[c:10][c:11]2[c:12][c:13][c:14]'
        f'({ewg})[c:15][c:16]2[c:17][c:18]1'
    )


_AR_NAS_KETONE_F = _ar_nas_ketone('F')
_AR_NAS_SULFONE_F = _ar_nas_sulfone('F')
_AR_NAS_NITRILE_F = _ar_nas_nitrile('F')
_AR_NAS_NITRILE_ORTHO_F = _ar_nas_nitrile_ortho('F')
_AR_NAS_KETONE_BIPHENYL_F = _ar_nas_ketone_biphenyl('F')
_AR_NAS_SULFONE_BIPHENYL_F = _ar_nas_sulfone_biphenyl('F')
_AR_NAS_NITRILE_BIPHENYL_F = _ar_nas_nitrile_biphenyl('F')
_AR_NAS_KETONE_X = _ar_nas_ketone('Cl,Br,I')
_AR_NAS_SULFONE_X = _ar_nas_sulfone('Cl,Br,I')
_AR_NAS_NITRILE_X = _ar_nas_nitrile('Cl,Br,I')
_AR_NAS_NITRILE_ORTHO_X = _ar_nas_nitrile_ortho('Cl,Br,I')
_AR_NAS_KETONE_BIPHENYL_X = _ar_nas_ketone_biphenyl('Cl,Br,I')
_AR_NAS_SULFONE_BIPHENYL_X = _ar_nas_sulfone_biphenyl('Cl,Br,I')
_AR_NAS_NITRILE_BIPHENYL_X = _ar_nas_nitrile_biphenyl('Cl,Br,I')
_AR_NAS_KETONE_NAPHTHYL_F = _ar_nas_naphthyl('F', _KETONE_EWG)
_AR_NAS_SULFONE_NAPHTHYL_F = _ar_nas_naphthyl('F', _SULFONE_EWG)
_AR_NAS_NITRILE_NAPHTHYL_F = _ar_nas_naphthyl('F', _NITRILE_EWG)
_AR_NAS_KETONE_NAPHTHYL_X = _ar_nas_naphthyl('Cl,Br,I', _KETONE_EWG)
_AR_NAS_SULFONE_NAPHTHYL_X = _ar_nas_naphthyl('Cl,Br,I', _SULFONE_EWG)
_AR_NAS_NITRILE_NAPHTHYL_X = _ar_nas_naphthyl('Cl,Br,I', _NITRILE_EWG)


# A cyclic anhydride is required: an acyclic anhydride must not be turned into
# an acyclic imide.  The recursive ring check accepts the five-membered
# succinic/maleic/phthalic anhydride ring and the six-membered glutaric ring.
# This is the net cyclodehydration step.  The open poly(amic acid) is not
# represented as a separate residue; the ring check ensures that only a cyclic
# anhydride can form the cyclic imide.
_ANHYDRIDE_PATH = '[#6,#7,#8,#16]'
_ANHYDRIDE = (
    '[#6&X3&R:1](=[#8:2])'
    '[#8&X2&R;'
    '$([#8]1[#6&X3](=[#8])'
    f'{_ANHYDRIDE_PATH}{_ANHYDRIDE_PATH}[#6&X3]1=[#8]),'
    '$([#8]1[#6&X3](=[#8])'
    f'{_ANHYDRIDE_PATH}{_ANHYDRIDE_PATH}{_ANHYDRIDE_PATH}'
    '[#6&X3]1=[#8])]'
    '[#6&X3&R:4](=[#8:5])'
)
_IMIDE = (
    _ANHYDRIDE + '.' + _amine(6, primary_only=True) +
    '>>[#6&X3:1](=[#8:2])[N&X3:6][#6&X3:4](=[#8:5])'
)


# Isocyanate + epoxide cycloaddition.  The isocyanate N attacks the less
# hindered epoxide carbon; the epoxide O stays on the more substituted carbon
# and becomes the ring O of the oxazolidinone.  The templates are ordered by
# selectivity: CH2 attack first, then CH attack, then the symmetric and
# fully-substituted fallbacks.
_ISOCYANATE = '[N&X2:4]=[C&X2:5]=[#8,#16;X1:6]'
_OXAZOLIDINONE_PRODUCT = (
    '[C:1]1[O:2][C:5](=[#8,#16:6])[N:4][C:3]1'
)


def _oxazolidinone(c1_h: str, c3_h: str) -> str:
    return (
        f'[C&X4;{c1_h};R:1]1[O&X2&R:2][C&X4;{c3_h};R:3]1.'
        + _ISOCYANATE + '>>' + _OXAZOLIDINONE_PRODUCT
    )


# Base-catalysed epoxide ROP: the nucleophile attacks the less hindered
# epoxide carbon; the epoxide O becomes the terminal alcohol on the more
# substituted carbon.  Water is allowed as the explicit initiator.
_EPOXIDE_PRODUCT = '[#8,#16;X1,X2:1][C&X4:2][C&X4:4][O&X2:3]'


def _epoxide_rop(c2_h: str, c4_h: str) -> str:
    return (
        _nucleophile(1, water=True) +
        f'.[C&X4;{c2_h};R:2]1[O&X2&R:3][C&X4;{c4_h};R:4]1>>'
        + _EPOXIDE_PRODUCT
    )


# Carbonate linkage by phosgene, the industrial route.  The carbonyl C and O
# come from phosgene and are mapped through, so the reaction is balanced
# (2 R-OH + Cl2C=O -> R-O-C(=O)-O-R + 2 HCl) and no oxygen is invented.  Two
# alcohol slots are needed because the coupling consumes both O-H groups of
# the carbonate.  A real process can also stop at a chloroformate or
# half-carbonate ester; the step-growth driver does not retain those transient
# end groups.
_CARBONATE = (
    _nucleophile(1) + '.[C&X3:3](=[O:4])([Cl,Br])[Cl,Br].' +
    _nucleophile(2) +
    '>>[#8,#16;X2:1][C&X3:3](=[O:4])[#8,#16;X2:2]'
)


STEP_GROWTH_COUPLINGS = {
    # Direct esterification/acylation (acid or acyl halide) first, then
    # methyl-ester transesterification (the DMT/industrial polyester route).
    # The methyl-ester template is deliberately methyl-specific: a general
    # `[O&X2]` ester query would also match backbone esters and would allow
    # the chain to cleave itself.
    'ester': [
        (
            '[C&X3:1](=[O:6])([#6:3])[O&X2&H1,Cl,Br,I].' +
            _nucleophile(2) +
            '>>[C&X3:1](=[O:6])([#6:3])[O,S&X2:2]'
        ),
        (
            '[C&X3:1](=[O:6])([#6:3])[O&X2][CX4;H3:7].' +
            _nucleophile(2) +
            '>>[C&X3:1](=[O:6])([#6:3])[O,S&X2:2]'
        ),
    ],
    'amide': (
        '[C&X3:1](=[O:6])([#6:3])[O&X2&H1,Cl,Br,I].' +
        _amine(2) +
        '>>[C&X3:1](=[O:6])([#6:3])[N&X3:2]'
    ),
    'urethane': (
        '[N&X2:1]=[C&X2:2]=[#8,#16;X1:3].' +
        _nucleophile(4) +
        '>>[N&X3:1][C&X3:2](=[#8,#16;X1:3])[O,S&X2:4]'
    ),
    'aryl_ether_nas': [
        _AR_NAS_KETONE_F,
        _AR_NAS_SULFONE_F,
        _AR_NAS_NITRILE_F,
        _AR_NAS_NITRILE_ORTHO_F,
        _AR_NAS_KETONE_BIPHENYL_F,
        _AR_NAS_SULFONE_BIPHENYL_F,
        _AR_NAS_NITRILE_BIPHENYL_F,
        _AR_NAS_KETONE_NAPHTHYL_F,
        _AR_NAS_SULFONE_NAPHTHYL_F,
        _AR_NAS_NITRILE_NAPHTHYL_F,
        _AR_NAS_KETONE_X,
        _AR_NAS_SULFONE_X,
        _AR_NAS_NITRILE_X,
        _AR_NAS_NITRILE_ORTHO_X,
        _AR_NAS_KETONE_BIPHENYL_X,
        _AR_NAS_SULFONE_BIPHENYL_X,
        _AR_NAS_NITRILE_BIPHENYL_X,
        _AR_NAS_KETONE_NAPHTHYL_X,
        _AR_NAS_SULFONE_NAPHTHYL_X,
        _AR_NAS_NITRILE_NAPHTHYL_X,
    ],
    'imide': _IMIDE,
    'oxazolidinone': [
        _oxazolidinone('H0,H1', 'H2'),
        _oxazolidinone('H0,H1', 'H1'),
        _oxazolidinone('H2', 'H2'),
        _oxazolidinone('H0', 'H0'),
    ],
    'epoxide_rop': [
        _epoxide_rop('H2', 'H0,H1'),
        _epoxide_rop('H1', 'H0,H1'),
        _epoxide_rop('H2', 'H2'),
        _epoxide_rop('H0', 'H0'),
    ],
    'carbonate': _CARBONATE,
    # Oxidative phenol coupling (poly(phenylene oxide)): the phenolic O-H of
    # one ring couples to the para C-H of another ring.  Both ortho positions
    # of the first ring must be blocked (H0), otherwise unsubstituted phenol
    # would be accepted and would branch through its free ortho positions.
    # The second ring's O-H and the first ring's para C-H remain real terminal
    # groups.
    'phenol_oxidative': (
        '[c:1]1([O&X2&H1:2])[cH0:10][c:11][c:12][c:13][cH0:14]1.'
        '[c:3]1([O&X2:4])[c:5][c:6][c&H1:7][c:8][c:9]1>>'
        '[c:1]1([O&X2:2][c:7]2[c:6][c:5][c:3]([O&X2:4])[c:9][c:8]2)'
        '[c:10][c:11][c:12][c:13][c:14]1'
    ),
}


# ROP/ROCOP propagation templates.  These couple the opened monomer's real
# end groups.  They are deliberately narrower than the generic step-growth
# couplings: they relax only the guard that is wrong for the opened end group.
RING_OPENING_PROPAGATIONS = {
    # A lactone opens to a hydroxyacid.  The electrophile is the carboxylic
    # acid carbonyl (its carbon has a carbon neighbour, so a carbamic acid
    # carbonyl attached to N is not selected).  The nucleophile may be an
    # alcohol or an acyloxy O-H attached to a carbamic carbonyl, but not a
    # carboxylic acid/ester/carbonate or heteroatom-acid O-H.
    'lactone_rop': [
        (
            '[C&X3:1](=[O:6])([#6:3])[O&X2&H1,Cl,Br,I].' +
            _nucleophile(2) +
            '>>[C&X3:1](=[O:6])([#6:3])[O,S&X2:2]'
        ),
        (
            '[C&X3:1](=[O:6])([#6:3])[O&X2&H1,Cl,Br,I].' +
            _carbamic_oh(2) +
            '>>[C&X3:1](=[O:6])([#6:3])[O&X2:2]'
        ),
    ],
    # A lactam opens to an aminoacid.  The amine may be attached to an
    # amidine/thioimidate carbon (C=N), but not to an amide/imide carbonyl
    # or thiocarbonyl (C=O/S), sulfonyl, or phosphoryl group.
    'lactam_rop': (
        '[C&X3:1](=[O:6])([#6:3])[O&X2&H1,Cl,Br,I].' +
        _amine(2) +
        '>>[C&X3:1](=[O:6])([#6:3])[N&X3:2]'
    ),
    # ROCOP of an anhydride + epoxide opens to a hydroxyacid with an ester
    # linkage.  Ordinary ester coupling then inserts the next opened pair.
    'anhydride_epoxide_rop': STEP_GROWTH_COUPLINGS['ester'],
}


# Balanced water-initiated ROCOP opening.  The cyclic anhydride and epoxide
# alone contain four oxygens for the succinic anhydride/ethylene oxide example;
# the opened hydroxyacid has five because one water molecule is consumed.  The
# anhydride ring O becomes the ester O, the water O becomes the remaining acid
# O-H, and the epoxide O becomes the terminal alcohol.  Mapping all three
# oxygen sources keeps the reaction mass-balanced (unlike the topological
# dataset key, which only describes the repeat unit).
ROCOP_OPENING_SMARTS = (
    '[#6,#7,#8,#16;R;r5,r6:1]'
    '[C&X3,c;R:10](=[O&X1:11])[O&X2,o;R;r5,r6:3]'
    '[C&X3,c;R:12](=[O&X1:13])[#6,#7,#8,#16;R;r5,r6:8].'
    '[C&X4&R:4]1[O&X2&R:5][C&X4&R:6]1.'
    '[O&X2;H2:9]>>'
    '([#6,#7,#8,#16:1][C&X3:10](=[O&X1:11])[O&X2:9].'
    '[#6,#7,#8,#16:8][C&X3:12](=[O&X1:13])[O&X2:3]'
    '[C&X4:4][C&X4:6][O&X2:5])'
)


def _couplings(*names):
    """Return the coupling reactions for the named chemistries.

    A chemistry may map to one SMARTS string or to a list of alternatives
    (for example, regioselective epoxide templates).  The returned list is
    always flat so the drivers receive a list of strings.
    """
    reactions = []
    for name in names:
        value = STEP_GROWTH_COUPLINGS[name]
        reactions.extend(value if isinstance(value, list) else [value])
    return reactions


def _propagations(*names):
    """Return the ring-opening propagation reactions for *names*."""
    reactions = []
    for name in names:
        value = RING_OPENING_PROPAGATIONS[name]
        reactions.extend(value if isinstance(value, list) else [value])
    return reactions


def reaction_config(key):
    """Return driver kwargs for a configured ``rxn_map`` entry.

    ``ring_open``/``ring_opening`` entries use the map key itself as the
    ring-opening SMARTS, and ``romp`` entries use it as the ring-bond SMARTS,
    so the key is not duplicated in the configuration dictionary.
    """
    config = rxn_map[key]
    if config is None:
        raise ValueError(f"No extender recipe configured for {key!r}")
    params = dict(config["params"])
    if params.pop("ring_open", False):
        params["ring_open_smarts"] = key
    if config.get("type") == "ring_opening":
        params.setdefault("ring_open_smarts", key)
    if config.get("type") == "romp":
        params.setdefault("ring_bond_smarts", key)
    return params


# Terminal-group policy
# ---------------------
# * addition entries: the polymer end groups come from the chain-growth
#   init_smarts/term_smarts.  The wildcards in the repeat-unit key are only
#   main-chain cut points; they are not end groups.
# * step_growth entries: a coupling reaction consumes one functional group
#   from the chain and one from a monomer.  RDKit carries the rest of the
#   monomer through unchanged, so the polymer keeps the monomer's exact
#   unreacted end groups (-COOH, -COCl, -OH, -SH, -NH2, -N=C=O, -N=C=S,
#   Ar-F/Ar-Cl, epoxide, anhydride, phenol, ...).  No end group is synthesized
#   from a wildcard.
# * ring_open/ring_opening entries: ring opening infers only the real open
#   hydroxyacid/aminoacid (and, for ROCOP, anhydride + epoxide hydroxyacid)
#   form of the cyclic monomer.  The free ends then come from that opened
#   monomer.  The open-form conversion is a chemical reaction, not a
#   wildcard-derived polymer cap; it is covered by tests, including the former
#   thioimidate-lactam and cyclic-carbamate/lactone gaps.
# * romp entries: the ring is cut programmatically at the map key's C=C
#   double bond (the metathesis event), not at the allylic single bond.  The
#   initiator and terminator fragments are explicit recipe parameters; ring
#   opening does not synthesize a heavy atom.
# * entries mapped to None would be left unconfigured rather than guessed.  The
#   registry has no such entries.

rxn_map = {
    # -- polyolefin --
    # Addition/chain-growth: initiation and termination are hydrogen caps, so
    # the chain has the monomer's own heavy atoms plus no extra end carbon.
    # The R15/R16 keys are comonomer pairs: the propagation template inserts
    # one unit of each per step, giving an alternating ABAB model, not a
    # conjugated-diene 1,4-addition.
    '[C&X3&R:1]=[C&X3&R:2]>>*-[C&X4&R:1][C&X4&R:2]-*': {
        "polymer_class": 'polyolefin',
        "mechanism": "A2",
        "type": "addition",
        "params": {
            'init_smarts': '[C&X3&R:1]=[C&X3&R:2]>>[C&X4&R:1][C&X4&R:2]-*',
            'prop_smarts': '[C:0][C&X4:1]-*.[C&X3&R:2]=[C&X3&R:3]>>[C:0][C:1][C&X4&R:2][C&X4&R:3]-*',
            'term_smarts': '[C:0][C:1]-*>>[C:0][C:1]',
        },
    },
    '[C&X3;H2,H1,H0:1]=[C;H2,H1,H0:2]>>*-[C:1][C:2]-*': {
        "polymer_class": 'polyolefin',
        "mechanism": "A1",
        "type": "addition",
        "params": {
            'init_smarts': '[C&X3;H2:1]=;!@[C&X3;H2,H1,H0:2]>>[C:1][C:2]-*',
            'prop_smarts': '[C:0][C&X4:1]-*.[C&X3;H2:2]=;!@[C&X3;H2,H1,H0:3]>>[C:0][C:1][C:2][C:3]-*',
            'term_smarts': '[C:0][C:1]-*>>[C:0][C:1]',
        },
    },
    '[C&X3;H2,H1,H0;!R:1]=[C&X3;H2,H1,H0;!R:2].[C&X3;H1,H0;R:3]=[C&X3;H1,H0;R:4]>>*-[C&X4:1][C&X4:2][C&X4:3][C&X4:4]-*': {
        "polymer_class": 'polyolefin',
        "mechanism": "A3",
        "type": "addition",
        "params": {
            'init_smarts': '[C&X3;H2:1]=;!@[C&X3;H2,H1,H0:2].[C&X3;H1,H0;R:3]=[C&X3;H1,H0;R:4]>>[C&X4:1][C&X4:2][C&X4:3][C&X4:4]-*',
            'prop_smarts': '[C:0][C&X4:1]-*.[C&X3;H2:2]=;!@[C&X3;H2,H1,H0:3].[C&X3;H1,H0;R:4]=[C&X3;H1,H0;R:5]>>[C:0][C:1][C&X4:2][C&X4:3][C&X4:4][C&X4:5]-*',
            'term_smarts': '[C:0][C:1]-*>>[C:0][C:1]',
        },
    },
    '[C&X3;H2,H1,H0;!R:1]=[C&X3;H2,H1,H0;!R:2].[C&X3;H2,H1,H0;!R:3]=[C&X3;H2,H1,H0;!R:4]>>*-[C&X4:1][C&X4:2][C&X4:3][C&X4:4]-*': {
        "polymer_class": 'polyolefin',
        "mechanism": "A4",
        "type": "addition",
        "params": {
            'init_smarts': '[C&X3;H2:1]=;!@[C&X3;H2,H1,H0:2].[C&X3;H2:3]=;!@[C&X3;H2,H1,H0:4]>>[C&X4:1][C&X4:2][C&X4:3][C&X4:4]-*',
            'prop_smarts': '[C:0][C&X4:1]-*.[C&X3;H2:2]=;!@[C&X3;H2,H1,H0:3].[C&X3;H2:4]=;!@[C&X3;H2,H1,H0:5]>>[C:0][C:1][C&X4:2][C&X4:3][C&X4:4][C&X4:5]-*',
            'term_smarts': '[C:0][C:1]-*>>[C:0][C:1]',
        },
    },
    # -- polyester --
    # Ester coupling: the unreacted acid (-COOH/-COCl) and alcohol/thiol
    # (-OH/-SH) ends are carried over from the monomers.  The compiled
    # coupling also has a methyl-ester transesterification tier for the
    # industrial DMT route.
    '([C&X3:1](=O)[O&X2&H1,Cl,Br].[C&X3:2](=O)[O&X2&H1,Cl,Br]).([O,S;X2;H1;!$([O,S]C=*):3].[O,S;X2;H1;!$([O,S]C=*):4])>>(*-[C&X3:1]=O.[C&X3:2](=O)-[O,S;X2;!$([O,S]C=*):3].[O,S;X2;!$([O,S]C=*):4]-*)': {
        "polymer_class": 'polyester',
        "mechanism": "C1",
        "type": "step_growth",
        "params": {
            "coupling_smarts": _couplings('ester'),
        },
    },
    # Ester coupling between two hydroxyacids: the free -OH and -COOH ends
    # remain the hydroxyacid monomers' real terminal groups.
    '([O&X2&H1&!$(OC=*):1].[C&X3:2](=O)[O&X2&H1]).([O&X2&H1&!$(OC=*):3].[C&X3:4](=O)[O&X2&H1])>>(*-[O&X2:1].[C&X3:2](=O)[O&X2:3].[C&X3:4](=O)-*)': {
        "polymer_class": 'polyester',
        "mechanism": "C1",
        "type": "step_growth",
        "params": {
            "coupling_smarts": _couplings('ester'),
        },
    },
    # Ester coupling of a single hydroxyacid: the terminal -OH and -COOH
    # are the monomer's own groups.
    '([O&X2&H1&!$(OC=*):1].[C&X3:2](=O)[O&X2&H1])>>(*-[O&X2:1].[C&X3:2](=O)-*)': {
        "polymer_class": 'polyester',
        "mechanism": "C1",
        "type": "step_growth",
        "params": {
            "coupling_smarts": _couplings('ester'),
        },
    },
    # Carbonate linkage by phosgene.  The dataset key describes the CO-based
    # repeat-unit topology, but a direct 2 R-OH + CO condensation would lose H2
    # without an oxidant.  The configured coupling therefore uses phosgene, the
    # industrial carbonate route: 2 R-OH + Cl2C=O -> R-O-C(=O)-O-R + 2 HCl.
    '([O,S;X2;H1;!$([O,S]C=*):1].[O,S;X2;H1;!$([O,S]C=*):2]).[C&-]#[O&+]>>(*-[O,S;X2;!$([O,S]C=*):1].[O,S;X2;!$([O,S]C=*):2][C&X3](=O)-*)': {
        "polymer_class": 'polyester',
        "mechanism": "C6",
        "type": "step_growth",
        "params": {
            "coupling_smarts": _couplings('carbonate'),
            "reagents": ['ClC(=O)Cl'],
        },
    },
    # Lactone ring-opening polymerisation (ROP).  The map key opens the
    # lactone to a real hydroxyacid; the ROP propagation template then inserts
    # an opened hydroxyacid.  The first residue is water-initiated and has one
    # more oxygen than the cyclic monomer; later residues add only the
    # monomer's own atoms.
    '[C&X3&R:1](=[O&X1])[O&X2&R:2]>>(*-[C&X3:1]=[O&X1].[O&X2:2]-*)': {
        "polymer_class": 'polyester',
        "mechanism": "R1",
        "type": "ring_opening",
        "params": {
            "prop_smarts": _propagations('lactone_rop'),
        },
    },
    # Anhydride + epoxide ring-opening copolymerisation (ROCOP).  The balanced
    # opening consumes one water molecule; the ring O becomes the ester O, the
    # water O becomes the acid O-H, and the epoxide O becomes the terminal
    # alcohol.  The opened unit is then propagated with ester coupling.  The
    # anhydride ring is five- or six-membered and its flanking atoms may be
    # C, N, O, or S.
    '[C,c;R:1][C&X3,c;R](=[O&X1])[O&X2,o;R][C&X3,c;R](=[O&X1])[C,c;R:2].[C&X4&R:3]1[O&X2&R:4][C&X4&R:5]1>>([C,c:1][C&X3](=[O&X1])-*.[C,c:2][C&X3](=[O&X1])[O&X2][C&X4:3][C&X4:5][O&X2:4]-*)': {
        "polymer_class": 'polyester',
        "mechanism": "R3",
        "type": "ring_opening",
        "params": {
            "prop_smarts": _propagations('anhydride_epoxide_rop'),
            "ring_open_smarts": ROCOP_OPENING_SMARTS,
            "opening_reagents": ['O'],
        },
    },
    # -- polyether --
    # Epoxide ROP: the epoxide monomer has no pre-existing terminal group.
    # Water is supplied as the initiator reagent; the -OH ends come from
    # water and the opened epoxide oxygen.  This initiator choice is explicit,
    # not a wildcard-derived cap.
    '[C&X4;H2,H1,H0;R:1]1[O&R][C&R:2]1>>*-[C&X4:1][C&X4:2]O-*': {
        "polymer_class": 'polyether',
        "mechanism": "R4",
        "type": "step_growth",
        "params": {
            "coupling_smarts": _couplings('epoxide_rop'),
            "reagents": ['O'],
        },
    },
    # Aryl-ether coupling (PEEK/PES/PAEN NAS): one activated Ar-X reacts; the
    # other Ar-X and the phenol -OH remain the monomers' exact terminal groups.
    # The compiled coupling covers ketone, sulfone, and nitrile activation in
    # same-ring, phenyl-spaced, and naphthyl-spaced topologies, with F before
    # Cl/Br/I.
    '[c:1]1[c:2][c:3](F)[c:4][c:5][c:6]1[C&X3](=[O&X1])[c:7]1[c:8][c:9][c:10](F)[c:11][c:12]1.([O&X2&H1&!$([O,S]C=*):13].[O&X2&H1&!$([O,S]C=*):14])>>([c:1]1[c:2][c:3](-*)[c:4][c:5][c:6]1[C&X3](=[O&X1])[c:7]1[c:8][c:9][c:10]([O&X2&!$([O,S]C=*):13])[c:11][c:12]1.[O&X2&!$([O,S]C=*):14]-*)': {
        "polymer_class": 'polyether',
        "mechanism": "C7",
        "type": "step_growth",
        "params": {
            "coupling_smarts": _couplings('aryl_ether_nas'),
        },
    },
    # Aryl-ether coupling (PES NAS): one Ar-F/Ar-Cl reacts; the other
    # Ar-F/Ar-Cl and the phenol -OH remain exact monomer ends.  The compiled
    # coupling is shared with the PEEK/PAEN templates above.
    '[c:1]1[c:2][c:3]([F,Cl,Br,I])[c:4][c:5][c:6]1[S&X4](=[O&X1])(=[O&X1])[c:7]1[c:8][c:9][c:10]([F,Cl,Br,I])[c:11][c:12]1.([O&X2&H1&!$([O,S]C=*):13].[O&X2&H1&!$([O,S]C=*):14])>>([c:1]1[c:2][c:3](-*)[c:4][c:5][c:6]1[S&X4](=[O&X1])(=[O&X1])[c:7]1[c:8][c:9][c:10]([O&X2&!$([O,S]C=*):13])[c:11][c:12]1.[O&X2&!$([O,S]C=*):14]-*)': {
        "polymer_class": 'polyether',
        "mechanism": "C7",
        "type": "step_growth",
        "params": {
            "coupling_smarts": _couplings('aryl_ether_nas'),
        },
    },
    # Oxidative phenol coupling (poly(phenylene oxide)): the phenolic O-H of
    # one ring couples to the para C-H of another ring.  The compiled coupling
    # requires both ortho positions of the first ring to be blocked.  The
    # second ring's O-H is carried over as the real terminal group; no end
    # group is guessed.
    'c1([O&H1:1])[c:2][c:3][c&H1:4][c:5][c:6]1>>c1([O&X2:1]-*)[c:2][c:3][c:4](-*)[c:5][c:6]1': {
        "polymer_class": 'polyether',
        "mechanism": "C8",
        "type": "step_growth",
        "params": {
            "coupling_smarts": _couplings('phenol_oxidative'),
        },
    },
    # -- polyamide --
    # Amide coupling: unreacted -COOH/-COCl and -NH2 ends are carried over
    # from the monomers.
    '([C&X3:1](=O)[O&X2&H1,Cl,Br].[C&X3:2](=O)[O&X2&H1,Cl,Br]).([N&X3;H2,H1;!$(NC=*):3].[N&X3;H2,H1;!$(NC=*):4])>>(*-[C&X3:1]=O.[C&X3:2](=O)-[N&X3&!$(NC=*):3].[N&X3&!$(NC=*):4]-*)': {
        "polymer_class": 'polyamide',
        "mechanism": "C2",
        "type": "step_growth",
        "params": {
            "coupling_smarts": _couplings('amide'),
        },
    },
    # Amide coupling of two aminoacids: the free -NH2 and -COOH ends remain
    # the aminoacid monomers' real terminal groups.
    '([N&X3;H2,H1;!$(NC=*):1].[C&X3:2](=O)[O&X2&H1]).([N&X3;H2,H1;!$(NC=*):3].[C&X3:4](=O)[O&X2&H1])>>(*-[N&X3&!$(NC=*):1].[C&X3:2](=O)[N&X3&!$(NC=*):3].[C&X3:4](=O)-*)': {
        "polymer_class": 'polyamide',
        "mechanism": "C2",
        "type": "step_growth",
        "params": {
            "coupling_smarts": _couplings('amide'),
        },
    },
    # Amide coupling of a single aminoacid: the terminal -NH2 and -COOH
    # are the monomer's own groups.
    '([N&X3;H2,H1;!$(OC=*):1].[C&X3:2](=O)[O&X2&H1])>>(*-[N&X3:1].[C&X3:2](=O)-*)': {
        "polymer_class": 'polyamide',
        "mechanism": "C2",
        "type": "step_growth",
        "params": {
            "coupling_smarts": _couplings('amide'),
        },
    },
    # Lactam ring-opening polymerisation (ROP).  The map key opens the
    # lactam to a real aminoacid; the ROP propagation template then inserts an
    # opened aminoacid.  The first residue is water-initiated and has one more
    # oxygen than the cyclic lactam; later residues add only the monomer's own
    # atoms.
    '[C&X3&R:1](=[O&X1])[N&X3&R:2]>>(*-[C&X3:1]=[O&X1].[N&X3:2]-*)': {
        "polymer_class": 'polyamide',
        "mechanism": "R2",
        "type": "ring_opening",
        "params": {
            "prop_smarts": _propagations('lactam_rop'),
        },
    },
    # -- polyimide --
    # Imide coupling: one anhydride + one amine react; the other amine
    # (-NH2) and the unreacted anhydride remain exact monomer ends.
    '([C&X3,c;R:1](=[O&X1])[O&X2,o;R][C&X3,c;R:2]=[O&X1].[C&X3,c;R:3](=[O&X1])[O&X2,o;R][C&X3,c;R:4]=[O&X1]).([C,c:5][N&X3&H2&!$(N[C,S]=*)].[C,c:6][N&X3&H2&!$(N[C,S]=*)])>>([C&X3,c;R:1](=[O&X1])[N&X3&R]([C,c:5])[C&X3,c;R:2]=[O&X1].[C,c:6]-*.[C&X3,c;R:3](=[O&X1])[N&X3&R](-*)[C&X3&R:4]=[O&X1])': {
        "polymer_class": 'polyimide',
        "mechanism": "C4",
        "type": "step_growth",
        "params": {
            "coupling_smarts": _couplings('imide'),
        },
    },
    # -- polyurethane --
    # Urethane/thiourethane coupling: one -N=C=O/-N=C=S reacts; the other
    # isocyanate/isothiocyanate and the diol/thiol -OH/-SH remain exact
    # monomer ends.
    '([N&X2:1]=[C&X2]=[O&X1,S&X1:2].[N&X2:3]=[C&X2:4]=[O&X1,S&X1:5]).([O&X2,S&X2;H1;!$([O,S]C=*):6].[O&X2,S&X2;H1;!$([O,S]C=*):7])>>(*-[C&X3](=[O&X1,S&X1:2])[N&X3:1].[N&X3:3][C&X3:4](=[O&X1,S&X1:5])[O&X2,S&X2;!$([O,S]C=*):6].[O&X2,S&X2;!$([O,S]C=*):7]-*)': {
        "polymer_class": 'polyurethane',
        "mechanism": "C3",
        "type": "step_growth",
        "params": {
            "coupling_smarts": _couplings('urethane'),
        },
    },
    # -- polyoxazolidone --
    # Oxazolidinone coupling: one epoxide + one isocyanate react; the other
    # epoxide and -N=C=O remain exact monomer ends.
    '([C&X4;H2,H1,H0;R:1]1[O&X2&R:2][C&X4;H1,H0;R:3]1.[C&X4;H2,H1,H0;R:4]1[O&X2&R:5][C&X4;H1,H0;R:6]1).([O&X1,S&X1:7]=[C&X2:8]=[N&X2:9].[O&X1,S&X1:10]=[C&X2:11]=[N&X2:12][C,c:13])>>([C&X4&R:6]1[O&X2&R:5][C&X2&R:8](=[O&X1,S&X1:7])[N&X3&R:9][C&X4&R:4]1.[C&X4&R:3]1[O&X2&R:2][C&X2&R:11](=[O&X1,S&X1:10])[N&X3&R:12](-*)[C&X4&R:1]1.[C,c:13]-*)': {
        "polymer_class": 'polyoxazolidone',
        "mechanism": "C5",
        "type": "step_growth",
        "params": {
            "coupling_smarts": _couplings('oxazolidinone'),
        },
    },
    # ROMP is chain-growth metathesis.  The reactant side names the ring
    # alkene and an allylic atom so the variable-size ring can be found; the
    # driver cuts the C=C double bond (not the allylic single bond) and joins
    # repeat units through new double bonds.  The product side is therefore
    # written as the two double-bonded cut ends, not as a retained ring alkene.
    '[C&X3&H1&R:1]=[C&X3&H1&R:2][C&X4&R:3]>>(*=[C&X3&H1:1].[C&X3&H1:2]=*)': {
        "polymer_class": 'polyolefin',
        "mechanism": "R5",
        "type": "romp",
        "params": {
            "initiator": 'C',
            "terminator": 'C',
        },
    },
}

# Keep the recorded dataset key as a compatibility alias.  The generalized
# key exposes the non-ring guard directly and is the one used by the
# practical-example manifest; both entries share the same compiled recipe.
rxn_map[
    '[C&X3;H2:1]=;!@[C&X3;H2,H1,H0:2]>>*-[C:1][C:2]-*'
] = rxn_map[
    '[C&X3;H2,H1,H0:1]=[C;H2,H1,H0:2]>>*-[C:1][C:2]-*'
]
