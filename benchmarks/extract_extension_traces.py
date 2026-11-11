"""Extract and document short extension traces for every driver family.

The conventional page uses well-known polymers that are easy to misrepresent
if repeat-unit wildcards are simply joined:

* polypropylene (PP) for chain growth, using the generic vinyl addition
  template with propene;
* poly(ethylene terephthalate) (PET) for step growth, using the polyester
  AABB template with terephthalic acid and ethylene glycol;
* poly(ethylene oxide) (PEO) for epoxide ring-opening propagation, using
  ethylene oxide and water;
* polycaprolactone (PCL) for lactone ring-opening polymerisation;
* Nylon-6 for lactam ring-opening polymerisation;
* poly(ethylene succinate) (PES) for anhydride + epoxide ring-opening
  copolymerisation (ROCOP);
* poly(2,6-dimethyl-1,4-phenylene oxide) (PPO) for oxidative phenol coupling;
* poly(norbornene) (PNB) for ring-opening metathesis polymerisation (ROMP).

Each trace is short enough to inspect by eye but long enough to show the
enchainment, the real end groups, and the difference between the mechanisms.

Other modes:

* ``--mode simple`` uses propene and glycolic acid as minimal template
  examples;
* ``--mode shortest`` uses the shortest surviving candidates from
  ``tests/test_data/extender_test_data_all.csv``;
* ``--mode mechanisms`` emits one persistent example for each of the 17
  mechanism IDs (A1–R5) in ``mechanism_table.tex``, using a known practical
  polymer for each mechanism.  The mechanism code is stored on the
  corresponding ``rxn_map`` entry, and the generated
  ``mechanism_examples.md`` records the exact reaction SMARTS, practical
  polymer, monomers, structures and image paths.  A mechanism with no
  identified practical polymer is marked ``needs candidate`` instead of
  inventing one.

Usage::

    python benchmarks/extract_extension_traces.py
    python benchmarks/extract_extension_traces.py --mode simple
    python benchmarks/extract_extension_traces.py --mode shortest
    python benchmarks/extract_extension_traces.py --mode mechanisms

Outputs::

    benchmarks/extension_traces.md
    benchmarks/traces/{chain,step,peo,pcl,nylon6,pes,ppo,romp}_*.svg
    benchmarks/mechanism_examples.md          (mechanisms mode)
    benchmarks/traces/mechanisms/practical/*.svg (mechanisms mode)
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem.Draw import MolsToGridImage, rdMolDraw2D

_PROJECT = Path(__file__).resolve().parent.parent
if str(_PROJECT) not in sys.path:
    sys.path.insert(0, str(_PROJECT))

from polyscript.extender import (  # noqa: E402
    ChainGrowthReaction,
    RingOpeningMetathesisReaction,
    RingOpeningReaction,
    StepGrowthReaction,
)
from polyscript.extender.map import reaction_config, rxn_map  # noqa: E402

RDLogger.DisableLog("rdApp.*")

_CSV = _PROJECT / "tests" / "test_data" / "extender_test_data_all.csv"
_BENCH = _PROJECT / "benchmarks"
_TRACES = _BENCH / "traces"
_MECH_TRACES = _TRACES / "mechanisms" / "practical"
_MD = _BENCH / "extension_traces.md"
_MECH_MD = _BENCH / "mechanism_examples.md"
_PROPAGATION_STEPS = 2

_VINYL = "[C&X3;H2,H1,H0:1]=[C;H2,H1,H0:2]>>*-[C:1][C:2]-*"
_GLYCOLIC_AB = (
    "([O&X2&H1&!$(OC=*):1].[C&X3:2](=O)[O&X2&H1])>>"
    "(*-[O&X2:1].[C&X3:2](=O)-*)"
)
_PET_AABB = (
    "([C&X3:1](=O)[O&X2&H1,Cl,Br].[C&X3:2](=O)[O&X2&H1,Cl,Br])."
    "([O,S;X2;H1;!$([O,S]C=*):3].[O,S;X2;H1;!$([O,S]C=*):4])>>"
    "(*-[C&X3:1]=O.[C&X3:2](=O)-[O,S;X2;!$([O,S]C=*):3]."
    "[O,S;X2;!$([O,S]C=*):4]-*)"
)

_CONVENTIONAL = {
    "addition": {
        "reaction": _VINYL,
        "monomers": ["C=CC"],
        "polymer_name": "Polypropylene (PP)",
        "challenge": (
            "The PP repeat unit has a CH2 (head) and a CH(CH3) (tail). "
            "Joining the two wildcard ends blindly can give head-to-head or "
            "tail-to-tail links instead of the head-to-tail enchainment "
            "produced by the vinyl template. It also discards the "
            "initiator/terminator-derived end groups and cannot encode "
            "tacticity. The trace keeps the head-to-tail repeat and ends "
            "with hydrogen caps."
        ),
    },
    "step_growth": {
        "reaction": _PET_AABB,
        "monomers": ["O=C(O)c1ccc(C(=O)O)cc1", "OCCO"],
        "polymer_name": "Poly(ethylene terephthalate) (PET)",
        "challenge": (
            "The PET repeat unit has two chemically different cut points: "
            "the alcohol-side O and the acid-side C(=O). Joining those "
            "wildcards without tracking which is which can invert the ester "
            "direction or create a non-ester link, and it drops the real "
            "-OH/-COOH ends and the water byproduct."
        ),
        "sequence_note": (
            "The trace is diacid-diol-diacid-diol: both monomers appear "
            "twice, and the ends are one -COOH and one -OH. The three-residue "
            "intermediate (diol-diacid-diol) has two -OH ends because its "
            "acid groups are consumed."
        ),
    },
}


def _key(prefix: str, rxn_type: str | None = None) -> str:
    """Return the unique ``rxn_map`` key starting with *prefix*."""
    matches = [
        key
        for key, config in rxn_map.items()
        if config
        and key.startswith(prefix)
        and (rxn_type is None or config.get("type") == rxn_type)
    ]
    if len(matches) != 1:
        raise KeyError((prefix, rxn_type, matches))
    return matches[0]


_EPOXIDE_ROP = _key("[C&X4;H2,H1,H0;R:1]1[O&R]")
_LACTONE_ROP = _key(
    "[C&X3&R:1](=[O&X1])[O&X2&R:2]", "ring_opening"
)
_LACTAM_ROP = _key(
    "[C&X3&R:1](=[O&X1])[N&X3&R:2]", "ring_opening"
)
_ROCOP = _key("[C,c;R:1]", "ring_opening")
_PHENOL_OXIDATIVE = _key("c1([O&H1:1])")
_ROMP = _key("[C&X3&H1&R:1]", "romp")


# One persistent example per mechanism ID in ``mechanism_table.tex``.  The
# ``reaction`` values are resolved from ``rxn_map`` (never copied by hand) and
# each registry entry carries the same ``mechanism`` code, so the table, the
# registry and the generated illustration page cannot drift apart.  Each spec
# names a known practical polymer and its monomers; if no practical example is
# identified, use ``"status": "needs candidate"`` and leave ``monomers``
# empty so the page leaves a clear placeholder.
_MECHANISMS = [
    {
        "id": "A1",
        "family": "Addition",
        "recipe": "single non-ring vinyl (CH2 attachment; H caps)",
        "reaction": _key(
            "[C&X3;H2:1]=;!@[C&X3;H2,H1,H0:2]", "addition"
        ),
        "monomers": ["C=CC"],
        "length": 2,
        "polymer": "Polypropylene (PP)",
        "status": "commercial",
        "note": (
            "Commodity polyolefin; each propagation inserts one propene "
            "unit, and termination leaves the hydrogen caps."
        ),
    },
    {
        "id": "A2",
        "family": "Addition",
        "recipe": "single cyclic vinyl (H caps)",
        "reaction": _key(
            "[C&X3&R:1]=[C&X3&R:2]", "addition"
        ),
        "monomers": ["C1CC2C=CC1C2"],
        "length": 2,
        "polymer": "Vinyl-addition poly(norbornene)",
        "status": "known/specialty",
        "prefix": "A2_norbornene",
        "note": (
            "Vinyl-addition polymerisation retains the norbornane ring, "
            "giving the saturated-backbone cyclic olefin polymer; this is "
            "distinct from the ROMP product in R5. The bridged ring in the "
            "images is norbornane, not cyclopentane."
        ),
    },
    {
        "id": "A3",
        "family": "Addition",
        "recipe": "acyclic + cyclic pair (CH2 attachment; H caps)",
        "reaction": _key(
            "[C&X3;H2,H1,H0;!R:1]=[C&X3;H2,H1,H0;!R:2]."
            "[C&X3;H1,H0;R:3]",
            "addition",
        ),
        "monomers": ["C=C", "C1CC2C=CC1C2"],
        "length": 1,
        "polymer": "Ethylene–norbornene cyclic olefin copolymer (COC)",
        "status": "commercial",
        "note": (
            "Commercial COC (TOPAS-type). This topological trace models "
            "the comonomer pair as alternating insertion, consuming one "
            "ethylene and one cyclic olefin per propagation step."
        ),
    },
    {
        "id": "A4",
        "family": "Addition",
        "recipe": "two non-ring vinyl pair (CH2 attachment; H caps)",
        "reaction": _key(
            "[C&X3;H2,H1,H0;!R:1]=[C&X3;H2,H1,H0;!R:2]."
            "[C&X3;H2,H1,H0;!R:3]",
            "addition",
        ),
        "monomers": ["C=C", "C=CC"],
        "length": 1,
        "polymer": "Ethylene–propylene rubber (EPR)",
        "status": "commercial",
        "paired": True,
        "note": (
            "Commercial EPR/EPDM backbone. The paired trace models the "
            "comonomer pair as alternating insertion, forcing one ethylene "
            "and one propylene into each propagation step; the default "
            "frontier cap alone would pick the same monomer twice."
        ),
    },
    {
        "id": "C1",
        "family": "Coupling",
        "recipe": (
            "AA+BB; HA; acid/acyl halide or methyl-ester "
            "transesterification"
        ),
        "reaction": _key(
            "([C&X3:1](=O)[O&X2&H1,Cl,Br]."
            "[C&X3:2](=O)[O&X2&H1,Cl,Br]).([O,S;X2;H1;",
            "step_growth",
        ),
        "monomers": [
            "O=C(O)c1ccc(C(=O)O)cc1",
            "OCCO",
        ],
        "length": 2,
        "polymer": "Poly(ethylene terephthalate) (PET)",
        "status": "commercial",
        "note": (
            "Commodity polyester; the acid and alcohol ends are consumed "
            "one at a time, so the chain keeps its real -COOH and -OH "
            "terminal groups."
        ),
    },
    {
        "id": "C2",
        "family": "Coupling",
        "recipe": (
            "AA+BB; HA AA+BB; HA AB (amidine OK; amide/SO2/P=O blocked)"
        ),
        "reaction": _key(
            "([C&X3:1](=O)[O&X2&H1,Cl,Br]."
            "[C&X3:2](=O)[O&X2&H1,Cl,Br]).([N&X3;H2,H1;",
            "step_growth",
        ),
        "monomers": [
            "O=C(O)CCCCC(=O)O",
            "NCCCCCCN",
        ],
        "length": 2,
        "polymer": "Nylon-6,6",
        "status": "commercial",
        "note": (
            "Commodity polyamide from adipic acid and "
            "hexamethylenediamine; the amide coupling leaves acid and "
            "amine ends."
        ),
    },
    {
        "id": "C3",
        "family": "Coupling",
        "recipe": "AA+BB (NCO + guarded OH/SH)",
        "reaction": _key(
            "([N&X2:1]=[C&X2]=[O&X1,S&X1:2]", "step_growth"
        ),
        "monomers": ["O=C=NCCCCCCN=C=O", "OCCCCO"],
        "length": 2,
        "polymer": "Polyurethane (HDI + 1,4-butanediol)",
        "status": "commercial",
        "note": (
            "Polyurethane elastomer/coating chemistry; the urethane "
            "coupling consumes one isocyanate and one alcohol, and the "
            "unreacted groups remain the chain ends."
        ),
    },
    {
        "id": "C4",
        "family": "Coupling",
        "recipe": (
            "AA+BB (cyclic anhydride + diamine; 5-/6-membered only)"
        ),
        "reaction": _key(
            "([C&X3,c;R:1](=[O&X1])", "step_growth"
        ),
        "monomers": [
            "O=C1OC(=O)c2cc3c(=O)oc(=O)c3cc12",
            "Nc1ccc(Oc2ccc(N)cc2)cc1",
        ],
        "length": 2,
        "polymer": "Kapton-type polyimide (PMDA + ODA)",
        "status": "commercial",
        "note": (
            "Kapton-type polyimide from pyromellitic dianhydride and "
            "4,4'-oxydianiline; the cyclic dianhydride and diamine form "
            "the net cyclodehydration product."
        ),
    },
    {
        "id": "C5",
        "family": "Coupling",
        "recipe": "AA+BB (epoxide + NCO; 4 regio tiers)",
        "reaction": _key(
            "([C&X4;H2,H1,H0;R:1]1[O&X2&R:2]", "step_growth"
        ),
        "monomers": [
            "CC(C)(c1ccc(OCC2CO2)cc1)c1ccc(OCC2CO2)cc1",
            "O=C=Nc1ccc(Cc2ccc(N=C=O)cc2)cc1",
        ],
        "length": 2,
        "polymer": "Polyoxazolidone (BADGE + MDI)",
        "status": "known/specialty",
        "note": (
            "Specialty high-performance polymer from bisphenol A "
            "diglycidyl ether and MDI; the isocyanate N attacks the less "
            "hindered epoxide carbon to form the oxazolidinone ring."
        ),
    },
    {
        "id": "C6",
        "family": "Coupling",
        "recipe": "AB (diol + phosgene)",
        "reaction": _key(
            "([O,S;X2;H1;!$([O,S]C=*):1]."
            "[O,S;X2;H1;!$([O,S]C=*):2]).",
            "step_growth",
        ),
        "monomers": ["CC(C)(c1ccc(O)cc1)c1ccc(O)cc1"],
        "length": 2,
        "polymer": "Bisphenol A polycarbonate (BPA-PC)",
        "status": "commercial",
        "note": (
            "Lexan-type polycarbonate made by the phosgene route; the "
            "carbonate carbonyl and oxygen come from the phosgene reagent, "
            "so no atom is invented."
        ),
    },
    {
        "id": "C7",
        "family": "Coupling",
        "recipe": (
            "PEEK; PES; PAEN (ketone/sulfone/nitrile; same-ring or "
            "phenyl-spaced; F before Cl/Br/I)"
        ),
        "reaction": _key(
            "[c:1]1[c:2][c:3](F)[c:4]", "step_growth"
        ),
        "monomers": [
            "O=C(c1ccc(F)cc1)c1ccc(F)cc1",
            "Oc1ccc(O)cc1",
        ],
        "length": 2,
        "polymer": "Poly(ether ether ketone) (PEEK)",
        "status": "commercial",
        "note": (
            "Victrex-type PEEK made by nucleophilic aromatic substitution; "
            "the para-ketone activates the aryl fluoride, and the second "
            "Ar-F remains an exact terminal group."
        ),
    },
    {
        "id": "C8",
        "family": "Coupling",
        "recipe": "AB (ortho-blocked phenol)",
        "reaction": _key("c1([O&H1:1])", "step_growth"),
        "monomers": ["Cc1cccc(C)c1O"],
        "length": 2,
        "polymer": "Poly(2,6-dimethyl-1,4-phenylene oxide) (PPO)",
        "status": "commercial",
        "note": (
            "Noryl-type PPO made by oxidative phenol coupling; the "
            "phenolic O-H couples to the para C-H of the next ring, "
            "leaving one free phenolic end."
        ),
    },
    {
        "id": "R1",
        "family": "Ring-opening",
        "recipe": "lactone ROP; water init.; carbamic-acid end allowed",
        "reaction": _key(
            "[C&X3&R:1](=[O&X1])[O&X2&R:2]", "ring_opening"
        ),
        "monomers": ["O=C1CCCCCO1"],
        "length": 2,
        "polymer": "Polycaprolactone (PCL)",
        "status": "commercial",
        "note": (
            "Biodegradable polyester; water opens caprolactone to "
            "6-hydroxyhexanoic acid, and later cycles insert the opened "
            "hydroxyacid with the ROP-specific ester coupling."
        ),
    },
    {
        "id": "R2",
        "family": "Ring-opening",
        "recipe": (
            "lactam ROP; water init.; amidine OK, amide/SO2/P=O blocked"
        ),
        "reaction": _key(
            "[C&X3&R:1](=[O&X1])[N&X3&R:2]", "ring_opening"
        ),
        "monomers": ["O=C1CCCCCN1"],
        "length": 2,
        "polymer": "Nylon-6",
        "status": "commercial",
        "note": (
            "Commodity polyamide; water opens caprolactam to "
            "6-aminohexanoic acid, and the ROP-specific amide coupling "
            "allows the real aminoacid ends."
        ),
    },
    {
        "id": "R3",
        "family": "Ring-opening",
        "recipe": (
            "ROCOP (5-/6-membered anhydride, C/N/O/S backbone + "
            "epoxide; water init.)"
        ),
        "reaction": _key(
            "[C,c;R:1][C&X3,c;R](=[O&X1])[O&X2,o;R]",
            "ring_opening",
        ),
        "monomers": ["O=C1CCC(=O)O1", "C1CO1"],
        "length": 2,
        "polymer": "Poly(ethylene succinate) (PES)",
        "status": "known",
        "note": (
            "Biodegradable ROCOP polyester; the balanced opening consumes "
            "one water molecule, the anhydride O becomes the ester O, and "
            "the epoxide O becomes the alcohol end."
        ),
    },
    {
        "id": "R4",
        "family": "Ring-opening",
        "recipe": "epoxide ROP (direct; 4 regio tiers)",
        "reaction": _key(
            "[C&X4;H2,H1,H0;R:1]1[O&R]", "step_growth"
        ),
        "monomers": ["C1CO1"],
        "length": 3,
        "polymer": "Poly(ethylene oxide) (PEO)",
        "status": "commercial",
        "note": (
            "Commodity polyether; water opens the epoxide directly at the "
            "less hindered carbon, supplying one -OH end and the epoxide O "
            "as the other. Water is the initiator, not an epoxide residue, "
            "so length 3 shows two epoxide units (diethylene glycol)."
        ),
    },
    {
        "id": "R5",
        "family": "Ring-opening",
        "recipe": "ROMP; cuts C=C; joins via C=C",
        "reaction": _key(
            "[C&X3&H1&R:1]=[C&X3&H1&R:2][C&X4&R:3]", "romp"
        ),
        "monomers": ["C1CC2C=CC1C2"],
        "length": 2,
        "polymer": "Poly(norbornene) (PNB, ROMP)",
        "status": "commercial",
        "note": (
            "Norsorex-type ROMP polymer; the metathesis event cuts the ring "
            "C=C, so the new double bonds sit at the junctions between "
            "repeat units."
        ),
    },
]
# If a mechanism has no identified practical polymer, set
# ``"status": "needs candidate"``, leave ``monomers`` empty, and explain the
# gap in ``note``.  The page then writes a placeholder for you to fill in.


# Executable checks for the monomer-scope fixes.  Each entry names a
# mechanism, the monomers to run, and whether a product or a rejection is
# expected.  The generator runs them and fails if the observed outcome does
# not match ``expect``.
_GENERALIZATION_CHECKS = [
    {
        "mechanism": "A1",
        "monomers": ["C=CC"],
        "length": 2,
        "expect": "product",
        "note": "terminal vinyl propagates head-to-tail",
    },
    {
        "mechanism": "A1",
        "monomers": ["C=C1CC(C)(c2cnc[nH]2)C1"],
        "length": 1,
        "expect": "product",
        "note": "exocyclic C=C is acyclic and propagates",
    },
    {
        "mechanism": "A1",
        "monomers": ["C1CC2C=CC1C2"],
        "length": 1,
        "expect": "reject",
        "note": "endocyclic C=C does not leak into A1",
    },
    {
        "mechanism": "A1",
        "monomers": ["C=C=C"],
        "length": 1,
        "expect": "reject",
        "note": "allene/ketene-like cumulene does not match",
    },
    {
        "mechanism": "C1",
        "monomers": ["O=C(O)c1ccc(C(=O)O)cc1", "OCCO"],
        "length": 2,
        "expect": "product",
        "note": "acid + diol esterification",
    },
    {
        "mechanism": "C1",
        "monomers": ["COC(=O)c1ccc(C(=O)OC)cc1", "OCCO"],
        "length": 2,
        "expect": "product",
        "note": "DMT methyl-ester transesterification",
    },
    {
        "mechanism": "C7",
        "monomers": ["O=C(c1ccc(F)cc1)c1ccc(F)cc1", "Oc1ccc(O)cc1"],
        "length": 2,
        "expect": "product",
        "note": "PEEK para-ketone activation",
    },
    {
        "mechanism": "C7",
        "monomers": ["N#Cc1ccc(F)cc1", "Oc1ccc(O)cc1"],
        "length": 2,
        "expect": "product",
        "note": "para-nitrile activation",
    },
    {
        "mechanism": "C7",
        "monomers": ["Fc1cccc(F)c1C#N", "Oc1ccc(O)cc1"],
        "length": 2,
        "expect": "product",
        "note": "PAEN 2,6-difluorobenzonitrile (ortho nitrile)",
    },
    {
        "mechanism": "C7",
        "monomers": [
            "O=C(c1ccc(-c2ccc(F)cc2)cc1)c1ccc(-c2ccc(F)cc2)cc1",
            "Oc1ccc(O)cc1",
        ],
        "length": 2,
        "expect": "product",
        "note": "biphenyl-spaced ketone activation",
    },
    {
        "mechanism": "C7",
        "monomers": [
            "O=C(c1ccc2cc(F)ccc2c1)c1ccc2cc(F)ccc2c1",
            "Oc1ccc(O)cc1",
        ],
        "length": 2,
        "expect": "product",
        "note": "naphthyl-spaced ketone activation",
    },
    {
        "mechanism": "C8",
        "monomers": ["Cc1cccc(C)c1O"],
        "length": 3,
        "expect": "product",
        "note": "2,6-disubstituted phenol couples linearly",
    },
    {
        "mechanism": "C8",
        "monomers": ["Oc1ccccc1"],
        "length": 3,
        "expect": "reject",
        "note": "unblocked ortho positions are rejected",
    },
    {
        "mechanism": "R3",
        "monomers": ["O=C1CCC(=O)O1", "C1CO1"],
        "length": 2,
        "expect": "product",
        "note": "succinic anhydride ROCOP",
    },
    {
        "mechanism": "R3",
        "monomers": ["O=C1COC(=O)O1", "C1CO1"],
        "length": 2,
        "expect": "product",
        "note": "diglycolic (heteroatom) anhydride ROCOP",
    },
]


def _mechanism_key(mechanism_id: str) -> str:
    """Return the registry key used by the practical example for *mechanism_id*."""
    return next(
        spec["reaction"] for spec in _MECHANISMS
        if spec["id"] == mechanism_id
    )


def _run_generalization_checks() -> list[dict]:
    """Run the monomer-scope checks and return their observed outcomes."""
    results = []
    for check in _GENERALIZATION_CHECKS:
        reaction = _mechanism_key(check["mechanism"])
        config = rxn_map[reaction]
        driver = _build(config["type"], reaction)
        try:
            products = driver.extend(check["monomers"], check["length"])
            if products:
                outcome = "product"
                result = products[0]
            else:
                outcome = "reject"
                result = "no product"
        except Exception as exc:  # noqa: BLE001 - a rejection is an outcome
            outcome = "reject"
            result = f"rejected ({type(exc).__name__})"
        results.append(
            {
                **check,
                "outcome": outcome,
                "result": result,
                "status": "pass" if outcome == check["expect"] else "FAIL",
            }
        )
    failures = [row for row in results if row["status"] != "pass"]
    if failures:
        details = "; ".join(
            f"{row['mechanism']} {row['note']}: {row['outcome']} "
            f"({row['result']})"
            for row in failures
        )
        raise ValueError(f"Generalization checks failed: {details}")
    return results


# Additional conventional traces.  They are appended after PP and PET in the
# default page; simple/shortest modes keep the original two examples.
_CONVENTIONAL_EXTRA = [
    {
        "kind": "step_growth",
        "reaction": _EPOXIDE_ROP,
        "monomers": ["C1CO1"],
        "polymer_name": "Poly(ethylene oxide) (PEO)",
        "section": "Epoxide ring-opening trace",
        "prefix": "peo",
        "length": 4,
        "challenge": (
            "The epoxide has no pre-existing terminal group. The configured "
            "recipe supplies water as the initiator; the first coupling "
            "opens the epoxide to ethylene glycol, and every later cycle "
            "inserts another opened epoxide between the growing -OH end and "
            "the new -OH end. A wildcard join would lose the water-derived "
            "end group and the directional C-O repeat."
        ),
        "mechanism_note": (
            "This is an epoxide ring-opening entry driven by the "
            "`step_growth` driver: water is an explicit reagent, and the "
            "first coupling is the ring-opening/initiation step."
        ),
        "ring_opening": True,
        "initiator_label": "H2O",
        "pool_note": (
            "Cycle 1 is the epoxide + water pool. Water opens the epoxide "
            "to ethylene glycol, so residue 2 has the same heavy-atom count "
            "as the pool and each later residue adds one opened epoxide."
        ),
        "linkage_note": (
            "The trace uses the `O` branch of the `[O,S]` coupling, so the "
            "linkage is a normal ether."
        ),
        "residue_note": (
            "Each residue is a complete molecule; the -OH ends come from "
            "water and from the opened epoxide oxygen."
        ),
    },
    {
        "kind": "ring_opening",
        "reaction": _LACTONE_ROP,
        "monomers": ["O=C1CCCCCO1"],
        "polymer_name": "Polycaprolactone (PCL)",
        "section": "Lactone ring-opening trace",
        "prefix": "pcl",
        "length": 4,
        "challenge": (
            "The lactone has no free -OH or -COOH before ring opening. "
            "Water is the initiator: opening the first lactone adds water "
            "across the acyl-O bond, so the first residue is "
            "6-hydroxyhexanoic acid and has one more oxygen than the "
            "cyclic lactone. The opened hydroxyacid is then propagated "
            "with a ROP-specific ester coupling. A wildcard join would "
            "lose the initiator-derived acid end."
        ),
        "mechanism_note": (
            "Cycle 1 is the opened hydroxyacid. The water initiator O is "
            "the acid O-H, so residue 1 has one more oxygen than the "
            "cyclic lactone. Each later residue adds only the monomer's "
            "own atoms. This is a topological open-then-condense model: "
            "it reproduces the terminated chain directly and does not "
            "track an active-end dummy."
        ),
        "initiator_label": "H2O",
    },
    {
        "kind": "ring_opening",
        "reaction": _LACTAM_ROP,
        "monomers": ["O=C1CCCCCN1"],
        "polymer_name": "Nylon-6",
        "section": "Lactam ring-opening trace",
        "prefix": "nylon6",
        "length": 4,
        "challenge": (
            "The lactam has no free amine or acid before ring opening. "
            "Water is the initiator: opening the first lactam adds water "
            "across the acyl-N bond, so the first residue is "
            "6-aminohexanoic acid and has one more oxygen than the cyclic "
            "lactam. The opened aminoacid is then propagated with a "
            "ROP-specific amide coupling. The generic amide guard "
            "`!$(NC=*)` is too broad for this mechanism."
        ),
        "mechanism_note": (
            "Cycle 1 is the opened aminoacid. The water initiator O is "
            "the acid O-H, so residue 1 has one more oxygen than the "
            "cyclic lactam. Each later residue adds only the monomer's "
            "own atoms. This is a topological open-then-condense model: "
            "it reproduces the terminated chain directly and does not "
            "track an active-end dummy."
        ),
        "initiator_label": "H2O",
    },
    {
        "kind": "ring_opening",
        "reaction": _ROCOP,
        "monomers": ["O=C1CCC(=O)O1", "C1CO1"],
        "polymer_name": "Poly(ethylene succinate) (PES)",
        "section": "Anhydride + epoxide ROCOP trace",
        "prefix": "pes",
        "length": 4,
        "challenge": (
            "Neither monomer has a free -OH or -COOH before opening. "
            "Water is the initiator: the balanced opening consumes one "
            "water molecule. The anhydride ring O becomes the ester O, "
            "the water O becomes the acid O-H, and the epoxide O becomes "
            "the terminal alcohol. Residue 1 therefore has one more "
            "oxygen than the anhydride + epoxide pair. One cycle consumes "
            "one anhydride + one epoxide pair."
        ),
        "mechanism_note": (
            "Cycle 1 is the opened anhydride + epoxide pair. The water "
            "initiator O is the acid O-H, so residue 1 has one more "
            "oxygen than the two monomers together. Each later cycle "
            "inserts another pair and adds only the monomers' own atoms. "
            "This is a topological open-then-condense model: it reproduces "
            "the terminated chain directly and does not track an "
            "alternating active-end dummy."
        ),
        "initiator_label": "H2O",
    },
    {
        "kind": "step_growth",
        "reaction": _PHENOL_OXIDATIVE,
        "monomers": ["Cc1cccc(C)c1O"],
        "polymer_name": "Poly(2,6-dimethyl-1,4-phenylene oxide) (PPO)",
        "section": "Oxidative phenol coupling trace",
        "prefix": "ppo",
        "length": 4,
        "challenge": (
            "This is a C-O coupling, not a ring-opening reaction. The "
            "phenolic O-H of one ring couples to the para C-H of another "
            "ring; the second ring's O-H is carried over as the real "
            "terminal group. A wildcard join would discard that phenol end "
            "and could couple at the wrong ring position."
        ),
        "mechanism_note": (
            "Each cycle couples one phenolic O-H to the para C-H of the "
            "next ring, leaving one free phenolic -OH at the chain end."
        ),
        "linkage_note": (
            "The linkage is a diaryl ether (-Ar-O-Ar-), not an ester."
        ),
    },
    {
        "kind": "romp",
        "reaction": _ROMP,
        "monomers": ["C1CC2C=CC1C2"],
        "polymer_name": "Poly(norbornene) (PNB)",
        "section": "ROMP trace",
        "prefix": "romp",
        "length": 2,
        "challenge": (
            "ROMP is chain-growth metathesis, not step growth. The map key "
            "cuts the metathesis C=C double bond (not the allylic single "
            "bond), so the double bonds sit at the junctions between repeat "
            "units. The ring path is variable-length, so the driver performs "
            "the cut programmatically instead of using a generic "
            "two-reactant SMARTS. A wildcard join would not encode the "
            "junction double bonds or the initiator/terminator caps."
        ),
        "mechanism_note": (
            "The `*` marks the growing carbene end; termination replaces it "
            "with the configured terminator fragment."
        ),
    },
]

_SIMPLE = {
    "addition": {
        "reaction": _VINYL,
        "monomers": ["C=CC"],
        "polymer_name": "Polypropylene (minimal propene example)",
        "challenge": "",
    },
    "step_growth": {
        "reaction": _GLYCOLIC_AB,
        "monomers": ["OCC(=O)O"],
        "polymer_name": "Poly(glycolic acid) (PGA, minimal hydroxyacid example)",
        "challenge": "",
    },
}


def _monomers(row) -> list[str]:
    return [
        value
        for value in (row["mon1"], row["mon2"])
        if isinstance(value, str) and value.strip()
    ]


def _heavy_atoms(smiles: str) -> int:
    mol = Chem.MolFromSmiles(smiles)
    return mol.GetNumHeavyAtoms() if mol is not None else 10**6


def _non_dummy_heavy_atoms(mol) -> int:
    return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() > 0)


def _rings(mol) -> int:
    return mol.GetRingInfo().NumRings()


def _target_length(kind: str, mode: str) -> int:
    """Number of residues/propagation steps to show.

    Chain growth and ROMP always show two propagation steps.  The
    conventional PET example uses four residues (three couplings) so both
    monomers appear twice; the other step-growth/ring-opening modes keep three
    residues.
    """
    if kind in {"addition", "romp"}:
        return _PROPAGATION_STEPS
    if mode == "conventional":
        return 4
    return _PROPAGATION_STEPS + 1


def _build(kind: str, reaction: str, max_products: int = 1):
    config = rxn_map[reaction]
    if kind == "addition":
        return ChainGrowthReaction(
            **config["params"], max_products=max_products
        )
    if kind == "romp":
        return RingOpeningMetathesisReaction(
            **reaction_config(reaction), max_products=max_products
        )
    if kind == "ring_opening":
        return RingOpeningReaction(
            **reaction_config(reaction), max_products=max_products
        )
    return StepGrowthReaction(
        **reaction_config(reaction), max_products=max_products
    )


def _candidates(kind: str):
    frame = pd.read_csv(_CSV)
    candidates = []
    for _, row in frame.iterrows():
        config = rxn_map.get(row["reaction_smarts"])
        if config is None or config.get("type") != kind:
            continue
        monomers = _monomers(row)
        if not monomers:
            continue
        heavy = sum(_heavy_atoms(monomer) for monomer in monomers)
        length = sum(len(monomer) for monomer in monomers)
        candidates.append(
            (
                heavy,
                length,
                row["reaction_smarts"],
                row["mon1"],
                row["mon2"],
                monomers,
                config,
            )
        )
    candidates.sort(
        key=lambda item: (
            item[0],
            item[1],
            item[2],
            str(item[3]),
            str(item[4]),
        )
    )
    return candidates


def _selection(
    kind: str,
    reaction: str,
    monomers: list[str],
    description: str,
    polymer_name: str,
    challenge: str,
    length_arg: int,
    sequence_note: str = "",
    **extra,
) -> dict:
    config = rxn_map[reaction]
    if not extra.get("paired"):
        driver = _build(kind, reaction)
        products = driver.extend(monomers, length_arg)
        if not products or any("*" in smiles for smiles in products):
            raise RuntimeError(
                f"{description} did not produce a valid product"
            )
    selection = {
        "kind": kind,
        "heavy": sum(_heavy_atoms(monomer) for monomer in monomers),
        "length": sum(len(monomer) for monomer in monomers),
        "reaction": reaction,
        "mon1": monomers[0],
        "mon2": monomers[1] if len(monomers) > 1 else None,
        "monomers": monomers,
        "config": config,
        "description": description,
        "polymer_name": polymer_name,
        "challenge": challenge,
        "sequence_note": sequence_note,
        "length_arg": length_arg,
    }
    selection.update(extra)
    return selection


def _spec_selection(spec: dict) -> dict:
    """Build a selection from one conventional/extra trace spec."""
    return _selection(
        spec["kind"],
        spec["reaction"],
        spec["monomers"],
        f"Conventional {spec['polymer_name']} example.",
        spec["polymer_name"],
        spec.get("challenge", ""),
        spec.get("length", _PROPAGATION_STEPS),
        spec.get("sequence_note", ""),
        section=spec.get("section", ""),
        prefix=spec.get("prefix", spec["kind"]),
        mechanism_note=spec.get("mechanism_note", ""),
        pool_note=spec.get("pool_note", ""),
        linkage_note=spec.get("linkage_note", ""),
        ring_opening=spec.get("ring_opening", False),
        residue_note=spec.get("residue_note", ""),
        initiator_label=spec.get("initiator_label", ""),
    )


def _mechanism_selection(spec: dict) -> dict:
    """Build a trace selection for one entry of ``_MECHANISMS``."""
    config = rxn_map[spec["reaction"]]
    if config.get("mechanism") != spec["id"]:
        raise ValueError(
            f"Mechanism mismatch for {spec['id']}: registry has "
            f"{config.get('mechanism')!r} for {spec['reaction']!r}"
        )
    kind = config["type"]
    selection = _selection(
        kind,
        spec["reaction"],
        spec["monomers"],
        (
            f"Mechanism {spec['id']} example: {spec['recipe']}."
        ),
        spec["polymer"],
        "",
        spec["length"],
        section=f"Mechanism {spec['id']}",
        prefix=spec.get("prefix", spec["id"]),
        mechanism_id=spec["id"],
        mechanism_family=spec["family"],
        recipe=spec["recipe"],
        mechanism_note=spec["note"],
        status=spec.get("status", "commercial"),
        paired=spec.get("paired", False),
        ring_opening=kind == "ring_opening"
        or bool(config["params"].get("ring_open")),
    )
    return selection


def _select(kind: str, mode: str) -> dict:
    length_arg = _target_length(kind, mode)
    if mode == "conventional":
        reference = _CONVENTIONAL[kind]
        return _selection(
            kind,
            reference["reaction"],
            reference["monomers"],
            f"Conventional {reference['polymer_name']} example.",
            reference["polymer_name"],
            reference["challenge"],
            length_arg,
            reference.get("sequence_note", ""),
            section=(
                "Chain-growth trace"
                if kind == "addition"
                else "Step-growth trace"
            ),
            prefix="chain" if kind == "addition" else "step",
        )
    if mode == "simple":
        reference = _SIMPLE[kind]
        return _selection(
            kind,
            reference["reaction"],
            reference["monomers"],
            f"Minimal {reference['polymer_name']} example.",
            reference["polymer_name"],
            reference["challenge"],
            length_arg,
            section=(
                "Chain-growth trace"
                if kind == "addition"
                else "Step-growth trace"
            ),
            prefix="chain" if kind == "addition" else "step",
        )

    for candidate in _candidates(kind):
        heavy, length, reaction, mon1, mon2, monomers, config = candidate
        try:
            driver = _build(kind, reaction)
            products = driver.extend(monomers, length_arg)
        except Exception:
            continue
        if products and all("*" not in smiles for smiles in products):
            return {
                "kind": kind,
                "heavy": heavy,
                "length": length,
                "reaction": reaction,
                "mon1": mon1,
                "mon2": mon2,
                "monomers": monomers,
                "config": config,
                "description": (
                    "Shortest surviving candidate in the extended test set."
                ),
                "polymer_name": (
                    f"{config['polymer_class']} (shortest candidate)"
                ),
                "challenge": "",
                "sequence_note": "",
                "length_arg": length_arg,
                "section": (
                    "Chain-growth trace"
                    if kind in {"addition", "romp"}
                    else "Step-growth trace"
                ),
                "prefix": (
                    "chain" if kind in {"addition", "romp"} else "step"
                ),
            }
    raise RuntimeError(f"No surviving {kind} candidate found")


def _first_reaction_product(rxn, reactants):
    """Return the first single connected, sanitised product of *rxn*."""
    for result in rxn.RunReactants(reactants):
        if not result:
            continue
        molecule = result[0]
        if molecule is None or molecule.GetNumAtoms() == 0:
            continue
        try:
            Chem.SanitizeMol(molecule)
        except Exception:
            continue
        if len(Chem.GetMolFrags(molecule)) != 1:
            continue
        return molecule
    return None


def _paired_chain_trace(selection: dict) -> dict:
    """Trace a two-comonomer addition step with one unit of each.

    The generic chain-growth driver enumerates the monomer pool with
    replacement, so a frontier cap of one can pick the same monomer in both
    slots.  A paired mechanism must show the actual comonomer pair, so this
    helper fills the two propagation slots with ``mon1`` and ``mon2``
    explicitly.
    """
    driver = _build("addition", selection["reaction"])
    mon1, mon2 = [
        _molecule(smiles) for smiles in selection["monomers"]
    ]
    chain = _first_reaction_product(driver.init_rxn, (mon1, mon2))
    if chain is None:
        raise RuntimeError(
            f"Mechanism {selection['mechanism_id']} initiation failed"
        )
    intermediates = [[Chem.MolToSmiles(chain)]]
    for _ in range(selection["length_arg"]):
        chain = _first_reaction_product(
            driver.prop_rxn, (chain, mon1, mon2)
        )
        if chain is None:
            raise RuntimeError(
                f"Mechanism {selection['mechanism_id']} propagation failed"
            )
        intermediates.append([Chem.MolToSmiles(chain)])
    final = _first_reaction_product(driver.term_rxn, (chain,))
    if final is None:
        raise RuntimeError(
            f"Mechanism {selection['mechanism_id']} termination failed"
        )
    return {
        "intermediates": intermediates,
        "final": [Chem.MolToSmiles(final)],
    }


def _trace(kind: str, selection: dict) -> dict:
    if selection.get("paired"):
        return _paired_chain_trace(selection)
    driver = _build(kind, selection["reaction"])
    result = driver.extend(
        selection["monomers"],
        selection["length_arg"],
        return_intermediate=True,
    )
    if not isinstance(result, dict):
        raise TypeError("trace extraction requires return_intermediate=True")
    return result


def _draw(mol, path: Path, legend: str) -> None:
    atoms = mol.GetNumAtoms()
    width = min(1200, max(520, 34 * atoms))
    height = 360
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    options = drawer.drawOptions()
    options.padding = 0.08
    options.bondLineWidth = 2.2
    options.addStereoAnnotation = True
    rdMolDraw2D.PrepareAndDrawMolecule(
        drawer, mol, legend=legend.replace("\u2013", "-").replace(
            "\u2014", "-"
        )
    )
    drawer.FinishDrawing()
    path.write_text(drawer.GetDrawingText())


def _molecule(smiles: str):
    return Chem.MolFromSmiles(smiles)


def _draw_result(
    selection: dict, result: dict, out_dir: Path | None = None
) -> None:
    """Draw every intermediate/final structure for one trace."""
    target = out_dir or _TRACES
    kind = selection["kind"]
    prefix = selection.get("prefix") or (
        "chain" if kind in {"addition", "romp"} else "step"
    )
    if kind in {"addition", "romp"}:
        for index, smiles in enumerate(result["intermediates"]):
            label = "init" if index == 0 else f"propagation {index}"
            _draw(
                _molecule(smiles[0]),
                target / f"{prefix}_{index:02d}.svg",
                label,
            )
        for smiles in result["final"]:
            _draw(
                _molecule(smiles),
                target / f"{prefix}_final.svg",
                "terminated",
            )
        return

    for index, smiles in enumerate(result["intermediates"]):
        for offset, smile in enumerate(smiles):
            suffix = (
                "" if len(smiles) == 1 else chr(ord("a") + offset)
            )
            label = (
                f"residue {index + 1}"
                if len(smiles) == 1
                else f"monomer {offset + 1}"
            )
            _draw(
                _molecule(smile),
                target / f"{prefix}_{index + 1:02d}{suffix}.svg",
                label,
            )


def _trace_entries(selection: dict, result: dict) -> list[tuple[str, str]]:
    """Return ``(label, SMILES)`` pairs for every visible trace structure."""
    kind = selection["kind"]
    entries: list[tuple[str, str]] = []
    if kind in {"addition", "romp"}:
        for index, smiles in enumerate(result["intermediates"]):
            if not smiles:
                continue
            label = "init" if index == 0 else f"propagation {index}"
            entries.append((label, smiles[0]))
        for smiles in result["final"]:
            entries.append(("terminated", smiles))
        return entries

    for index, smiles in enumerate(result["intermediates"]):
        if len(smiles) == 1:
            entries.append((f"residue {index + 1}", smiles[0]))
        else:
            for offset, smile in enumerate(smiles):
                entries.append((f"monomer {offset + 1}", smile))
    return entries


def _draw_grid(
    entries: list[tuple[str, str]],
    path: Path,
    mols_per_row: int = 2,
) -> None:
    """Draw a labelled SVG panel of trace structures."""
    molecules = []
    legends = []
    for label, smiles in entries:
        molecule = _molecule(smiles)
        if molecule is None:
            continue
        molecules.append(molecule)
        legends.append(label)
    if not molecules:
        return
    svg = MolsToGridImage(
        molecules,
        molsPerRow=mols_per_row,
        subImgSize=(390, 290),
        legends=legends,
        useSVG=True,
    )
    if isinstance(svg, bytes):
        svg = svg.decode("utf-8")
    path.write_text(svg)


def _row(step: str, smiles, include_dummy: bool = False) -> str:
    smiles_list = [smiles] if isinstance(smiles, str) else list(smiles)
    heavy = 0
    rings = 0
    for smile in smiles_list:
        mol = _molecule(smile)
        if mol is None:
            heavy = "?"
            rings = "?"
            break
        heavy += (
            _non_dummy_heavy_atoms(mol)
            if include_dummy
            else mol.GetNumHeavyAtoms()
        )
        rings += _rings(mol)
    smiles_text = " + ".join(f"`{smile}`" for smile in smiles_list)
    return f"| {step} | {heavy} | {rings} | {smiles_text} |"


def _write_challenge(selection: dict, lines: list[str]) -> None:
    if selection.get("challenge"):
        lines.append(
            "**Why wildcard joining is inaccurate.** "
            f"{selection['challenge']}"
        )
        lines.append("")


def _write_chain_section(
    selection: dict, result: dict, lines: list[str]
) -> None:
    reaction = selection["reaction"]
    config = selection["config"]
    params = config["params"]
    monomers = " + ".join(
        f"`{value}`" for value in selection["monomers"]
    )
    lines.extend(
        [
            f"## Chain-growth trace — {selection['polymer_name']}",
            "",
            selection["description"],
            "",
            "| field | value |",
            "| --- | --- |",
            f"| polymer | {selection['polymer_name']} |",
            f"| polymer class | {config['polymer_class']} |",
            "| driver | `ChainGrowthReaction` |",
            f"| monomer(s) | {monomers} |",
            f"| heavy atoms | {selection['heavy']} |",
            f"| reaction SMARTS | `{reaction}` |",
            f"| `init_smarts` | `{params['init_smarts']}` |",
            f"| `prop_smarts` | `{params['prop_smarts']}` |",
            f"| `term_smarts` | `{params['term_smarts']}` |",
            "",
            (
                "In chain-growth intermediates the `*` is the dummy atom "
                "at the active growing end. It is removed by termination in "
                "the final row. Heavy-atom counts below exclude the dummy."
            ),
            "",
        ]
    )
    _write_challenge(selection, lines)
    lines.extend(
        [
            "| step | heavy atoms | rings | SMILES |",
            "| --- | --- | --- | --- |",
        ]
    )
    for index, smiles in enumerate(result["intermediates"]):
        label = (
            "init (n=0)"
            if index == 0
            else f"propagation {index} (n={index})"
        )
        lines.append(_row(label, smiles[0], include_dummy=True))
    for smiles in result["final"]:
        lines.append(_row("terminated", smiles))
    lines.append("")

    intermediate_heavies = [
        _non_dummy_heavy_atoms(_molecule(smiles[0]))
        for smiles in result["intermediates"]
    ]
    final_heavy = (
        _molecule(result["final"][0]).GetNumHeavyAtoms()
        if result["final"]
        else None
    )
    if len(intermediate_heavies) > 1 and final_heavy is not None:
        lines.append(
            "Growth accounting: initiation gives "
            f"{intermediate_heavies[0]} heavy atoms; each propagation adds "
            f"{intermediate_heavies[1] - intermediate_heavies[0]} (one "
            "monomer unit); termination adds "
            f"{final_heavy - intermediate_heavies[-1]} (the hydrogen cap)."
        )
        lines.append("")
    lines.append("### Structures")
    lines.append("")
    for index, smiles in enumerate(result["intermediates"]):
        structure = (
            "initiation" if index == 0 else f"propagation {index}"
        )
        lines.append(f"#### {structure}")
        lines.append("")
        lines.append(f"![{structure}](traces/chain_{index:02d}.svg)")
        lines.append("")
    lines.append("#### termination")
    lines.append("")
    lines.append("![termination](traces/chain_final.svg)")
    lines.append("")


def _write_romp_section(
    selection: dict, result: dict, lines: list[str]
) -> None:
    reaction = selection["reaction"]
    config = selection["config"]
    params = config["params"]
    ring_bond = reaction_config(reaction).get(
        "ring_bond_smarts", reaction
    )
    ring_bond = str(ring_bond).split(">")[0]
    monomers = " + ".join(
        f"`{value}`" for value in selection["monomers"]
    )
    section = selection.get("section") or "ROMP trace"
    lines.extend(
        [
            f"## {section} — {selection['polymer_name']}",
            "",
            selection["description"],
            "",
            "| field | value |",
            "| --- | --- |",
            f"| polymer | {selection['polymer_name']} |",
            f"| polymer class | {config['polymer_class']} |",
            "| driver | `RingOpeningMetathesisReaction` |",
            f"| monomer(s) | {monomers} |",
            f"| heavy atoms | {selection['heavy']} |",
            f"| reaction SMARTS | `{reaction}` |",
            f"| ring-bond SMARTS | `{ring_bond}` |",
            f"| initiator | `{params['initiator']}` |",
            f"| terminator | `{params['terminator']}` |",
            "",
            selection.get("mechanism_note")
            or (
                "The `*` marks the growing carbene end; termination "
                "replaces it with the configured terminator fragment."
            ),
            "",
            (
                "In ROMP intermediates the `*` is the dummy atom at the "
                "active growing end. It is removed by termination in the "
                "final row. Heavy-atom counts below exclude the dummy."
            ),
            "",
        ]
    )
    _write_challenge(selection, lines)
    lines.extend(
        [
            "| step | heavy atoms | rings | SMILES |",
            "| --- | --- | --- | --- |",
        ]
    )
    for index, smiles in enumerate(result["intermediates"]):
        label = (
            "init (n=0)"
            if index == 0
            else f"propagation {index} (n={index})"
        )
        lines.append(_row(label, smiles[0], include_dummy=True))
    for smiles in result["final"]:
        lines.append(_row("terminated", smiles))
    lines.append("")

    intermediate_heavies = [
        _non_dummy_heavy_atoms(_molecule(smiles[0]))
        for smiles in result["intermediates"]
    ]
    final_heavy = (
        _molecule(result["final"][0]).GetNumHeavyAtoms()
        if result["final"]
        else None
    )
    if len(intermediate_heavies) > 1 and final_heavy is not None:
        lines.append(
            "Growth accounting: initiation gives "
            f"{intermediate_heavies[0]} heavy atoms; each propagation adds "
            f"{intermediate_heavies[1] - intermediate_heavies[0]} (one "
            "opened cycloalkene); termination adds "
            f"{final_heavy - intermediate_heavies[-1]} (the terminator)."
        )
        lines.append("")

    prefix = selection.get("prefix") or "romp"
    lines.append("### Structures")
    lines.append("")
    for index, smiles in enumerate(result["intermediates"]):
        structure = (
            "initiation" if index == 0 else f"propagation {index}"
        )
        lines.append(f"#### {structure}")
        lines.append("")
        lines.append(
            f"![{structure}](traces/{prefix}_{index:02d}.svg)"
        )
        lines.append("")
    lines.append("#### termination")
    lines.append("")
    lines.append(f"![termination](traces/{prefix}_final.svg)")
    lines.append("")


def _write_step_section(
    selection: dict, result: dict, lines: list[str]
) -> None:
    reaction = selection["reaction"]
    config = selection["config"]
    params = config["params"]
    monomers = " + ".join(
        f"`{value}`" for value in selection["monomers"]
    )
    is_ring_open = bool(selection.get("ring_opening")) or (
        config.get("type") in {"ring_open", "ring_opening"}
        or bool(params.get("ring_open"))
    )
    propagation_smarts = params.get("coupling_smarts") or params.get(
        "prop_smarts", []
    )
    couplings = ", ".join(
        f"`{value}`" for value in propagation_smarts
    )
    reagents = (
        ", ".join(f"`{value}`" for value in params.get("reagents", []))
        or "-"
    )
    opening_reagents = params.get("opening_reagents", [])
    opening_reagents_text = (
        ", ".join(f"`{value}`" for value in opening_reagents)
        if opening_reagents
        else None
    )
    if selection.get("initiator_label"):
        # The initiator is already shown in its own field; do not repeat the
        # same water/reagent under a second label.
        reagents = "-"
        opening_reagents_text = None
    contains_sulfur = any(
        "S" in monomer for monomer in selection["monomers"]
    )
    driver_name = (
        "RingOpeningReaction"
        if config.get("type") == "ring_opening"
        else "StepGrowthReaction"
    )
    section = selection.get("section") or (
        "Ring-opening trace" if is_ring_open else "Step-growth trace"
    )
    smarts_label = (
        "propagation SMARTS"
        if config.get("type") == "ring_opening"
        else "coupling SMARTS"
    )
    mechanism_note = selection.get("mechanism_note") or (
        (
            "This is a ring-opening entry, so cycle 1 is the opened "
            "monomer rather than the raw cyclic monomer. The "
            "`reaction SMARTS` field is the ring-opening key; the "
            "`propagation SMARTS` field is the ROP coupling compiled by "
            "the driver."
        )
        if is_ring_open
        else (
            "This is a condensation entry: the chain alternates the two "
            "monomers and keeps their real acid/alcohol end groups. The "
            "`reaction SMARTS` field is the registry key; the "
            "`coupling SMARTS` field is the recipe compiled by the driver."
        )
    )
    linkage_note = selection.get("linkage_note") or (
        (
            "This monomer contains sulfur; the `[O,S]` coupling can "
            "therefore produce a thioester linkage."
        )
        if contains_sulfur
        else (
            "This example uses the `O` branch of the `[O,S]` coupling, so "
            "the linkage is a normal ester."
        )
    )
    lines.extend(
        [
            f"## {section} — {selection['polymer_name']}",
            "",
            selection["description"],
            "",
            "| field | value |",
            "| --- | --- |",
            f"| polymer | {selection['polymer_name']} |",
            f"| polymer class | {config['polymer_class']} |",
            f"| driver | `{driver_name}` |",
            f"| monomer(s) | {monomers} |",
            f"| heavy atoms | {selection['heavy']} |",
            f"| reaction SMARTS | `{reaction}` |",
            f"| {smarts_label} | {couplings} |",
            f"| ring opening | {'yes' if is_ring_open else 'no'} |",
            f"| reagents | {reagents} |",
            *(
                [f"| initiator | {selection['initiator_label']} |"]
                if selection.get("initiator_label")
                else []
            ),
            *(
                [f"| ring-opening reagent(s) | {opening_reagents_text} |"]
                if opening_reagents_text
                else []
            ),
            "",
            mechanism_note,
            "",
            selection.get("residue_note")
            or (
                "Each residue is a complete molecule with the monomer's "
                "real terminal groups; there is no wildcard in the "
                "intermediates."
            ),
            "",
            linkage_note,
            "",
        ]
    )
    _write_challenge(selection, lines)
    if selection.get("sequence_note"):
        lines.append(selection["sequence_note"])
        lines.append("")
    lines.extend(
        [
            "| step | heavy atoms | rings | SMILES |",
            "| --- | --- | --- | --- |",
        ]
    )
    for index, smiles in enumerate(result["intermediates"]):
        label = (
            "monomers"
            if index == 0 and len(smiles) > 1
            else f"residue {index + 1}"
        )
        lines.append(_row(label, smiles))
    lines.append("")

    first = result["intermediates"][0]
    if len(first) == 1 and len(result["intermediates"]) > 1:
        residue_heavies = [
            _molecule(smiles[0]).GetNumHeavyAtoms()
            for smiles in result["intermediates"]
        ]
        first_label = "opened monomer" if is_ring_open else "first residue"
        lines.append(
            f"Growth accounting: the {first_label} has "
            f"{residue_heavies[0]} heavy atoms; each additional residue adds "
            f"{residue_heavies[1] - residue_heavies[0]}."
        )
        lines.append("")
    else:
        lines.append(
            selection.get("pool_note")
            or (
                "Cycle 1 is the monomer pool; its heavy-atom count is the "
                "sum of both monomers. The first condensation loses a small "
                "molecule, so residue 2 can have fewer heavy atoms than the "
                "pool."
            )
        )
        lines.append("")

    lines.append("### Structures")
    lines.append("")
    prefix = selection.get("prefix") or "step"
    for index, smiles in enumerate(result["intermediates"]):
        label = (
            "monomers"
            if index == 0 and len(smiles) > 1
            else f"residue {index + 1}"
        )
        lines.append(f"#### {label}")
        lines.append("")
        for offset, _smile in enumerate(smiles):
            suffix = (
                ""
                if len(smiles) == 1
                else chr(ord("a") + offset)
            )
            alt = label if len(smiles) == 1 else f"{label} {offset + 1}"
            lines.append(
                f"![{alt}](traces/{prefix}_{index + 1:02d}{suffix}.svg)"
            )
        lines.append("")


_DRIVER_LABELS = {
    "chain": "chain growth",
    "step": "step growth",
    "peo": "epoxide ROP",
    "pcl": "lactone ROP",
    "nylon6": "lactam ROP",
    "pes": "ROCOP",
    "ppo": "oxidative coupling",
    "romp": "ROMP",
}


def _write_section(selection: dict, result: dict, lines: list[str]) -> None:
    if selection["kind"] == "addition":
        _write_chain_section(selection, result, lines)
    elif selection["kind"] == "romp":
        _write_romp_section(selection, result, lines)
    else:
        _write_step_section(selection, result, lines)


def _summary_row(selection: dict) -> str:
    label = _DRIVER_LABELS.get(
        selection.get("prefix", ""), selection["kind"]
    )
    monomers = " + ".join(
        f"`{value}`" for value in selection["monomers"]
    )
    initiator = selection.get("initiator_label") or "-"
    return (
        f"| {label} | {selection['polymer_name']} | {monomers} | "
        f"{initiator} | {selection['heavy']} |"
    )


_MECH_SHORT = {
    "A1": "PP",
    "A2": "PNB-addition",
    "A3": "COC",
    "A4": "EPR",
    "C1": "PET",
    "C2": "Nylon-6,6",
    "C3": "PU",
    "C4": "Kapton",
    "C5": "POX",
    "C6": "PC",
    "C7": "PEEK",
    "C8": "PPO",
    "R1": "PCL",
    "R2": "Nylon-6",
    "R3": "PES",
    "R4": "PEO",
    "R5": "PNB-ROMP",
}


def _ascii_label(text: str) -> str:
    """Keep RDKit SVG legends to characters it can render reliably."""
    return (
        text.replace("\u2013", "-")
        .replace("\u2014", "-")
        .replace("\u2082", "2")
    )


_COMPILED_FIELDS = (
    "init_smarts",
    "prop_smarts",
    "term_smarts",
    "coupling_smarts",
    "ring_open_smarts",
    "ring_bond_smarts",
    "initiator",
    "terminator",
    "reagents",
    "opening_reagents",
)


def _compiled_recipe(reaction: str) -> dict[str, list[str]]:
    """Return the driver kwargs that are compiled for *reaction*."""
    params = reaction_config(reaction)
    compiled: dict[str, list[str]] = {}
    for name in _COMPILED_FIELDS:
        value = params.get(name)
        if value is None:
            continue
        if isinstance(value, str):
            values = [value]
        else:
            values = [str(item) for item in value]
        compiled[name] = values
    return compiled


def _format_compiled_recipe(
    compiled: dict[str, list[str]],
) -> list[str]:
    """Render the compiled driver kwargs as a Markdown bullet list."""
    if not compiled:
        return []
    lines = ["**Compiled recipe (what the driver runs):**", ""]
    for name, values in compiled.items():
        if len(values) == 1:
            lines.append(f"- `{name}`: `{values[0]}`")
        else:
            lines.append(f"- `{name}`:")
            lines.extend(f"  - `{value}`" for value in values)
    lines.append("")
    return lines


def _mechanism_record(
    selection: dict,
    result: dict,
    entries: list[tuple[str, str]],
    trace_images: list[Path],
    representative: Path,
    panel: Path,
) -> dict:
    config = selection["config"]
    return {
        "id": selection["mechanism_id"],
        "family": selection["mechanism_family"],
        "recipe": selection["recipe"],
        "reaction_smarts": selection["reaction"],
        "compiled": _compiled_recipe(selection["reaction"]),
        "registry_keys": sorted(
            key
            for key, config in rxn_map.items()
            if config is not None
            and config.get("mechanism") == selection["mechanism_id"]
        ),
        "registry_type": config["type"],
        "polymer_class": config["polymer_class"],
        "polymer": selection["polymer_name"],
        "monomers": selection["monomers"],
        "residues_or_steps": selection["length_arg"],
        "product": result["final"][0] if result["final"] else None,
        "status": selection.get("status", "commercial"),
        "note": selection.get("mechanism_note", ""),
        "paired": bool(selection.get("paired")),
        "images": {
            "representative": representative.relative_to(
                _BENCH
            ).as_posix(),
            "trace_panel": panel.relative_to(_BENCH).as_posix(),
            "trace": [
                path.relative_to(_BENCH).as_posix()
                for path in trace_images
            ],
        },
        "trace": [
            {"step": label, "smiles": smiles}
            for label, smiles in entries
        ],
    }


def _write_mechanism_markdown(records: list[dict]) -> None:
    lines = [
        "# Mechanism examples (A1–R5)",
        "",
        (
            "One persistent example for each of the 17 distinct mechanisms in "
            "`mechanism_table.tex`, using the exact `rxn_map` recipe assigned "
            "to that mechanism. Every registry entry carries a "
            "`\"mechanism\"` field (A1–R5). Each example uses a known "
            "practical polymer; the trace panel shows the mechanism "
            "walkthrough. If no practical polymer is identified for a "
            "mechanism, the entry is marked `needs candidate` for follow-up."
        ),
        "",
        (
            "**Reading the SMARTS fields.** `reaction SMARTS` is the exact "
            "registry key used to select the polymer family. For coupling and "
            "ring-opening entries it is a topology key, not necessarily the "
            "template that is executed. The `Compiled recipe` block below "
            "each example is the SMARTS the driver actually compiles, "
            "including atom maps and chemoselectivity guards."
        ),
        "",
        "Regenerate with:",
        "",
        "```sh",
        "python benchmarks/extract_extension_traces.py --mode mechanisms",
        "```",
        "",
        (
            "| ID | Family | Recipe / topology | Practical polymer | "
            "Status | Monomers |"
        ),
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for record in records:
        monomers = " + ".join(
            f"`{value}`" for value in record["monomers"]
        ) or "-"
        lines.append(
            f"| {record['id']} | {record['family']} | "
            f"{record['recipe']} | {record['polymer']} | "
            f"{record['status']} | {monomers} |"
        )
    lines.extend(
        [
            "",
            (
                "The overview panel collects the representative structure "
                "from every mechanism:"
            ),
            "",
            "![mechanism overview](traces/mechanisms/practical/overview.svg)",
            "",
        ]
    )
    for record in records:
        monomers = " + ".join(
            f"`{value}`" for value in record["monomers"]
        ) or "-"
        product = (
            f"`{record['product']}`"
            if record.get("product")
            else "-"
        )
        lines.extend(
            [
                f"## {record['id']} — {record['polymer']}",
                "",
                f"**Family:** {record['family']}  ",
                f"**Recipe / topology:** {record['recipe']}",
                "",
                "| field | value |",
                "| --- | --- |",
                f"| mechanism | `{record['id']}` |",
                f"| family | {record['family']} |",
                (
                    "| reaction SMARTS | "
                    f"`{record['reaction_smarts']}` |"
                ),
                f"| registry type | `{record['registry_type']}` |",
                f"| polymer class | `{record['polymer_class']}` |",
                f"| practical polymer | {record['polymer']} |",
                f"| status | {record['status']} |",
                f"| monomer(s) | {monomers} |",
                (
                    "| residues / propagation steps | "
                    f"{record['residues_or_steps']} |"
                ),
                f"| generated product | {product} |",
                "",
                record.get("note", ""),
                "",
                *_format_compiled_recipe(record.get("compiled", {})),
            ]
        )
        registry_keys = record.get("registry_keys", [])
        if len(registry_keys) > 1:
            lines.append(
                f"**All registry keys assigned to {record['id']}:**"
            )
            lines.append("")
            lines.extend(f"- `{key}`" for key in registry_keys)
            lines.append("")
        if record.get("images"):
            lines.extend(
                [
                    (
                        f"![{record['id']} representative]"
                        f"({record['images']['representative']})"
                    ),
                    "",
                    "**Mechanism trace:**",
                    "",
                    (
                        f"![{record['id']} trace]"
                        f"({record['images']['trace_panel']})"
                    ),
                    "",
                    "| step | SMILES |",
                    "| --- | --- |",
                ]
            )
            for step in record["trace"]:
                lines.append(
                    f"| {step['step']} | `{step['smiles']}` |"
                )
            lines.append("")
        else:
            lines.extend(
                [
                    (
                        "**Candidate needed.** No practical polymer is "
                        "identified for this mechanism yet; provide its "
                        "monomers and regenerate."
                    ),
                    "",
                ]
            )
    checks = _run_generalization_checks()
    lines.extend(
        [
            "## Generalization checks",
            "",
            (
                "These small reactions are executed during generation. The "
                "build fails if any observed product/rejection differs from "
                "the expected outcome."
            ),
            "",
            (
                "| Mechanism | Check | Monomers | Length | Expected | "
                "Outcome | Result |"
            ),
            "| --- | --- | --- | --- | --- | --- | --- |",
        ]
    )
    for row in checks:
        monomers = " + ".join(
            f"`{monomer}`" for monomer in row["monomers"]
        )
        result = (
            f"`{row['result']}`"
            if row["outcome"] == "product"
            else row["result"]
        )
        lines.append(
            f"| {row['mechanism']} | {row['note']} | {monomers} | "
            f"{row['length']} | {row['expect']} | {row['outcome']} | "
            f"{result} |"
        )
    lines.extend(
        [
            "",
            "## Generalization scope",
            "",
            (
                "The examples on this page are selected by an explicit "
                "`mechanism` tag, not by first-match SMARTS discovery, so two "
                "templates that match the same monomer do not reroute it. The "
                "compiled recipes already include the atom maps and "
                "chemoselectivity guards that the registry keys do not show. "
                "The remaining limits are monomer-scope choices, not driver "
                "defects."
            ),
            "",
            "| Area | Classification | Detail |",
            "| --- | --- | --- |",
            (
                "| Alcohol/thiol guards | implemented | `_nucleophile` "
                "rejects acyloxy, sulfonic/boronic/phosphoric, "
                "oxime/hydroxylamine, peroxide, and enol-ether O-H. |"
            ),
            (
                "| Amine guards | implemented | `_amine` rejects "
                "amide/imide/thioamide, sulfonamide, and phosphoramide N-H, "
                "but allows amidine/thioimidate. A hydrazide N-H is not "
                "excluded. |"
            ),
            (
                "| C3/C4/C6/R3 atom maps | implemented | The compiled "
                "recipes map every reactive centre. C6 consumes phosgene "
                "`ClC(=O)Cl`; C4 closes a 5-/6-membered imide ring. |"
            ),
            (
                "| R1/R2 cross-matching | not applicable | The drivers "
                "receive only the entry selected by its `mechanism` tag. "
                "The keys also match cyclic anhydrides/imides, but keys are "
                "not used for discovery. |"
            ),
            (
                "| A1 acyclic scope | implemented | The compiled template "
                "requires a non-ring C=C (`!@`) and `X3` on both carbons. "
                "Endocyclic alkenes and allenes/ketenes are rejected; "
                "exocyclic vinyls still propagate. |"
            ),
            (
                "| A1-A4 regio | implemented for terminal vinyls | Each "
                "acyclic unit attaches at its CH2 carbon, giving "
                "head-to-tail enchainment. Internal alkenes without a CH2 "
                "are not matched. |"
            ),
            (
                "| C7 PEEK/PES/PAEN | implemented | NAS queries cover "
                "Ar-F/Cl/Br/I with ketone, sulfone, or nitrile activation, "
                "either on the halide ring or on a directly attached phenyl "
                "or naphthyl spacer. Nitrile activation covers both para and "
                "ortho isomers, including the difunctional "
                "2,6-difluorobenzonitrile PAEN monomer. |"
            ),
            (
                "| C8 PPO | implemented | Both ortho positions of the "
                "coupling ring must be blocked; unsubstituted phenol is "
                "rejected and 2,6-disubstituted phenol couples linearly. |"
            ),
            (
                "| C1 polyester | implemented | Acid/acyl-halide "
                "esterification and methyl-ester transesterification (the "
                "DMT route) are supported. Ethyl-ester transesterification "
                "is not enabled because it cannot be distinguished safely "
                "from backbone esters. |"
            ),
            (
                "| R3 ROCOP | implemented | The balanced opening maps all "
                "atoms and accepts 5-/6-membered anhydrides with C, N, O, or "
                "S ring atoms flanking the carbonyls. |"
            ),
            (
                "| A1-A4 dienes | not modelled | Each propagation inserts "
                "one C=C. Conjugated 1,4-diene addition is outside the "
                "current templates. |"
            ),
            (
                "| C5/R4 epoxide regio | implemented | Four compiled tiers "
                "include CH2/CH2 and CH0/CH0, so ethylene oxide and fully "
                "substituted epoxides match. |"
            ),
            "",
        ]
    )
    _MECH_MD.write_text("\n".join(lines) + "\n")


def _pending_mechanism_record(spec: dict) -> dict:
    """Build a page entry for a mechanism whose candidate is still needed."""
    config = rxn_map[spec["reaction"]]
    return {
        "id": spec["id"],
        "family": spec["family"],
        "recipe": spec["recipe"],
        "reaction_smarts": spec["reaction"],
        "registry_type": config["type"],
        "polymer_class": config["polymer_class"],
        "polymer": spec["polymer"],
        "monomers": [],
        "residues_or_steps": 0,
        "product": None,
        "status": "needs candidate",
        "note": spec["note"],
        "compiled": _compiled_recipe(spec["reaction"]),
        "paired": bool(spec.get("paired")),
        "registry_keys": sorted(
            key
            for key, cfg in rxn_map.items()
            if cfg is not None
            and cfg.get("mechanism") == spec["id"]
        ),
        "images": {},
        "trace": [],
    }


def _write_mechanism_examples() -> None:
    """Generate one persistent example and image set per mechanism ID."""
    ids = [spec["id"] for spec in _MECHANISMS]
    if len(ids) != 17 or len(set(ids)) != 17:
        raise ValueError(f"Expected 17 unique mechanism IDs, got {ids}")
    registry_ids = {
        config.get("mechanism")
        for config in rxn_map.values()
        if config is not None
    }
    if registry_ids != set(ids):
        raise ValueError(
            "rxn_map mechanism codes do not match the table IDs: "
            f"{sorted(registry_ids)} vs {sorted(ids)}"
        )

    _MECH_TRACES.mkdir(parents=True, exist_ok=True)
    for path in _MECH_TRACES.glob("*.svg"):
        path.unlink()

    records = []
    overview_molecules = []
    overview_legends = []
    for spec in _MECHANISMS:
        if (
            spec.get("status") == "needs candidate"
            or not spec.get("monomers")
        ):
            records.append(_pending_mechanism_record(spec))
            continue

        selection = _mechanism_selection(spec)
        result = _trace(selection["kind"], selection)
        _draw_result(selection, result, _MECH_TRACES)
        entries = _trace_entries(selection, result)

        prefix = selection["prefix"]
        trace_images = sorted(
            path
            for path in _MECH_TRACES.glob(f"{prefix}_*.svg")
            if not path.stem.endswith("_trace")
        )
        panel = _MECH_TRACES / f"{prefix}_trace.svg"
        short = _MECH_SHORT.get(spec["id"], spec["id"])
        _draw_grid(
            [
                (f"{short} {label}", smiles)
                for label, smiles in entries
            ],
            panel,
        )

        if not result["final"]:
            raise RuntimeError(
                f"Mechanism {spec['id']} produced no final structure"
            )
        representative = _MECH_TRACES / f"{prefix}.svg"
        representative_mol = _molecule(result["final"][0])
        if representative_mol is None:
            raise RuntimeError(
                f"Mechanism {spec['id']} has no drawable representative"
            )
        _draw(
            representative_mol,
            representative,
            _ascii_label(f"{spec['id']} {spec['polymer']}"),
        )
        overview_molecules.append(representative_mol)
        overview_legends.append(
            _ascii_label(
                f"{spec['id']} {_MECH_SHORT.get(spec['id'], '')}"
            )
        )

        records.append(
            _mechanism_record(
                selection,
                result,
                entries,
                trace_images,
                representative,
                panel,
            )
        )

    overview = _MECH_TRACES / "overview.svg"
    svg = MolsToGridImage(
        overview_molecules,
        molsPerRow=4,
        subImgSize=(400, 300),
        legends=overview_legends,
        useSVG=True,
    )
    if isinstance(svg, bytes):
        svg = svg.decode("utf-8")
    overview.write_text(svg)

    _write_mechanism_markdown(records)
    print(f"Wrote {_MECH_MD}")
    print(
        f"Wrote {len(list(_MECH_TRACES.glob('*.svg')))} SVG files "
        f"for {len(records)} mechanisms"
    )


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--mode",
        choices=("conventional", "simple", "shortest", "mechanisms"),
        default="conventional",
        help=(
            "conventional: PP, PET, PEO, PCL, Nylon-6, PES, PPO and "
            "poly(norbornene) (default); simple: minimal template examples; "
            "shortest: shortest surviving extended-set candidates; "
            "mechanisms: one persistent example per A1–R5 mechanism"
        ),
    )
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    if args.mode == "mechanisms":
        _write_mechanism_examples()
        return

    _TRACES.mkdir(parents=True, exist_ok=True)
    for path in _TRACES.glob("*.svg"):
        path.unlink()

    chain_selection = _select("addition", args.mode)
    step_selection = _select("step_growth", args.mode)
    trace_selections = [chain_selection, step_selection]
    if args.mode == "conventional":
        trace_selections.extend(
            _spec_selection(spec) for spec in _CONVENTIONAL_EXTRA
        )

    results = []
    for selection in trace_selections:
        result = _trace(selection["kind"], selection)
        _draw_result(selection, result)
        results.append((selection, result))

    titles = {
        "conventional": "# Conventional polymer extension traces",
        "simple": (
            "# Minimal template traces (two propagation steps)"
        ),
        "shortest": (
            "# Shortest-candidate traces (two propagation steps)"
        ),
    }
    intros = {
        "conventional": (
            "Eight well-known polymers that are easy to misrepresent by "
            "simply joining repeat-unit wildcards, covering every driver "
            "family and the ring-opening/step-growth coupling variants: "
            "polypropylene (chain growth), poly(ethylene terephthalate) "
            "(step growth), poly(ethylene oxide) (epoxide ring opening), "
            "polycaprolactone (lactone ROP), Nylon-6 (lactam ROP), "
            "poly(ethylene succinate) (anhydride + epoxide ROCOP), "
            "poly(2,6-dimethyl-1,4-phenylene oxide) (oxidative phenol "
            "coupling), and poly(norbornene) (ROMP). The same reaction "
            "templates as the benchmark are used. Chain growth and ROMP "
            "show two propagation steps; the step-growth and ring-opening "
            "traces show three couplings (four residues), and the PET trace "
            "shows four residues so both monomers appear twice."
        ),
        "simple": (
            "Small, chemically clear reference molecules with the same "
            "reaction templates: propene for chain growth and glycolic "
            "acid for step growth."
        ),
        "shortest": (
            "The shortest surviving candidates from the extended test set, "
            "shown with the same reaction templates."
        ),
    }
    lines = [
        titles[args.mode],
        "",
        intros[args.mode],
        "",
        "Regenerate with:",
        "",
        "```sh",
        "python benchmarks/extract_extension_traces.py",
        "python benchmarks/extract_extension_traces.py --mode simple",
        "python benchmarks/extract_extension_traces.py --mode shortest",
        "```",
        "",
        "| driver | polymer | monomer(s) | initiator | heavy atoms |",
        "| --- | --- | --- | --- | --- |",
    ]
    lines.extend(_summary_row(selection) for selection in trace_selections)
    lines.extend(
        [
            "",
            (
                "_SMILES are RDKit canonical strings. In chain-growth and "
                "ROMP intermediates the wildcard `*` marks the active end; "
                "step-growth and ring-opening intermediates are complete "
                "molecules with real monomer end groups._"
            ),
            "",
        ]
    )
    if args.mode == "conventional":
        lines.extend(
            [
                (
                    "Chain growth and ROMP: two propagation steps. Step "
                    "growth: three couplings (four residues) so both PET "
                    "monomers appear twice. Ring-opening traces start at the "
                    "opened residue; the initiator O is the extra oxygen in the "
                    "first residue."
                ),
                "",
            ]
        )
    for selection, result in results:
        _write_section(selection, result, lines)

    _MD.write_text("\n".join(lines) + "\n")
    print(f"Wrote {_MD}")
    print(f"Wrote {len(list(_TRACES.glob('*.svg')))} SVG structures")


if __name__ == "__main__":
    main()
