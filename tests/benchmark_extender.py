"""Benchmark the polymer-extender drivers on the full evaluation set.

The benchmark answers three questions for every iterative cycle ``c``:

1. **Effectiveness** - how many candidates still extend successfully?
2. **Where it fails** - which reaction, row, monomer and cycle broke, and
   with which failure mode (no product, invalid product, invented end group,
   invented atom, no growth, ...)?
3. **What changed** - how do the numbers compare with a saved baseline?

Input is ``tests/test_data/extender_test_data_all.csv`` produced by
``tests/extract_extender_test_data.py``.  Every configured reaction is run at
each requested cycle.  A cycle is one addition/ROMP propagation step for
``addition``/``romp`` entries (cycle 1 is the terminated initiation product)
and one opened monomer residue for ``step_growth``/``ring_opening`` entries.

Checks
------
Structural
    Product parses and sanitises, is a single connected molecule, contains no
    dummy atoms, and grows in heavy-atom count from the previous successful
    cycle.

End-group provenance (step growth)
    The product may only contain terminal functional groups that are already
    present in the (ring-opened) monomers or in the configured reagents.
    This is the check that catches a wildcard-derived cap silently replacing
    an acid chloride with an acid, an isocyanate with an amine, etc.

Atom provenance (step growth, optional)
    Monomer and reagent heavy atoms are isotope-labelled before the run.  Any
    unlabelled heavy atom in the product was created by a reaction template.
    A coupling is allowed only the product-template atoms that are not mapped
    back to a reactant.  The configured ester/amide/lactone templates map
    their carbonyl oxygen through, so the allowance is normally zero.  This
    catches invented caps and invented atoms.

Determinism (optional)
    Running the same row and cycle twice must return the same canonical
    products.

Coverage
    Reactions and rows are partitioned into configured / unconfigured and
    with-data / without-data, so the cost of an unconfigured mechanism is
    visible next to the pass rates.

Usage::

    python tests/benchmark_extender.py
    python tests/benchmark_extender.py --cycles 1,2,3,5,10 --limit 20
    python tests/benchmark_extender.py --baseline benchmarks/results/<run>/summary.json
    python tests/benchmark_extender.py --fail-on-new-failures

Exit codes::

    0  benchmark completed and no optional gate fired
    1  ``--fail-under`` was not met by at least one requested cycle
    2  ``--fail-on-regression`` found a regression
    3  ``--fail-on-new-failures`` found a failure that is not a known gap
"""

from __future__ import annotations

import argparse
import json
import statistics
import sys
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger

try:
    import resource
except ImportError:  # pragma: no cover - non-Unix fallback
    resource = None

_PROJECT = Path(__file__).resolve().parent.parent
if str(_PROJECT) not in sys.path:
    sys.path.insert(0, str(_PROJECT))

from polyscript.extender import (  # noqa: E402
    ChainGrowthReaction,
    RingOpeningMetathesisReaction,
    RingOpeningReaction,
    StepGrowthReaction,
)
from polyscript.extender.base import prepare_monomers  # noqa: E402
from polyscript.extender.map import (  # noqa: E402
    RING_OPENING_PROPAGATIONS,
    STEP_GROWTH_COUPLINGS,
    reaction_config,
    rxn_map,
)

RDLogger.DisableLog("rdApp.*")

_CSV_PATH = (
    _PROJECT / "tests" / "test_data" / "extender_test_data_all.csv"
)
_DEFAULT_RESULTS = _PROJECT / "benchmarks" / "results"
_DEFAULT_CYCLES = (1, 2, 3, 5, 10)

# ---------------------------------------------------------------------------
# failure taxonomy
# ---------------------------------------------------------------------------
OK = ""
FAIL_CONFIG = "config_error"
FAIL_EXCEPTION = "exception"
FAIL_NO_PRODUCT = "no_product"
FAIL_INVALID = "invalid_product"
FAIL_WILDCARD = "wildcard_product"
FAIL_DISCONNECTED = "disconnected_product"
FAIL_INVENTED_END_GROUP = "invented_end_group"
FAIL_INVENTED_ATOM = "invented_atom"
FAIL_NO_GROWTH = "no_growth"
FAIL_NONDETERMINISTIC = "nondeterministic"

_HARD_FAILURES = {
    FAIL_CONFIG,
    FAIL_EXCEPTION,
    FAIL_NO_PRODUCT,
    FAIL_INVALID,
    FAIL_WILDCARD,
    FAIL_DISCONNECTED,
    FAIL_INVENTED_END_GROUP,
    FAIL_INVENTED_ATOM,
    FAIL_NO_GROWTH,
    FAIL_NONDETERMINISTIC,
}

# ---------------------------------------------------------------------------
# terminal-group catalogue
# ---------------------------------------------------------------------------
# These are *terminal* groups only.  Backbone groups that a coupling creates
# (ester, amide, urethane, imide, ...) are deliberately absent, because they
# are not terminal chemistry and would make every successful coupling look
# like an invented end group.
END_GROUP_SMARTS = {
    "carboxylic_acid": "[#6&X3](=[#8&X1])[#8&X2&H1]",
    "acid_halide": "[#6&X3](=[#8&X1])[F,Cl,Br,I]",
    "alcohol": "[#8&X2&H1;!$(O[#6&X3]=[#8&X1])]",
    "phenol": "[c][#8&X2&H1]",
    "thiol": "[#16&X2&H1;!$(S[#6&X3]=[#8&X1])]",
    "primary_amine": "[#7&X3;H2;!$(N[#6&X3]=[#8,#16&X1])]",
    "secondary_amine": "[#7&X3;H1;!$(N[#6&X3]=[#8,#16&X1])]",
    "isocyanate": "[#7&X2]=[#6&X2]=[#8&X1]",
    "isothiocyanate": "[#7&X2]=[#6&X2]=[#16&X1]",
    "aryl_halide": "[c][F,Cl,Br,I]",
    "epoxide": "[#6&X4;R]1[#8&X2;R][#6&X4;R]1",
    "anhydride": (
        "[#6&X3](=[#8&X1])[#8&X2][#6&X3]=[#8&X1]"
    ),
}
_END_GROUP_QUERIES = {
    name: Chem.MolFromSmarts(smarts)
    for name, smarts in END_GROUP_SMARTS.items()
}

# Legitimate new terminal groups created by a coupling rather than carried
# over from a monomer.
#
# * ``epoxide_rop`` opens the epoxide with water, so the product's two -OH
#   ends do not exist in the raw epoxide monomer: one O is the water initiator
#   and the other is the epoxide ring O after opening.  Neither is an invented
#   cap; the exception documents the ring-opening event the step-growth
#   provenance model cannot infer from the intact monomer.
# * ``lactone_rop`` on a cyclic carbamate/carbonic anhydride (an NCA-like
#   oxazinedione) opens to an acyloxy O-H.  Coupling that O-H to the next
#   opened acid necessarily forms a mixed anhydride linkage and retains the
#   carbamate carbonyl.  This is the documented topological match for that
#   monomer family; no new terminal group is invented.
ALLOWED_NEW_END_GROUPS = {
    "epoxide_rop": {"alcohol"},
    "lactone_rop": {"anhydride"},
}

_COUPLING_NAMES = {}
for _name, _smarts in STEP_GROWTH_COUPLINGS.items():
    for _alternative in (
        _smarts if isinstance(_smarts, list) else [_smarts]
    ):
        _COUPLING_NAMES[_alternative] = _name
for _name, _smarts in RING_OPENING_PROPAGATIONS.items():
    # ``anhydride_epoxide_rop`` reuses the ordinary ester template, so keep
    # the canonical step-growth name for that shared SMARTS.
    for _alternative in (
        _smarts if isinstance(_smarts, list) else [_smarts]
    ):
        _COUPLING_NAMES.setdefault(_alternative, _name)

# Ring-opening propagation makes the two previously known chemistry gaps
# pass.  Keep the mapping explicit so a regression is reported as a new
# failure rather than being silently excused.
KNOWN_GAPS = {}

# Configurations whose products are checked for terminal-group and atom
# provenance.  ROMP is chain growth and takes its end groups from its
# initiator/terminator fragments, so it uses the addition semantics.
_STEP_LIKE_TYPES = {"step_growth", "ring_opening"}
_ADDITION_LIKE_TYPES = {"addition", "romp"}


# ---------------------------------------------------------------------------
# small helpers
# ---------------------------------------------------------------------------
def _monomers(row) -> list[str]:
    values = []
    for column in ("mon1", "mon2"):
        value = row.get(column)
        if isinstance(value, str) and value.strip():
            values.append(value)
    return values


def _known_gap(monomers: list[str]) -> str | None:
    for monomer in monomers:
        if monomer in KNOWN_GAPS:
            return KNOWN_GAPS[monomer]
    return None


def _end_group_counts(mol) -> dict[str, int]:
    counts = {}
    for name, query in _END_GROUP_QUERIES.items():
        matches = len(mol.GetSubstructMatches(query))
        if matches:
            counts[name] = matches
    return counts


def _end_group_labels(mol) -> set[str]:
    return set(_end_group_counts(mol))


def _dedupe_mols(mols: list) -> list:
    unique = []
    seen = set()
    for mol in mols:
        smiles = Chem.MolToSmiles(mol)
        if smiles in seen:
            continue
        seen.add(smiles)
        unique.append(mol)
    return unique


def _strip_isotopes(mol):
    rw = Chem.RWMol(mol)
    for atom in rw.GetAtoms():
        atom.SetIsotope(0)
    return rw.GetMol()


def _label_heavy_atoms(mol, base: int):
    rw = Chem.RWMol(mol)
    for index, atom in enumerate(rw.GetAtoms()):
        if atom.GetAtomicNum() > 1:
            atom.SetIsotope(base + index)
    return rw.GetMol()


def _unlabeled_heavy_atoms(mol) -> Counter:
    return Counter(
        atom.GetSymbol()
        for atom in mol.GetAtoms()
        if atom.GetAtomicNum() > 1 and atom.GetIsotope() == 0
    )


def _percentile(values: list[float], fraction: float) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    position = round(fraction * (len(ordered) - 1))
    return ordered[position]


def _median(values: list[float]) -> float | None:
    if not values:
        return None
    return statistics.median(values)


def _json_compact(value) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True)


def _parse_cycles(text: str) -> tuple[int, ...]:
    cycles = sorted({int(part) for part in text.split(",") if part.strip()})
    if not cycles or cycles[0] < 1:
        raise argparse.ArgumentTypeError("cycles must be positive integers")
    return tuple(cycles)


# ---------------------------------------------------------------------------
# driver construction and provenance
# ---------------------------------------------------------------------------
def _build_driver(reaction: str, config: dict, max_products: int):
    rxn_type = config["type"]
    if rxn_type == "addition":
        return ChainGrowthReaction(
            **config["params"], max_products=max_products
        )
    if rxn_type == "step_growth":
        return StepGrowthReaction(
            **reaction_config(reaction), max_products=max_products
        )
    if rxn_type == "ring_opening":
        return RingOpeningReaction(
            **reaction_config(reaction), max_products=max_products
        )
    if rxn_type == "romp":
        return RingOpeningMetathesisReaction(
            **reaction_config(reaction), max_products=max_products
        )
    raise ValueError(f"Unknown reaction type: {rxn_type!r}")


def _build_probe_driver(reaction: str, config: dict, max_products: int):
    """Driver used for isotope provenance.

    Ring-opening is done by the caller before labelling, so the probe driver
    must not open the monomer a second time.  Reagents are moved into the
    labelled pool instead of being injected by the driver.
    """
    params = dict(reaction_config(reaction))
    params.pop("ring_open_smarts", None)
    params.pop("opening_reagents", None)
    prop_smarts = params.pop("prop_smarts", None)
    if prop_smarts is not None:
        params["coupling_smarts"] = prop_smarts
    params["reagents"] = []
    return StepGrowthReaction(**params, max_products=max_products)


def _step_growth_pool(primary, monomers: list[str], reagents: list[str]):
    """Return the real (ring-opened) monomer/reagent pool."""
    if hasattr(primary, "prepare_pool"):
        return _dedupe_mols(primary.prepare_pool(prepare_monomers(monomers)))
    pool = prepare_monomers(monomers)
    if primary.ring_open_rxn is not None:
        pool = [primary._open_monomer(Chem.Mol(mol)) for mol in pool]
    pool += prepare_monomers(reagents)
    return _dedupe_mols(pool)


def _labeled_pool(primary, monomers: list[str], reagents: list[str]):
    pool = _step_growth_pool(primary, monomers, reagents)
    return [
        _label_heavy_atoms(mol, 1000 * (index + 1))
        for index, mol in enumerate(pool)
    ]


def _allowed_new_atoms(probe) -> Counter:
    """Element-wise maximum number of unmapped product atoms per coupling."""
    per_coupling: Counter = Counter()
    for reaction in probe.coupling_rxns:
        local: Counter = Counter()
        for template_index in range(reaction.GetNumProductTemplates()):
            template = reaction.GetProductTemplate(template_index)
            for atom in template.GetAtoms():
                if atom.GetAtomMapNum() == 0 and atom.GetAtomicNum() > 1:
                    local[atom.GetSymbol()] += 1
        for element, count in local.items():
            per_coupling[element] = max(
                per_coupling[element], count
            )
    return per_coupling


def _coupling_names(reaction: str, config: dict) -> list[str]:
    names = []
    params = config.get("params", {})
    for smarts in list(params.get("coupling_smarts", [])) + list(
        params.get("prop_smarts", [])
    ):
        name = _COUPLING_NAMES.get(smarts)
        if name and name not in names:
            names.append(name)
    return names


def _allowed_new_end_groups(config: dict) -> set[str]:
    allowed: set[str] = set()
    for name in _coupling_names("", config):
        allowed |= ALLOWED_NEW_END_GROUPS.get(name, set())
    return allowed


# ---------------------------------------------------------------------------
# diagnosis
# ---------------------------------------------------------------------------
def _diagnose_step_growth(
    primary, reaction: str, config: dict, monomers: list[str], cycle: int
) -> str:
    """Return a short hint about why a step-growth cycle failed."""
    try:
        if primary.ring_open_rxn is not None:
            if primary.ring_open_rxn.GetNumReactantTemplates() == 1:
                query = primary.ring_open_rxn.GetReactantTemplate(0)
                unmatched = [
                    Chem.MolToSmiles(mol)
                    for mol in prepare_monomers(monomers)
                    if not mol.HasSubstructMatch(query)
                ]
                if unmatched:
                    return "ring_open_no_match: " + "; ".join(
                        smile[:60] for smile in unmatched
                    )

        if cycle <= 1:
            return "empty monomer pool after ring opening"

        try:
            previous = primary.extend(monomers, cycle - 1)
        except Exception as exc:  # noqa: BLE001 - diagnosis only
            return f"previous_cycle_failed: {type(exc).__name__}: {exc}"
        if not previous:
            return "previous_cycle_empty"

        chain = Chem.MolFromSmiles(previous[0])
        if chain is None:
            return "previous_cycle_product_unparsable"
        pool = _step_growth_pool(
            primary,
            monomers,
            config.get("params", {}).get("reagents", []),
        )
        for index, coupling in enumerate(primary.coupling_rxns):
            if primary._run_coupling(coupling, chain, pool):
                return (
                    f"coupling[{index}] matched but extension failed "
                    "(frontier/trim)"
                )
        return (
            "coupling_no_match: none of "
            f"{len(primary.coupling_rxns)} coupling template(s) matched "
            "the chain and monomer pool"
        )
    except Exception as exc:  # noqa: BLE001 - diagnosis must never raise
        return f"diagnosis_failed: {type(exc).__name__}: {exc}"


def _diagnose_chain_growth(
    driver, monomers: list[str], cycle: int
) -> str:
    """Return a short hint about which chain-growth stage failed."""
    try:
        chains = driver.chain_init(monomers)
        if not chains:
            return "init_no_match"
        for _ in range(max(0, cycle - 1)):
            chains = driver.chain_extend(monomers, chains)
            if not chains:
                return "propagation_no_match"
        terminated = driver.chain_term(monomers, chains)
        if not terminated:
            return "termination_no_match"
    except Exception as exc:  # noqa: BLE001 - diagnosis must never raise
        return f"diagnosis_failed: {type(exc).__name__}: {exc}"
    return "unknown"


# ---------------------------------------------------------------------------
# one row / one cycle
# ---------------------------------------------------------------------------
def _record(
    *,
    reaction,
    config,
    row,
    cycle,
    length_arg,
    ok,
    failure,
    detail,
    hint,
    warnings,
    products,
    heavy,
    rings,
    fragments,
    dummies,
    end_groups,
    new_end_groups,
    allowed_new_end_groups,
    missing_end_groups,
    unlabeled,
    allowed_new_atoms,
    excess_new_atoms,
    delta_heavy,
    delta_rings,
    runtime_ms,
    deterministic,
    known_gap,
) -> dict:
    return {
        "reaction_smarts": reaction,
        "polymer_class": config.get("polymer_class", ""),
        "reaction_type": config.get("type", ""),
        "source_file": row.get("source_file", ""),
        "row_index": row.get("row_index", ""),
        "mon1": row.get("mon1", ""),
        "mon2": "" if pd.isna(row.get("mon2", "")) else row.get("mon2", ""),
        "cycle": cycle,
        "length_arg": length_arg,
        "ok": bool(ok),
        "failure": failure,
        "detail": detail,
        "hint": hint,
        "warnings": list(warnings),
        "n_products": len(products) if products else 0,
        "product": products[0] if products else "",
        "heavy_atoms": heavy,
        "rings": rings,
        "fragments": fragments,
        "dummies": dummies,
        "end_groups": _json_compact(end_groups),
        "new_end_groups": _json_compact(sorted(new_end_groups)),
        "allowed_new_end_groups": _json_compact(
            sorted(allowed_new_end_groups)
        ),
        "missing_end_groups": _json_compact(sorted(missing_end_groups)),
        "unlabeled_atoms": _json_compact(dict(unlabeled)),
        "allowed_new_atoms": _json_compact(dict(allowed_new_atoms)),
        "excess_new_atoms": _json_compact(dict(excess_new_atoms)),
        "delta_heavy": delta_heavy,
        "delta_rings": delta_rings,
        "runtime_ms": round(runtime_ms, 3),
        "deterministic": deterministic,
        "known_gap": known_gap or "",
    }


def _probe_atom_provenance(
    primary,
    probe,
    monomers: list[str],
    reagents: list[str],
    cycle: int,
    per_coupling_atoms: Counter,
) -> dict | None:
    """Run the cycle with labelled atoms and return provenance findings."""
    if probe is None:
        return None

    try:
        pool = _labeled_pool(primary, monomers, reagents)
        products = probe.extend(pool, cycle)
    except Exception as exc:  # noqa: BLE001 - recorded as a warning
        return {"error": f"{type(exc).__name__}: {exc}"}

    if not products:
        return {"error": "probe_no_product"}

    mol = Chem.MolFromSmiles(products[0])
    if mol is None:
        return {"error": "probe_product_unparsable"}

    unlabeled = _unlabeled_heavy_atoms(mol)
    allowed = Counter(
        {
            element: count * max(0, cycle - 1)
            for element, count in per_coupling_atoms.items()
        }
    )
    excess = {
        element: count - allowed.get(element, 0)
        for element, count in unlabeled.items()
        if count > allowed.get(element, 0)
    }
    return {
        "unlabeled": dict(unlabeled),
        "allowed": dict(allowed),
        "excess": excess,
        "product": Chem.MolToSmiles(_strip_isotopes(mol)),
    }


def _run_row_cycle(
    *,
    reaction: str,
    config: dict,
    row,
    cycle: int,
    driver,
    probe_driver,
    reagents: list[str],
    provenance_labels: set[str],
    allowed_labels: set[str],
    per_coupling_atoms: Counter,
    max_products: int,
    atom_provenance: bool,
    check_determinism: bool,
) -> dict:
    monomers = _monomers(row)
    known_gap = _known_gap(monomers)
    length_arg = (
        cycle - 1
        if config.get("type") in _ADDITION_LIKE_TYPES
        else cycle
    )
    is_step_growth = config.get("type") in _STEP_LIKE_TYPES

    started = time.perf_counter()
    try:
        products = driver.extend(monomers, length_arg)
        runtime_ms = (time.perf_counter() - started) * 1000.0
    except Exception as exc:  # noqa: BLE001 - one row must not stop the run
        runtime_ms = (time.perf_counter() - started) * 1000.0
        hint = ""
        if config.get("type") in _STEP_LIKE_TYPES:
            hint = _diagnose_step_growth(
                driver, reaction, config, monomers, cycle
            )
        else:
            hint = _diagnose_chain_growth(driver, monomers, cycle)
        message = f"{type(exc).__name__}: {exc}"
        failure = (
            FAIL_NO_PRODUCT
            if "couplings failed" in str(exc).lower()
            else FAIL_EXCEPTION
        )
        return _record(
            reaction=reaction,
            config=config,
            row=row,
            cycle=cycle,
            length_arg=length_arg,
            ok=False,
            failure=failure,
            detail=message,
            hint=hint,
            warnings=[],
            products=[],
            heavy=None,
            rings=None,
            fragments=None,
            dummies=None,
            end_groups={},
            new_end_groups=set(),
            allowed_new_end_groups=set(),
            missing_end_groups=set(),
            unlabeled={},
            allowed_new_atoms={},
            excess_new_atoms={},
            delta_heavy=None,
            delta_rings=None,
            runtime_ms=runtime_ms,
            deterministic=None,
            known_gap=known_gap,
        )

    if not products:
        hint = (
            _diagnose_step_growth(
                driver, reaction, config, monomers, cycle
            )
            if config.get("type") in _STEP_LIKE_TYPES
            else _diagnose_chain_growth(driver, monomers, cycle)
        )
        return _record(
            reaction=reaction,
            config=config,
            row=row,
            cycle=cycle,
            length_arg=length_arg,
            ok=False,
            failure=FAIL_NO_PRODUCT,
            detail="driver returned no products",
            hint=hint,
            warnings=[],
            products=[],
            heavy=None,
            rings=None,
            fragments=None,
            dummies=None,
            end_groups={},
            new_end_groups=set(),
            allowed_new_end_groups=set(),
            missing_end_groups=set(),
            unlabeled={},
            allowed_new_atoms={},
            excess_new_atoms={},
            delta_heavy=None,
            delta_rings=None,
            runtime_ms=runtime_ms,
            deterministic=None,
            known_gap=known_gap,
        )

    deterministic = None
    if check_determinism:
        try:
            second = driver.extend(monomers, length_arg)
            deterministic = second == products
        except Exception:  # noqa: BLE001 - determinism is best effort
            deterministic = False

    first = products[0]
    mol = Chem.MolFromSmiles(first)
    if mol is None:
        return _record(
            reaction=reaction,
            config=config,
            row=row,
            cycle=cycle,
            length_arg=length_arg,
            ok=False,
            failure=FAIL_INVALID,
            detail=f"product is not a valid SMILES: {first[:160]}",
            hint="",
            warnings=[],
            products=[],
            heavy=None,
            rings=None,
            fragments=None,
            dummies=None,
            end_groups={},
            new_end_groups=set(),
            allowed_new_end_groups=set(),
            missing_end_groups=set(),
            unlabeled={},
            allowed_new_atoms={},
            excess_new_atoms={},
            delta_heavy=None,
            delta_rings=None,
            runtime_ms=runtime_ms,
            deterministic=deterministic,
            known_gap=known_gap,
        )

    dummies = sum(
        1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 0
    )
    fragments = len(Chem.GetMolFrags(mol))
    heavy = mol.GetNumHeavyAtoms()
    rings = mol.GetRingInfo().NumRings()

    # End-group provenance.  Addition end groups come from init/term SMARTS,
    # not from the monomers, so this check applies to step growth only.
    end_groups = _end_group_counts(mol)
    product_labels = set(end_groups)
    if is_step_growth:
        new_end_groups = (
            product_labels - provenance_labels - allowed_labels
        )
        allowed_new_end_groups = (
            product_labels - provenance_labels
        ) & allowed_labels
        missing_end_groups = provenance_labels - product_labels
    else:
        new_end_groups = set()
        allowed_new_end_groups = set()
        missing_end_groups = set()

    # Atom provenance (step growth only).
    unlabeled: dict = {}
    allowed_atoms: dict = {}
    excess_atoms: dict = {}
    atom_warning = ""
    if atom_provenance and is_step_growth:
        probe = _probe_atom_provenance(
            driver,
            probe_driver,
            monomers,
            reagents,
            cycle,
            per_coupling_atoms,
        )
        if probe is None:
            pass
        elif "error" in probe:
            atom_warning = probe["error"]
        else:
            unlabeled = probe["unlabeled"]
            allowed_atoms = probe["allowed"]
            excess_atoms = probe["excess"]

    failure = OK
    detail = ""
    hint = ""
    if dummies:
        failure = FAIL_WILDCARD
        detail = f"{dummies} dummy atom(s) remain in the product"
    elif fragments != 1:
        failure = FAIL_DISCONNECTED
        detail = f"product has {fragments} disconnected fragments"
    elif new_end_groups:
        failure = FAIL_INVENTED_END_GROUP
        detail = (
            "product contains terminal group(s) absent from the monomer "
            f"pool: {sorted(new_end_groups)}"
        )
        hint = (
            "monomer terminal groups: "
            f"{sorted(provenance_labels) or 'none'}"
        )
    elif excess_atoms:
        failure = FAIL_INVENTED_ATOM
        detail = (
            "product contains heavy atom(s) not inherited from monomers "
            f"or reagents: {excess_atoms}"
        )
        hint = f"unlabeled={unlabeled} allowed={allowed_atoms}"
    elif check_determinism and deterministic is False:
        failure = FAIL_NONDETERMINISTIC
        detail = "repeated extension returned different products"

    warnings = []
    if atom_warning:
        warnings.append(f"atom_provenance:{atom_warning}")
    if missing_end_groups and not product_labels:
        warnings.append("no_end_groups_after_growth")

    return _record(
        reaction=reaction,
        config=config,
        row=row,
        cycle=cycle,
        length_arg=length_arg,
        ok=failure not in _HARD_FAILURES,
        failure=failure,
        detail=detail,
        hint=hint,
        warnings=warnings,
        products=products,
        heavy=heavy,
        rings=rings,
        fragments=fragments,
        dummies=dummies,
        end_groups=end_groups,
        new_end_groups=new_end_groups,
        allowed_new_end_groups=allowed_new_end_groups,
        missing_end_groups=missing_end_groups,
        unlabeled=unlabeled,
        allowed_new_atoms=allowed_atoms,
        excess_new_atoms=excess_atoms,
        delta_heavy=None,
        delta_rings=None,
        runtime_ms=runtime_ms,
        deterministic=deterministic,
        known_gap=known_gap,
    )


# ---------------------------------------------------------------------------
# aggregation
# ---------------------------------------------------------------------------
def _summarise_cycle(records: list[dict]) -> dict:
    attempted = len(records)
    passed = sum(1 for record in records if record["ok"])
    failures = Counter(
        record["failure"] for record in records if not record["ok"]
    )
    warnings = Counter(
        warning
        for record in records
        for warning in record["warnings"]
    )
    runtimes = [record["runtime_ms"] for record in records]
    product_counts = [record["n_products"] for record in records]
    heavy = [
        record["heavy_atoms"]
        for record in records
        if record["heavy_atoms"] is not None
    ]
    deltas = [
        record["delta_heavy"]
        for record in records
        if record["delta_heavy"] is not None
    ]
    return {
        "attempted": attempted,
        "passed": passed,
        "failed": attempted - passed,
        "pass_rate": round(passed / attempted, 4) if attempted else None,
        "failures": dict(failures),
        "warnings": dict(warnings),
        "invented_end_groups": sum(
            record["failure"] == FAIL_INVENTED_END_GROUP
            for record in records
        ),
        "invented_atoms": sum(
            record["failure"] == FAIL_INVENTED_ATOM
            for record in records
        ),
        "known_gap_failures": sum(
            bool(record["known_gap"]) and not record["ok"]
            for record in records
        ),
        "new_failures": sum(
            not record["ok"] and not record["known_gap"]
            for record in records
        ),
        "median_runtime_ms": _median(runtimes),
        "p95_runtime_ms": _percentile(runtimes, 0.95),
        "mean_products": (
            round(sum(product_counts) / len(product_counts), 3)
            if product_counts
            else None
        ),
        "median_products": _median(product_counts),
        "max_products_observed": (
            max(product_counts) if product_counts else None
        ),
        "rows_with_multiple_products": sum(
            count > 1 for count in product_counts
        ),
        "median_heavy_atoms": _median(heavy),
        "median_delta_heavy": _median(deltas),
    }


def _summarise(records: list[dict], cycles: tuple[int, ...]) -> dict:
    by_cycle = {
        str(cycle): _summarise_cycle(
            [record for record in records if record["cycle"] == cycle]
        )
        for cycle in cycles
    }

    by_class_cycle: dict[str, dict[str, dict]] = defaultdict(dict)
    for cycle in cycles:
        for polymer_class in sorted(
            {record["polymer_class"] for record in records}
        ):
            subset = [
                record
                for record in records
                if record["cycle"] == cycle
                and record["polymer_class"] == polymer_class
            ]
            if subset:
                by_class_cycle[polymer_class][str(cycle)] = (
                    _summarise_cycle(subset)
                )

    # Per-row survival and first-failure onset.
    rows: dict[tuple, list[dict]] = defaultdict(list)
    for record in records:
        key = (
            record["reaction_smarts"],
            record["source_file"],
            record["row_index"],
        )
        rows[key].append(record)

    onset = Counter()
    survival = {str(cycle): 0 for cycle in cycles}
    fully_ok = 0
    for row_records in rows.values():
        ordered = sorted(row_records, key=lambda item: item["cycle"])
        failure_cycles = [
            record["cycle"] for record in ordered if not record["ok"]
        ]
        if failure_cycles:
            onset[str(min(failure_cycles))] += 1
        else:
            fully_ok += 1
        for cycle in cycles:
            if all(
                record["ok"]
                for record in ordered
                if record["cycle"] <= cycle
            ):
                survival[str(cycle)] += 1

    # Per-reaction detail.
    reaction_ids = {
        reaction: f"R{index:02d}"
        for index, reaction in enumerate(
            sorted({record["reaction_smarts"] for record in records}),
            start=1,
        )
    }
    by_reaction: dict[str, dict] = {}
    for reaction in sorted({record["reaction_smarts"] for record in records}):
        subset = [
            record
            for record in records
            if record["reaction_smarts"] == reaction
        ]
        config = rxn_map.get(reaction, {}) or {}
        row_keys = {
            (record["source_file"], record["row_index"])
            for record in subset
        }
        reaction_onsets = []
        for row_key in row_keys:
            row_records = [
                record
                for record in subset
                if (record["source_file"], record["row_index"])
                == row_key
            ]
            failure_cycles = [
                record["cycle"]
                for record in row_records
                if not record["ok"]
            ]
            reaction_onsets.append(
                min(failure_cycles) if failure_cycles else None
            )
        by_reaction[reaction] = {
            "id": reaction_ids[reaction],
            "reaction_smarts": reaction,
            "polymer_class": config.get("polymer_class", ""),
            "type": config.get("type", ""),
            "rows": len(row_keys),
            "fully_ok_rows": sum(
                onset_value is None for onset_value in reaction_onsets
            ),
            "first_failure_cycle": dict(Counter(
                str(value)
                for value in reaction_onsets
                if value is not None
            )),
            "per_cycle": {
                str(cycle): _summarise_cycle(
                    [
                        record
                        for record in subset
                        if record["cycle"] == cycle
                    ]
                )
                for cycle in cycles
            },
            "couplings": _coupling_names(reaction, config),
            "ring_open": bool(
                config.get("params", {}).get("ring_open")
            )
            or config.get("type")
            in {
                "ring_opening",
                "romp",
            },
            "reagents": list(
                config.get("params", {}).get("reagents", [])
            )
            + list(
                config.get("params", {}).get("opening_reagents", [])
            ),
            "allowed_new_end_groups": sorted(
                _allowed_new_end_groups(config)
            ),
        }

    failure_kinds = Counter(
        record["failure"] for record in records if not record["ok"]
    )
    warning_kinds = Counter(
        warning
        for record in records
        for warning in record["warnings"]
    )
    failure_examples = [
        {
            "reaction_id": reaction_ids.get(
                record["reaction_smarts"], ""
            ),
            "reaction_smarts": record["reaction_smarts"],
            "polymer_class": record["polymer_class"],
            "source_file": record["source_file"],
            "row_index": record["row_index"],
            "mon1": record["mon1"],
            "mon2": record["mon2"],
            "cycle": record["cycle"],
            "failure": record["failure"],
            "detail": record["detail"],
            "hint": record["hint"],
            "product": record["product"],
            "known_gap": record["known_gap"],
        }
        for record in records
        if not record["ok"]
    ]
    new_end_group_findings = [
        {
            "reaction_smarts": record["reaction_smarts"],
            "polymer_class": record["polymer_class"],
            "source_file": record["source_file"],
            "row_index": record["row_index"],
            "cycle": record["cycle"],
            "new_end_groups": json.loads(record["new_end_groups"]),
            "allowed": json.loads(record["allowed_new_end_groups"]),
            "mon1": record["mon1"],
            "mon2": record["mon2"],
            "product": record["product"],
        }
        for record in records
        if record["new_end_groups"] != "[]"
    ]
    atom_findings = [
        {
            "reaction_smarts": record["reaction_smarts"],
            "polymer_class": record["polymer_class"],
            "source_file": record["source_file"],
            "row_index": record["row_index"],
            "cycle": record["cycle"],
            "unlabeled": json.loads(record["unlabeled_atoms"]),
            "allowed": json.loads(record["allowed_new_atoms"]),
            "excess": json.loads(record["excess_new_atoms"]),
            "product": record["product"],
        }
        for record in records
        if record["excess_new_atoms"] != "{}"
    ]
    allowed_end_group_findings = [
        {
            "reaction_smarts": record["reaction_smarts"],
            "polymer_class": record["polymer_class"],
            "source_file": record["source_file"],
            "row_index": record["row_index"],
            "cycle": record["cycle"],
            "allowed_new_end_groups": json.loads(
                record["allowed_new_end_groups"]
            ),
            "product": record["product"],
        }
        for record in records
        if record["allowed_new_end_groups"] != "[]"
    ]

    return {
        "cycles": {str(cycle): by_cycle[str(cycle)] for cycle in cycles},
        "by_class_cycle": {
            polymer_class: {
                cycle: summary
                for cycle, summary in cycle_map.items()
            }
            for polymer_class, cycle_map in by_class_cycle.items()
        },
        "by_reaction": by_reaction,
        "failure_kinds": dict(failure_kinds),
        "warning_kinds": dict(warning_kinds),
        "failure_onset": dict(onset),
        "survival": survival,
        "rows_fully_ok": fully_ok,
        "rows_evaluated": len(rows),
        "failure_examples": failure_examples,
        "new_end_group_findings": new_end_group_findings,
        "allowed_end_group_findings": allowed_end_group_findings,
        "atom_provenance_findings": atom_findings,
    }


def _coverage(rows: pd.DataFrame) -> dict:
    reactions_with_data = set(rows["reaction_smarts"].unique())
    configured_reactions = {
        key for key, config in rxn_map.items() if config is not None
    }
    unconfigured_reactions = {
        key for key, config in rxn_map.items() if config is None
    }

    row_reactions = rows["reaction_smarts"]
    configured_mask = row_reactions.isin(configured_reactions)
    unconfigured_mask = row_reactions.isin(unconfigured_reactions)

    class_rows = Counter(rows["polymer_class"].astype(str))
    class_configured = Counter(
        rows.loc[configured_mask, "polymer_class"].astype(str)
    )
    rows_by_reaction = Counter(rows["reaction_smarts"])
    class_by_reaction = {
        reaction: (config or {}).get("polymer_class", "")
        for reaction, config in rxn_map.items()
    }
    for reaction, group in rows.groupby("reaction_smarts"):
        if not class_by_reaction.get(reaction):
            class_by_reaction[str(reaction)] = str(
                group["polymer_class"].iloc[0]
            )

    return {
        "registry_reactions": len(rxn_map),
        "configured_reactions": len(configured_reactions),
        "unconfigured_reactions": len(unconfigured_reactions),
        "reactions_with_data": len(reactions_with_data),
        "configured_with_data": len(
            configured_reactions & reactions_with_data
        ),
        "unconfigured_with_data": sorted(
            unconfigured_reactions & reactions_with_data
        ),
        "configured_without_data": sorted(
            configured_reactions - reactions_with_data
        ),
        "unconfigured_without_data": sorted(
            unconfigured_reactions - reactions_with_data
        ),
        "rows_total": int(len(rows)),
        "rows_configured": int(configured_mask.sum()),
        "rows_unconfigured": int(unconfigured_mask.sum()),
        "rows_by_class": dict(class_rows),
        "configured_rows_by_class": dict(class_configured),
        "rows_by_reaction": dict(rows_by_reaction),
        "class_by_reaction": class_by_reaction,
    }


# ---------------------------------------------------------------------------
# baseline
# ---------------------------------------------------------------------------
def _record_key(record: dict) -> str:
    return "|".join(
        (
            str(record["source_file"]),
            str(record["row_index"]),
            str(record["cycle"]),
        )
    )


def _compare_baseline(records: list[dict], baseline_path: Path) -> dict:
    path = Path(baseline_path)
    if path.is_dir():
        path = path / "summary.json"
    baseline = json.loads(path.read_text())
    old = baseline.get("record_index", {})
    current = {_record_key(record): record for record in records}

    regressions = []
    fixes = []
    new_rows = []
    regressions_by_cycle: Counter = Counter()
    fixes_by_cycle: Counter = Counter()
    for key, record in current.items():
        previous = old.get(key)
        if previous is None:
            new_rows.append(key)
            continue
        if record["ok"] and not previous["ok"]:
            fixes.append(key)
            fixes_by_cycle[str(record["cycle"])] += 1
        elif not record["ok"] and previous["ok"]:
            regressions.append(key)
            regressions_by_cycle[str(record["cycle"])] += 1

    removed = sorted(set(old) - set(current))
    return {
        "baseline": str(path),
        "regressions": regressions,
        "regressions_by_cycle": dict(regressions_by_cycle),
        "fixes": fixes,
        "fixes_by_cycle": dict(fixes_by_cycle),
        "new_rows": new_rows,
        "removed_rows": removed,
        "unchanged_failures": sorted(
            key
            for key, record in current.items()
            if not record["ok"]
            and key in old
            and not old[key]["ok"]
        ),
    }


# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------
def _format_rate(value) -> str:
    if value is None:
        return "-"
    return f"{100 * value:.1f}%"


def _format_number(value, digits: int = 2) -> str:
    if value is None:
        return "-"
    return f"{value:.{digits}f}"


def _format_table(headers: list[str], rows: list[list]) -> str:
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(str(cell) for cell in row) + " |")
    return "\n".join(lines)


def _short_reaction(reaction: str, width: int = 80) -> str:
    if len(reaction) <= width:
        return reaction
    return reaction[: width - 3] + "..."


def _dedupe_failures(failures: list[dict]) -> list[dict]:
    """Collapse the same failing row across cycles to its first failure."""
    grouped: dict[tuple, dict] = {}
    for failure in failures:
        key = (
            failure.get("reaction_id", ""),
            str(failure.get("mon1", "")),
            str(failure.get("mon2", "")),
            failure.get("failure", ""),
            failure.get("known_gap", ""),
        )
        current = grouped.get(key)
        if current is None or failure["cycle"] < current["cycle"]:
            grouped[key] = failure
    return sorted(
        grouped.values(),
        key=lambda item: (
            item.get("reaction_id", ""),
            item["cycle"],
            str(item.get("mon1", "")),
        ),
    )


def _write_report(summary: dict, output: Path) -> None:
    coverage = summary["coverage"]
    cycles = summary["config"]["cycles"]
    lines = [
        "# Extender benchmark",
        "",
        f"Generated: `{summary['generated_at']}`",
        "",
        f"Input: `{summary['config']['csv']}`",
        "",
        f"Cycles: `{', '.join(str(cycle) for cycle in cycles)}`",
        "",
        "A cycle is one addition/ROMP propagation step (cycle 1 is the "
        "terminated initiation product) or one step-growth/ring-opening "
        "monomer residue.",
        "",
        "## Coverage",
        "",
        _format_table(
            ["metric", "value"],
            [
                ["registry reactions", coverage["registry_reactions"]],
                ["configured reactions", coverage["configured_reactions"]],
                ["unconfigured reactions", coverage["unconfigured_reactions"]],
                ["reactions with data", coverage["reactions_with_data"]],
                ["configured with data", coverage["configured_with_data"]],
                [
                    "unconfigured with data",
                    len(coverage["unconfigured_with_data"]),
                ],
                [
                    "configured without data",
                    len(coverage["configured_without_data"]),
                ],
                [
                    "unconfigured without data",
                    len(coverage["unconfigured_without_data"]),
                ],
                ["rows total", coverage["rows_total"]],
                ["rows configured", coverage["rows_configured"]],
                ["rows unconfigured", coverage["rows_unconfigured"]],
            ],
        ),
        "",
        "### Rows by polymer class",
        "",
        _format_table(
            ["polymer class", "rows", "configured rows"],
            [
                [
                    polymer_class,
                    count,
                    coverage["configured_rows_by_class"].get(
                        polymer_class, 0
                    ),
                ]
                for polymer_class, count in sorted(
                    coverage["rows_by_class"].items()
                )
            ],
        ),
        "",
        "### Mechanisms without a recipe, or without data",
        "",
        _format_table(
            ["status", "polymer class", "rows", "reaction"],
            [
                [
                    "unconfigured",
                    coverage["class_by_reaction"].get(reaction, ""),
                    coverage["rows_by_reaction"].get(reaction, 0),
                    _short_reaction(reaction),
                ]
                for reaction in coverage["unconfigured_with_data"]
            ]
            + [
                [
                    "unconfigured, no data",
                    coverage["class_by_reaction"].get(reaction, ""),
                    0,
                    _short_reaction(reaction),
                ]
                for reaction in coverage["unconfigured_without_data"]
            ]
            + [
                [
                    "configured, no data",
                    coverage["class_by_reaction"].get(reaction, ""),
                    0,
                    _short_reaction(reaction),
                ]
                for reaction in coverage["configured_without_data"]
            ],
        ),
        "",
        "## Pass rate by cycle",
        "",
        _format_table(
            ["cycle", "attempted", "passed", "pass rate", "known gaps",
             "new failures", "p50 ms", "p95 ms", "median heavy atoms",
             "median delta heavy"],
            [
                [
                    cycle,
                    summary["cycles"][str(cycle)]["attempted"],
                    summary["cycles"][str(cycle)]["passed"],
                    _format_rate(
                        summary["cycles"][str(cycle)]["pass_rate"]
                    ),
                    summary["cycles"][str(cycle)]["known_gap_failures"],
                    summary["cycles"][str(cycle)]["new_failures"],
                    _format_number(
                        summary["cycles"][str(cycle)]["median_runtime_ms"]
                    ),
                    _format_number(
                        summary["cycles"][str(cycle)]["p95_runtime_ms"]
                    ),
                    _format_number(
                        summary["cycles"][str(cycle)]["median_heavy_atoms"],
                        1,
                    ),
                    _format_number(
                        summary["cycles"][str(cycle)]["median_delta_heavy"],
                        1,
                    ),
                ]
                for cycle in cycles
            ],
        ),
        "",
        "## Frontier size by cycle",
        "",
        _format_table(
            ["cycle", "mean products", "median products",
             "max products", "rows with >1 product"],
            [
                [
                    cycle,
                    _format_number(
                        summary["cycles"][str(cycle)]["mean_products"],
                        3,
                    ),
                    _format_number(
                        summary["cycles"][str(cycle)]["median_products"],
                        1,
                    ),
                    summary["cycles"][str(cycle)][
                        "max_products_observed"
                    ]
                    or 0,
                    summary["cycles"][str(cycle)][
                        "rows_with_multiple_products"
                    ],
                ]
                for cycle in cycles
            ],
        ),
        "",
        "## Failure taxonomy by cycle",
        "",
        _format_table(
            ["cycle", "failures"],
            [
                [
                    cycle,
                    _json_compact(
                        summary["cycles"][str(cycle)]["failures"]
                    ),
                ]
                for cycle in cycles
            ],
        ),
        "",
        "## First failure onset",
        "",
        _format_table(
            ["first failing cycle", "rows"],
            [
                [cycle, summary["failure_onset"].get(str(cycle), 0)]
                for cycle in cycles
            ]
            + [
                [
                    "never",
                    summary["rows_evaluated"]
                    - sum(summary["failure_onset"].values()),
                ]
            ],
        ),
        "",
        "## Row survival",
        "",
        _format_table(
            ["cycle", "rows still passing"],
            [
                [cycle, summary["survival"][str(cycle)]]
                for cycle in cycles
            ],
        ),
        "",
        f"Rows passing every requested cycle: "
        f"**{summary['rows_fully_ok']} / {summary['rows_evaluated']}**",
        "",
        "## Pass rate by polymer class and cycle",
        "",
    ]

    class_cycles = summary["by_class_cycle"]
    lines.append(
        _format_table(
            ["polymer class"]
            + [f"cycle {cycle}" for cycle in cycles],
            [
                [polymer_class]
                + [
                    _format_rate(
                        class_cycles[polymer_class][str(cycle)][
                            "pass_rate"
                        ]
                    )
                    if str(cycle) in class_cycles[polymer_class]
                    else "-"
                    for cycle in cycles
                ]
                for polymer_class in sorted(class_cycles)
            ],
        )
    )
    lines.append("")

    lines.append("## Per-reaction results")
    lines.append("")
    lines.append(
        _format_table(
            ["id", "polymer class", "type", "rows", "fully ok",
             "first failure", "couplings", "ring open", "reagents"],
            [
                [
                    info["id"],
                    info["polymer_class"],
                    info["type"],
                    info["rows"],
                    info["fully_ok_rows"],
                    _json_compact(info["first_failure_cycle"]) or "-",
                    ",".join(info["couplings"]) or "-",
                    "yes" if info["ring_open"] else "no",
                    ",".join(info["reagents"]) or "-",
                ]
                for _reaction, info in sorted(
                    summary["by_reaction"].items(),
                    key=lambda item: (
                        item[1]["polymer_class"],
                        item[0],
                    ),
                )
            ],
        )
    )
    lines.append("")

    lines.append("### Pass rate by reaction and cycle")
    lines.append("")
    lines.append(
        _format_table(
            ["id", "polymer class", "coupling", "rows"]
            + [f"c{cycle}" for cycle in cycles],
            [
                [
                    info["id"],
                    info["polymer_class"],
                    ",".join(info["couplings"]) or "-",
                    info["rows"],
                ]
                + [
                    _format_rate(
                        info["per_cycle"][str(cycle)]["pass_rate"]
                    )
                    for cycle in cycles
                ]
                for _reaction, info in sorted(
                    summary["by_reaction"].items(),
                    key=lambda item: (
                        item[1]["polymer_class"],
                        item[0],
                    ),
                )
            ],
        )
    )
    lines.append("")

    lines.append("### Reaction id map")
    lines.append("")
    lines.append(
        "`Rxx` ids are assigned by this benchmark, not by `rxn_map`: the "
        "configured reactions present in the run are sorted lexicographically "
        "by their reaction SMARTS and numbered from R01. The full SMARTS is "
        "also recorded in `records.csv` and `summary.json`."
    )
    lines.append("")
    lines.append(
        _format_table(
            ["id", "polymer class", "coupling", "ring open",
             "reaction SMARTS"],
            [
                [
                    info["id"],
                    info["polymer_class"],
                    ",".join(info["couplings"]) or "-",
                    "yes" if info["ring_open"] else "no",
                    f"`{_short_reaction(info.get('reaction_smarts', ''), 110)}`",
                ]
                for _reaction, info in sorted(
                    summary["by_reaction"].items(),
                    key=lambda item: item[1]["id"],
                )
            ],
        )
    )
    lines.append("")

    lines.append("## End-group provenance findings")
    lines.append("")
    allowed = summary.get("allowed_end_group_findings", [])
    violations = summary["new_end_group_findings"]
    lines.append(
        f"Violations: **{len(violations)}**; allowed new groups "
        f"(documented chemistry): **{len(allowed)}**."
    )
    lines.append("")
    if allowed:
        grouped = Counter(
            (
                finding["polymer_class"],
                finding["cycle"],
                tuple(finding["allowed_new_end_groups"]),
            )
            for finding in allowed
        )
        lines.append(
            _format_table(
                ["class", "cycle", "allowed new groups", "records"],
                [
                    [polymer_class, cycle, ",".join(groups), count]
                    for (polymer_class, cycle, groups), count in sorted(
                        grouped.items()
                    )
                ],
            )
        )
        lines.append("")

    lines.append("## Atom provenance findings")
    lines.append("")
    lines.append(
        f"Reactions with invented atoms: "
        f"**{len(summary['atom_provenance_findings'])}** records."
    )
    lines.append("")

    lines.append("## Failures")
    lines.append("")
    known = _dedupe_failures(
        [
            failure
            for failure in summary["failure_examples"]
            if failure["known_gap"]
        ]
    )
    unknown = _dedupe_failures(
        [
            failure
            for failure in summary["failure_examples"]
            if not failure["known_gap"]
        ]
    )
    lines.append(
        f"Known gaps: **{len(known)}**; new failures: "
        f"**{len(unknown)}**."
    )
    lines.append("")
    if unknown:
        lines.append(
            _format_table(
                ["id", "class", "first cycle", "failure", "mon1", "hint"],
                [
                    [
                        failure["reaction_id"],
                        failure["polymer_class"],
                        failure["cycle"],
                        failure["failure"],
                        str(failure["mon1"])[:50],
                        str(failure["hint"])[:80],
                    ]
                    for failure in unknown[:50]
                ],
            )
        )
        lines.append("")
    if known:
        lines.append("Known gaps:")
        lines.append("")
        lines.append(
            _format_table(
                ["id", "class", "first cycle", "mon1", "reason"],
                [
                    [
                        failure["reaction_id"],
                        failure["polymer_class"],
                        failure["cycle"],
                        str(failure["mon1"])[:50],
                        failure["known_gap"],
                    ]
                    for failure in known[:20]
                ],
            )
        )
        lines.append("")

    baseline = summary.get("baseline_diff")
    if baseline:
        lines.append("## Baseline comparison")
        lines.append("")
        lines.append(
            _format_table(
                ["metric", "count"],
                [
                    ["regressions", len(baseline["regressions"])],
                    ["fixes", len(baseline["fixes"])],
                    ["new rows", len(baseline["new_rows"])],
                    ["removed rows", len(baseline["removed_rows"])],
                    [
                        "unchanged failures",
                        len(baseline["unchanged_failures"]),
                    ],
                ],
            )
        )
        lines.append("")
        regression_cycles = baseline.get("regressions_by_cycle", {})
        fix_cycles = baseline.get("fixes_by_cycle", {})
        if regression_cycles or fix_cycles:
            lines.append(
                _format_table(
                    ["cycle", "regressions", "fixes"],
                    [
                        [
                            cycle,
                            regression_cycles.get(cycle, 0),
                            fix_cycles.get(cycle, 0),
                        ]
                        for cycle in sorted(
                            set(regression_cycles) | set(fix_cycles),
                            key=lambda value: int(value),
                        )
                    ],
                )
            )
            lines.append("")
        if baseline["regressions"]:
            lines.append("Regressed records:")
            lines.append("")
            for key in baseline["regressions"][:50]:
                lines.append(f"- `{key}`")
            lines.append("")

    (output / "report.md").write_text("\n".join(lines) + "\n")


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Benchmark the extender drivers over the full evaluation set."
        )
    )
    parser.add_argument(
        "--csv",
        type=Path,
        default=_CSV_PATH,
        help="evaluation CSV (default: tests/test_data/extender_test_data_all.csv)",
    )
    parser.add_argument(
        "--cycles",
        type=_parse_cycles,
        default=_DEFAULT_CYCLES,
        help="comma-separated cycles (default: 1,2,3,5,10)",
    )
    parser.add_argument(
        "--max-products",
        type=int,
        default=1,
        help="frontier cap passed to the driver (default: 1)",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="use at most N rows per reaction (default: all)",
    )
    parser.add_argument(
        "--class",
        dest="polymer_class",
        default=None,
        help="comma-separated polymer classes to include",
    )
    parser.add_argument(
        "--reaction-contains",
        default=None,
        help="only reactions whose SMARTS contain this substring",
    )
    parser.add_argument(
        "--no-atom-provenance",
        action="store_true",
        help="skip the isotope-labelled atom provenance run",
    )
    parser.add_argument(
        "--determinism",
        action="store_true",
        help="run every row twice and compare the products",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="result directory (default: benchmarks/results/<timestamp>)",
    )
    parser.add_argument(
        "--baseline",
        type=Path,
        default=None,
        help="previous summary.json (or result directory) to compare with",
    )
    parser.add_argument(
        "--fail-under",
        type=float,
        default=None,
        help=(
            "exit 1 if any requested cycle's pass rate is below "
            "this fraction"
        ),
    )
    parser.add_argument(
        "--fail-on-regression",
        action="store_true",
        help="exit 2 if the baseline comparison finds a regression",
    )
    parser.add_argument(
        "--fail-on-new-failures",
        action="store_true",
        help="exit 3 if any failure is not a known chemistry gap",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="only print the output directory",
    )
    parser.add_argument(
        "--mpi",
        action="store_true",
        help=(
            "require MPI (launch with `mpiexec -n N`, N > 1); "
            "MPI is used automatically when N > 1"
        ),
    )
    parser.add_argument(
        "--scaling-log",
        type=Path,
        default=None,
        help=(
            "append one row per cycle to this CSV "
            "(for max_products/concurrency plots)"
        ),
    )
    return parser.parse_args(argv)


def _init_mpi(args):
    """Return ``(MPI, comm, size, rank)``; fall back to sequential."""
    try:
        from mpi4py import MPI
    except Exception as exc:  # pragma: no cover - environment dependent
        if args.mpi:
            raise SystemExit(
                f"--mpi requested but mpi4py/MPI is unavailable: {exc}"
            )
        return None, None, 1, 0

    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
    if size <= 1 and args.mpi:
        raise SystemExit(
            "--mpi requires launching with `mpiexec -n N` where N > 1"
        )
    return MPI, comm, size, rank


def _load_rows(args) -> pd.DataFrame:
    """Read, filter, sort, and limit the evaluation CSV."""
    rows = pd.read_csv(args.csv)
    required = {
        "mon1",
        "mon2",
        "reaction_smarts",
        "polymer_class",
        "source_file",
        "row_index",
    }
    missing = required - set(rows.columns)
    if missing:
        raise SystemExit(f"CSV is missing columns: {sorted(missing)}")

    if args.polymer_class:
        classes = {
            value.strip()
            for value in args.polymer_class.split(",")
            if value.strip()
        }
        rows = rows[rows["polymer_class"].isin(classes)]
    if args.reaction_contains:
        rows = rows[
            rows["reaction_smarts"].str.contains(
                args.reaction_contains, regex=False
            )
        ]

    rows = rows.sort_values(
        ["reaction_smarts", "source_file", "row_index"]
    ).reset_index(drop=True)
    if args.limit is not None:
        rows = (
            rows.groupby("reaction_smarts", sort=True)
            .head(args.limit)
            .reset_index(drop=True)
        )
    return rows


def _row_payload(row) -> dict:
    """Convert a pandas row to a plain, picklable dict."""
    payload = {}
    for key, value in row.items():
        try:
            if bool(pd.isna(value)):
                payload[key] = None
                continue
        except (TypeError, ValueError):
            pass
        payload[key] = value
    if payload.get("row_index") is not None:
        try:
            payload["row_index"] = int(payload["row_index"])
        except (TypeError, ValueError):
            pass
    return payload


def _resolve_output(args) -> Path:
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S")
    output = args.output or (_DEFAULT_RESULTS / timestamp)
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    return output


def _reaction_context(reaction: str, config: dict, args) -> dict:
    """Compile everything a worker needs for one reaction."""
    context = {
        "config": config,
        "driver": None,
        "probe_driver": None,
        "reagents": config.get("params", {}).get("reagents", []),
        "allowed_labels": set(),
        "per_coupling_atoms": Counter(),
        "error": None,
    }
    try:
        context["driver"] = _build_driver(
            reaction, config, args.max_products
        )
    except Exception as exc:  # noqa: BLE001 - reported as config errors
        context["error"] = f"{type(exc).__name__}: {exc}"
        return context

    if config.get("type") in _STEP_LIKE_TYPES:
        try:
            context["probe_driver"] = _build_probe_driver(
                reaction, config, args.max_products
            )
            context["per_coupling_atoms"] = _allowed_new_atoms(
                context["probe_driver"]
            )
        except Exception:
            context["probe_driver"] = None
        context["allowed_labels"] = _allowed_new_end_groups(config)
    return context


def _error_records(
    reaction: str,
    config: dict,
    row,
    cycles: tuple[int, ...],
    failure: str,
    detail: str,
) -> list[dict]:
    """Build terminal records for a row that could not be processed."""
    length = lambda cycle: (  # noqa: E731 - local one-liner
        cycle - 1 if config.get("type") in _ADDITION_LIKE_TYPES else cycle
    )
    return [
        _record(
            reaction=reaction,
            config=config,
            row=row,
            cycle=cycle,
            length_arg=length(cycle),
            ok=False,
            failure=failure,
            detail=detail,
            hint="",
            warnings=[],
            products=[],
            heavy=None,
            rings=None,
            fragments=None,
            dummies=None,
            end_groups={},
            new_end_groups=set(),
            allowed_new_end_groups=set(),
            missing_end_groups=set(),
            unlabeled={},
            allowed_new_atoms={},
            excess_new_atoms={},
            delta_heavy=None,
            delta_rings=None,
            runtime_ms=0.0,
            deterministic=None,
            known_gap=_known_gap(_monomers(row)),
        )
        for cycle in cycles
    ]


def _process_row(
    reaction: str,
    config: dict,
    row,
    cycles: tuple[int, ...],
    args,
    context: dict,
) -> list[dict]:
    """Run every requested cycle for one candidate row."""
    if context.get("error"):
        return _error_records(
            reaction, config, row, cycles, FAIL_CONFIG, context["error"]
        )

    driver = context["driver"]
    probe_driver = context["probe_driver"]
    reagents = context["reagents"]
    allowed_labels = context["allowed_labels"]
    per_coupling_atoms = context["per_coupling_atoms"]
    is_step_growth = config.get("type") in _STEP_LIKE_TYPES

    monomers = _monomers(row)
    provenance_labels: set[str] = set()
    if is_step_growth:
        try:
            pool = _step_growth_pool(driver, monomers, reagents)
        except Exception:
            pool = []
        for pool_mol in pool:
            provenance_labels |= _end_group_labels(pool_mol)

    records = []
    previous = None
    for cycle in cycles:
        record = _run_row_cycle(
            reaction=reaction,
            config=config,
            row=row,
            cycle=cycle,
            driver=driver,
            probe_driver=probe_driver,
            reagents=reagents,
            provenance_labels=provenance_labels,
            allowed_labels=allowed_labels,
            per_coupling_atoms=per_coupling_atoms,
            max_products=args.max_products,
            atom_provenance=not args.no_atom_provenance,
            check_determinism=args.determinism,
        )
        if previous is not None and previous["ok"] and record["ok"]:
            record["delta_heavy"] = (
                record["heavy_atoms"] - previous["heavy_atoms"]
            )
            record["delta_rings"] = (
                record["rings"] - previous["rings"]
            )
            if record["delta_heavy"] <= 0:
                record["ok"] = False
                record["failure"] = FAIL_NO_GROWTH
                record["detail"] = (
                    "heavy-atom count did not grow: "
                    f"{previous['heavy_atoms']} -> "
                    f"{record['heavy_atoms']}"
                )
        if record["ok"]:
            previous = record
        records.append(record)
    return records


def _process_row_task(
    reaction: str, row: dict, args, cache: dict
) -> list[dict]:
    config = rxn_map.get(reaction)
    if config is None:
        return []
    if reaction not in cache:
        cache[reaction] = _reaction_context(reaction, config, args)
    return _process_row(
        reaction,
        config,
        row,
        tuple(args.cycles),
        args,
        cache[reaction],
    )


def _run_sequential(
    args, rows: pd.DataFrame, coverage: dict, output: Path
) -> int:
    cycles = tuple(args.cycles)
    cache: dict[str, dict] = {}
    records: list[dict] = []

    started = time.perf_counter()
    for reaction, group in rows.groupby("reaction_smarts", sort=True):
        config = rxn_map.get(reaction)
        if config is None:
            continue
        cache[reaction] = _reaction_context(reaction, config, args)
        for _, row in group.iterrows():
            records.extend(
                _process_row(
                    reaction,
                    config,
                    _row_payload(row),
                    cycles,
                    args,
                    cache[reaction],
                )
            )
    elapsed = time.perf_counter() - started

    return _finalize(
        args,
        rows,
        coverage,
        records,
        elapsed,
        output,
        mpi_info={"mode": "sequential", "n_ranks": 1, "n_workers": 1},
        worker_rss=[],
    )


def _mpi_worker(args, MPI, comm) -> None:
    """Worker loop: receive row tasks, process them, return records."""
    cache: dict[str, dict] = {}
    cycles = tuple(args.cycles)
    while True:
        task = comm.recv(source=0, tag=MPI.ANY_TAG)
        if task is None:
            break
        reaction, row = task
        try:
            records = _process_row_task(reaction, row, args, cache)
        except Exception as exc:  # noqa: BLE001 - keep the worker alive
            config = rxn_map.get(reaction)
            if config is None:
                records = []
            else:
                records = _error_records(
                    reaction,
                    config,
                    row,
                    cycles,
                    FAIL_EXCEPTION,
                    f"{type(exc).__name__}: {exc}",
                )
        comm.send(records, dest=0, tag=1)

    peak_rss_kb = (
        resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        if resource is not None
        else None
    )
    comm.gather(peak_rss_kb, root=0)


def _mpi_master_dispatch(
    MPI,
    comm,
    size: int,
    tasks: list[tuple],
    quiet: bool,
) -> tuple[list[dict], list[float]]:
    """Dynamically dispatch row tasks to worker ranks."""
    records: list[dict] = []
    active: set[int] = set()
    next_task = 0
    completed = 0

    def send_to(worker: int) -> bool:
        nonlocal next_task
        if next_task >= len(tasks):
            return False
        comm.send(tasks[next_task], dest=worker, tag=0)
        active.add(worker)
        next_task += 1
        return True

    for worker in range(1, size):
        send_to(worker)

    while active:
        status = MPI.Status()
        result = comm.recv(
            source=MPI.ANY_SOURCE, tag=1, status=status
        )
        source = status.Get_source()
        records.extend(result)
        completed += 1
        if not quiet and completed % 100 == 0:
            print(
                f"  completed {completed}/{len(tasks)} rows",
                flush=True,
            )
        if not send_to(source):
            active.remove(source)

    for worker in range(1, size):
        comm.send(None, dest=worker, tag=0)
    gathered = comm.gather(None, root=0)
    worker_rss = [value for value in gathered if value is not None]
    return records, worker_rss


def _run_mpi_master(
    args,
    MPI,
    comm,
    size: int,
    rows: pd.DataFrame,
    coverage: dict,
    output: Path,
) -> int:
    tasks: list[tuple] = []
    for reaction, group in rows.groupby("reaction_smarts", sort=True):
        if rxn_map.get(reaction) is None:
            continue
        for _, row in group.iterrows():
            tasks.append((reaction, _row_payload(row)))

    started = time.perf_counter()
    records, worker_rss = _mpi_master_dispatch(
        MPI, comm, size, tasks, args.quiet
    )
    elapsed = time.perf_counter() - started

    return _finalize(
        args,
        rows,
        coverage,
        records,
        elapsed,
        output,
        mpi_info={
            "mode": "mpi",
            "n_ranks": size,
            "n_workers": max(1, size - 1),
        },
        worker_rss=worker_rss,
    )


def _append_scaling_log(
    args, summary: dict, output: Path
) -> None:
    """Append one long-format row per cycle for later plotting."""
    if not args.scaling_log:
        return
    path = Path(args.scaling_log)
    path.parent.mkdir(parents=True, exist_ok=True)

    mpi = summary["config"]["mpi"]
    peak_mb = (
        round(summary["peak_rss_kb"] / 1024, 1)
        if summary.get("peak_rss_kb")
        else None
    )
    worker_peak = summary.get("worker_peak_rss_kb")
    worker_peak_mb = (
        round(worker_peak / 1024, 1) if worker_peak else None
    )
    rows = []
    for cycle in summary["config"]["cycles"]:
        info = summary["cycles"][str(cycle)]
        rows.append(
            {
                "timestamp": summary["generated_at"],
                "run_id": output.name,
                "mode": mpi["mode"],
                "n_ranks": mpi["n_ranks"],
                "n_workers": mpi["n_workers"],
                "max_products": args.max_products,
                "cycles": ",".join(
                    str(value)
                    for value in summary["config"]["cycles"]
                ),
                "cycle": cycle,
                "attempted": info["attempted"],
                "passed": info["passed"],
                "pass_rate": info["pass_rate"],
                "p50_ms": info["median_runtime_ms"],
                "p95_ms": info["p95_runtime_ms"],
                "mean_products": info["mean_products"],
                "wall_s": summary["runtime_s"],
                "peak_rss_mb": peak_mb,
                "worker_peak_rss_mb": worker_peak_mb,
            }
        )

    frame = pd.DataFrame(rows)
    frame.to_csv(
        path,
        mode="a",
        header=not path.exists(),
        index=False,
    )


def _finalize(
    args,
    rows: pd.DataFrame,
    coverage: dict,
    records: list[dict],
    elapsed: float,
    output: Path,
    mpi_info: dict,
    worker_rss: list[float],
) -> int:
    """Aggregate, write outputs, print the summary, and apply gates."""
    cycles = tuple(args.cycles)
    records.sort(
        key=lambda record: (
            record["reaction_smarts"],
            str(record["source_file"]),
            str(record["row_index"]),
            record["cycle"],
        )
    )

    summary = _summarise(records, cycles)
    summary["coverage"] = coverage
    summary["config"] = {
        "csv": str(args.csv),
        "cycles": list(cycles),
        "max_products": args.max_products,
        "limit": args.limit,
        "polymer_class": args.polymer_class,
        "reaction_contains": args.reaction_contains,
        "atom_provenance": not args.no_atom_provenance,
        "determinism": args.determinism,
        "mpi": mpi_info,
    }
    summary["generated_at"] = datetime.now(timezone.utc).isoformat()
    summary["runtime_s"] = round(elapsed, 3)
    summary["peak_rss_kb"] = (
        resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        if resource is not None
        else None
    )
    summary["worker_peak_rss_kb"] = (
        max(worker_rss) if worker_rss else None
    )
    summary["record_index"] = {
        _record_key(record): {
            "ok": record["ok"],
            "failure": record["failure"],
        }
        for record in records
    }
    summary["known_gaps"] = KNOWN_GAPS

    if args.baseline:
        summary["baseline_diff"] = _compare_baseline(
            records, args.baseline
        )

    frame = pd.DataFrame(records)
    frame.to_csv(output / "records.csv", index=False)

    (output / "summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False, default=str)
        + "\n"
    )
    _write_report(summary, output)
    _append_scaling_log(args, summary, output)

    if not args.quiet:
        print(f"Input:   {args.csv}")
        print(
            f"Mode:    {mpi_info['mode']} "
            f"(ranks={mpi_info['n_ranks']}, "
            f"workers={mpi_info['n_workers']})"
        )
        print(
            "Coverage: "
            f"{coverage['rows_total']} rows / "
            f"{coverage['reactions_with_data']} reactions "
            f"({coverage['rows_configured']} configured rows, "
            f"{coverage['rows_unconfigured']} unconfigured rows)"
        )
        print(
            "Registry: "
            f"{coverage['configured_reactions']} configured, "
            f"{coverage['unconfigured_reactions']} unconfigured; "
            f"{len(coverage['unconfigured_with_data'])} unconfigured with data, "
            f"{len(coverage['configured_without_data'])} configured without data"
        )
        header = (
            f"{'cycle':>5} {'attempted':>9} {'passed':>7} "
            f"{'pass_rate':>9} {'known':>5} {'new':>4} "
            f"{'invented_end':>12} {'invented_atom':>13} "
            f"{'mean_prod':>9} {'p50_ms':>8}"
        )
        print(header)
        print("-" * len(header))
        for cycle in cycles:
            info = summary["cycles"][str(cycle)]
            pass_rate = (
                "-"
                if info["pass_rate"] is None
                else f"{100 * info['pass_rate']:.1f}%"
            )
            print(
                f"{cycle:>5} {info['attempted']:>9} {info['passed']:>7} "
                f"{pass_rate:>9} {info['known_gap_failures']:>5} "
                f"{info['new_failures']:>4} "
                f"{info['invented_end_groups']:>12} "
                f"{info['invented_atoms']:>13} "
                f"{(info['mean_products'] or 0):>9.3f} "
                f"{(info['median_runtime_ms'] or 0):>8.2f}"
            )
        new_failures = sum(
            info["new_failures"] for info in summary["cycles"].values()
        )
        print(
            f"Rows passing every cycle: "
            f"{summary['rows_fully_ok']}/{summary['rows_evaluated']}"
        )
        print(f"New (non-known-gap) failure records: {new_failures}")
        if summary.get("baseline_diff"):
            diff = summary["baseline_diff"]
            print(
                f"Baseline: {len(diff['regressions'])} regressions, "
                f"{len(diff['fixes'])} fixes"
            )
        print(f"Runtime: {summary['runtime_s']} s")
        if summary.get("peak_rss_kb"):
            print(f"Peak RSS (master): {summary['peak_rss_kb'] / 1024:.1f} MiB")
        if summary.get("worker_peak_rss_kb"):
            print(
                "Peak RSS (worker): "
                f"{summary['worker_peak_rss_kb'] / 1024:.1f} MiB"
            )
    print(f"Results: {output}")

    if args.fail_on_new_failures and any(
        record["failure"] != OK and not record["known_gap"]
        for record in records
    ):
        return 3
    if (
        args.fail_on_regression
        and summary.get("baseline_diff")
        and summary["baseline_diff"]["regressions"]
    ):
        return 2
    if args.fail_under is not None:
        rates = [
            info["pass_rate"]
            for info in summary["cycles"].values()
            if info["pass_rate"] is not None
        ]
        if rates and min(rates) < args.fail_under:
            return 1
    return 0


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)

    if not args.csv.exists():
        raise SystemExit(
            f"Evaluation CSV not found: {args.csv}\n"
            "Run `python tests/extract_extender_test_data.py` first."
        )
    if args.max_products < 1:
        raise SystemExit("--max-products must be at least 1")

    MPI, comm, size, rank = _init_mpi(args)

    # Workers do not load the CSV; they only process tasks from the master.
    if size > 1 and rank != 0:
        _mpi_worker(args, MPI, comm)
        return 0

    if size > 1:
        try:
            rows = _load_rows(args)
        except BaseException:
            # Do not leave workers blocked on a master that failed to load.
            for worker in range(1, size):
                comm.send(None, dest=worker, tag=0)
            comm.gather(None, root=0)
            raise
    else:
        rows = _load_rows(args)

    coverage = _coverage(rows)
    output = _resolve_output(args)

    if size > 1:
        return _run_mpi_master(
            args, MPI, comm, size, rows, coverage, output
        )
    return _run_sequential(args, rows, coverage, output)


if __name__ == "__main__":
    raise SystemExit(main())
