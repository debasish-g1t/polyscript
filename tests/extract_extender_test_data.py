"""Build the full extender evaluation set from ``test_data/sample_outputs``.

For every reaction SMARTS key in :mod:`polyscript.extender.map.rxn_map`, this
script selects **all** unique candidate rows whose reaction SMARTS matches an
entry in ``rxn_map``, whose non-empty monomer count equals the reaction's
reactant-template count, and whose monomers match the reaction's reactant
components.

The selected rows are merged into a single
``test_data/extender_test_data_all.csv`` file.

Usage::

    python tests/extract_extender_test_data.py
"""

from __future__ import annotations

import sys
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

_PROJECT = Path(__file__).resolve().parent.parent
if str(_PROJECT) not in sys.path:
    sys.path.insert(0, str(_PROJECT))

from polyscript.extender.map import rxn_map  # noqa: E402

RDLogger.DisableLog("rdApp.*")

_TEST_DATA = _PROJECT / "tests" / "test_data"
_SAMPLE_DIR = _TEST_DATA / "sample_outputs"
_OUT_FILE = _TEST_DATA / "extender_test_data_all.csv"
_LEGACY_OUT_DIR = _TEST_DATA / "extender_test_data"

_POLYM_COLUMNS = ("mon1", "mon2", "polym", "polymer_class")
_CSV_COLUMNS = (
    "mon1",
    "mon2",
    "reaction_smarts",
    "polym",
    "polymer_class",
    "source_file",
    "row_index",
)


def iter_reaction_smarts(reactions) -> list[str]:
    """Return reaction SMARTS from the registry or a collection."""
    if isinstance(reactions, str):
        return [reactions]
    if isinstance(reactions, dict):
        return [key for key in reactions if isinstance(key, str)]
    if isinstance(reactions, (set, list, tuple)):
        out = []
        for value in reactions:
            out.extend(iter_reaction_smarts(value))
        return out
    raise TypeError(f"Unsupported rxn_map value: {type(reactions)!r}")


def extract_reaction(polym: str) -> str:
    """Return the reaction SMARTS embedded in a polymerizer ``polym`` string."""
    parts = str(polym).split("=>")
    if len(parts) < 3:
        raise ValueError(f"Unexpected polym format: {polym!r}")
    return "=>".join(parts[1:-1])


def split_reactant_components(reaction: str) -> list[str]:
    """Split the reactant side of a reaction SMARTS into valid components.

    Top-level dots and grouping parentheses are split, while ``[]`` bracket
    expressions (including recursive ``$(...)`` SMARTS) are preserved.
    """
    left = reaction.split(">>", 1)[0]

    def _split(text: str) -> list[str]:
        parts: list[str] = []
        buf: list[str] = []
        bracket_depth = 0
        paren_depth = 0
        for char in text:
            if char == "[":
                bracket_depth += 1
            elif char == "]":
                bracket_depth = max(0, bracket_depth - 1)
            elif bracket_depth == 0 and char == "(":
                paren_depth += 1
            elif bracket_depth == 0 and char == ")":
                paren_depth = max(0, paren_depth - 1)

            if char == "." and bracket_depth == 0 and paren_depth == 0:
                parts.append("".join(buf))
                buf = []
            else:
                buf.append(char)
        parts.append("".join(buf))

        result: list[str] = []
        for part in parts:
            part = part.strip()
            if part.startswith("(") and part.endswith(")"):
                inner = part[1:-1]
                depth = 0
                balanced = True
                for char in inner:
                    if char == "(":
                        depth += 1
                    elif char == ")":
                        depth -= 1
                        if depth < 0:
                            balanced = False
                            break
                if balanced and depth == 0:
                    result.extend(_split(inner))
                    continue
            if part:
                result.append(part)
        return result

    return _split(left)


def compile_components(reaction: str) -> list[Chem.Mol]:
    """Compile the reactant components of *reaction* into SMARTS queries."""
    queries = []
    for component in split_reactant_components(reaction):
        query = Chem.MolFromSmarts(component)
        if query is not None:
            queries.append(query)
    return queries


def reactant_template_count(reaction: str) -> int | None:
    """Return the number of reactant templates in *reaction*."""
    try:
        reaction_obj = AllChem.ReactionFromSmarts(reaction)
    except Exception:
        return None
    if reaction_obj is None:
        return None
    return reaction_obj.GetNumReactantTemplates()


def _is_empty(value: object) -> bool:
    """True for ``None``, NaN/NA, and the string forms of an empty value."""
    if value is None:
        return True
    try:
        if bool(pd.isna(value)):
            return True
    except (TypeError, ValueError):
        pass
    return str(value).strip() in {"", "None", "nan", "NaN"}


def candidate_matches(
    smiles: object, components: list[str] | list[Chem.Mol]
) -> bool:
    """True if *smiles* matches at least one reactant component."""
    if _is_empty(smiles):
        return False
    molecule = Chem.MolFromSmiles(str(smiles).strip())
    if molecule is None:
        return False
    for component in components:
        query = (
            component
            if isinstance(component, Chem.Mol)
            else Chem.MolFromSmarts(component)
        )
        if query is not None and molecule.HasSubstructMatch(query):
            return True
    return False


def load_samples() -> list[dict]:
    """Load every sample row with its source file and reaction SMARTS."""
    rows: list[dict] = []
    for source_file in sorted(_SAMPLE_DIR.rglob("*.parquet")):
        dataframe = pd.read_parquet(source_file)
        missing = [c for c in _POLYM_COLUMNS if c not in dataframe.columns]
        if missing:
            raise ValueError(f"{source_file} is missing columns: {missing}")
        for row_index, row in dataframe.iterrows():
            rows.append(
                {
                    "source_file": source_file,
                    "row_index": int(row_index),
                    "reaction": extract_reaction(row["polym"]),
                    "polymer_class": str(row["polymer_class"]),
                    "mon1": row["mon1"],
                    "mon2": row["mon2"],
                    "row": row,
                }
            )
    return rows


def unique_candidates(
    rows: list[dict],
    components: list[Chem.Mol],
    num_templates: int | None,
) -> list[dict]:
    """Return every unique, shape-compatible row for one reaction."""
    candidates: list[dict] = []
    seen: set[tuple[str, str]] = set()

    for row in rows:
        monomers = [
            monomer
            for monomer in (row["mon1"], row["mon2"])
            if not _is_empty(monomer)
        ]
        if num_templates is not None and len(monomers) != num_templates:
            continue
        if not monomers:
            continue
        if not all(
            candidate_matches(monomer, components) for monomer in monomers
        ):
            continue

        key = (str(row["mon1"]), str(row["mon2"]))
        if key in seen:
            continue
        seen.add(key)
        candidates.append(row)
    return candidates


def select_samples(rows: list[dict], known_reactions: set[str]) -> dict:
    """Select all unique candidates for every matching reaction."""
    reaction_rows: dict[str, list[dict]] = defaultdict(list)
    for row in rows:
        if row["reaction"] in known_reactions:
            reaction_rows[row["reaction"]].append(row)

    selected: list[dict] = []
    reaction_manifest: list[dict] = []
    for reaction in sorted(reaction_rows):
        candidates = unique_candidates(
            reaction_rows[reaction],
            compile_components(reaction),
            reactant_template_count(reaction),
        )
        if not candidates:
            continue
        selected.extend(candidates)
        reaction_manifest.append(
            {
                "smarts": reaction,
                "polymer_class": candidates[0]["polymer_class"],
                "sample_count": len(candidates),
                "samples": [
                    {
                        "source_file": str(
                            row["source_file"].relative_to(_SAMPLE_DIR)
                        ),
                        "row_index": row["row_index"],
                        "mon1": str(row["mon1"]),
                        "mon2": str(row["mon2"]),
                    }
                    for row in candidates
                ],
            }
        )

    return {
        "selected": selected,
        "matched_reactions": reaction_manifest,
        "uncovered_reactions": sorted(known_reactions - set(reaction_rows)),
    }


def write_dataset(
    selected: list[dict], output_file: Path = _OUT_FILE
) -> None:
    """Merge selected rows into the full evaluation CSV."""
    output_file = Path(output_file)

    if _LEGACY_OUT_DIR.exists():
        for path in _LEGACY_OUT_DIR.rglob("*"):
            if path.is_file():
                path.unlink()
        for path in sorted(_LEGACY_OUT_DIR.rglob("*"), reverse=True):
            if path.is_dir():
                path.rmdir()
        _LEGACY_OUT_DIR.rmdir()

    records = []
    for selected_row in selected:
        row = selected_row["row"]
        records.append(
            {
                "mon1": row["mon1"],
                "mon2": row["mon2"],
                "reaction_smarts": selected_row["reaction"],
                "polym": row["polym"],
                "polymer_class": row["polymer_class"],
                "source_file": str(
                    selected_row["source_file"].relative_to(_SAMPLE_DIR)
                ),
                "row_index": selected_row["row_index"],
            }
        )

    output_file.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(records, columns=_CSV_COLUMNS).to_csv(
        output_file, index=False
    )


def print_summary(result: dict, output_file: Path) -> None:
    """Print selection counts for the full set."""
    print(f"Matched reactions: {len(result['matched_reactions'])}")
    print(f"Selected samples:  {len(result['selected'])}")
    print(f"Output:            {output_file}")

    class_counts = Counter(
        row["polymer_class"] for row in result["selected"]
    )
    for polymer_class in sorted(class_counts):
        print(f"  {polymer_class:18s} {class_counts[polymer_class]:4d}")

    counts = [item["sample_count"] for item in result["matched_reactions"]]
    if counts:
        print(f"Per reaction:      min={min(counts)}, max={max(counts)}")


def main() -> None:
    known_reactions = set(iter_reaction_smarts(rxn_map))
    rows = load_samples()
    result = select_samples(rows, known_reactions)

    write_dataset(result["selected"])
    print_summary(result, _OUT_FILE)

    if result["uncovered_reactions"]:
        print("\nUncovered reactions (no matching samples):")
        for reaction in result["uncovered_reactions"]:
            print(f"  {reaction}")


if __name__ == "__main__":
    main()
