"""Keep mechanism_table.tex, rxn_map and the practical example page in sync."""

import importlib.util
import re
from pathlib import Path

from rdkit import Chem

from polyscript.extender.map import rxn_map

_ROOT = Path(__file__).resolve().parent.parent
_TABLE = _ROOT / "benchmarks" / "mechanism_table.tex"
_EXAMPLES = _ROOT / "benchmarks" / "mechanism_examples.md"
_SCRIPT = _ROOT / "benchmarks" / "extract_extension_traces.py"

_EXPECTED = (
    {f"A{i}" for i in range(1, 5)}
    | {f"C{i}" for i in range(1, 9)}
    | {f"R{i}" for i in range(1, 6)}
)

_STATUSES = {"commercial", "known", "known/specialty", "needs candidate"}

# Linkage that the generated product must contain for each mechanism.
_LINKAGE = {
    "C1": "[CX3](=O)[OX2]",
    "C2": "[CX3](=O)[NX3]",
    "C3": "[NX3][CX3](=O)[OX2]",
    "C4": "[#6](=O)[#7][#6](=O)",
    "C5": "[NX3]1[CX3](=O)[OX2][CX4][CX4]1",
    "C6": "[OX2][CX3](=O)[OX2]",
    "C7": "[c][OX2][c]",
    "C8": "[c][OX2][c]",
    "R1": "[CX3](=O)[OX2]",
    "R2": "[CX3](=O)[NX3]",
    "R3": "[CX3](=O)[OX2]",
    "R4": "[OX2][CX4]",
}


def _table_ids() -> set[str]:
    text = _TABLE.read_text()
    return set(
        re.findall(r"^([ACR]\d+) & ", text, flags=re.MULTILINE)
    )


def _assigned_keys() -> dict[str, list[str]]:
    assigned: dict[str, list[str]] = {}
    for key, config in rxn_map.items():
        if config is None:
            continue
        mechanism = config.get("mechanism")
        assert mechanism in _EXPECTED, (mechanism, key)
        assigned.setdefault(mechanism, []).append(key)
    return assigned


def _module():
    spec = importlib.util.spec_from_file_location(
        "extract_extension_traces", _SCRIPT
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _mechanism_specs() -> list[dict]:
    return _module()._MECHANISMS


def test_table_and_registry_cover_the_same_mechanisms():
    assert _table_ids() == _EXPECTED
    assert set(_assigned_keys()) == _EXPECTED


def test_every_registry_entry_has_a_mechanism_code():
    for key, config in rxn_map.items():
        if config is None:
            continue
        assert config.get("mechanism") in _EXPECTED, key


def test_illustration_page_covers_every_mechanism():
    text = _EXAMPLES.read_text()
    sections = set(
        re.findall(r"^## ([ACR]\d+) — ", text, flags=re.MULTILINE)
    )
    assert sections == _EXPECTED
    assert "dataset" not in text.lower()

    # Every full registry key assigned to a mechanism is recorded verbatim.
    for keys in _assigned_keys().values():
        for key in keys:
            assert key in text

    # Every practical example has a representative and a trace panel.
    links = re.findall(r"!\[[^\]]*\]\(([^)]+)\)", text)
    assert links
    for link in links:
        assert (_ROOT / "benchmarks" / link).is_file(), link
    pending = {
        spec["id"]
        for spec in _mechanism_specs()
        if spec.get("status") == "needs candidate"
    }
    for spec in _mechanism_specs():
        if spec["id"] in pending:
            continue
        prefix = spec.get("prefix", spec["id"])
        assert f"traces/mechanisms/practical/{prefix}.svg" in links
        assert f"traces/mechanisms/practical/{prefix}_trace.svg" in links


def test_every_mechanism_has_a_practical_example_or_a_placeholder():
    module = _module()
    text = _EXAMPLES.read_text()
    for spec in _mechanism_specs():
        assert spec["status"] in _STATUSES, spec["id"]
        if spec["status"] == "needs candidate":
            assert not spec.get("monomers"), spec["id"]
            assert f"## {spec['id']} " in text
            assert "Candidate needed." in text
            continue

        assert spec.get("monomers"), spec["id"]
        selection = module._mechanism_selection(spec)
        result = module._trace(selection["kind"], selection)
        assert result["final"], spec["id"]
        product = result["final"][0]
        assert "*" not in product
        molecule = Chem.MolFromSmiles(product)
        assert molecule is not None, spec["id"]

        linkage = _LINKAGE.get(spec["id"])
        if linkage:
            assert molecule.HasSubstructMatch(
                Chem.MolFromSmarts(linkage)
            ), (spec["id"], linkage, product)


def _svg_edges(path: Path) -> set[frozenset[int]]:
    text = path.read_text()
    return {
        frozenset((int(begin), int(end)))
        for begin, end in re.findall(
            r"class='bond-\d+ atom-(\d+) atom-(\d+)'", text
        )
    }


def _mol_edges(mol) -> set[frozenset[int]]:
    return {
        frozenset((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
        for bond in mol.GetBonds()
    }


def test_mechanism_svg_images_match_generated_structures():
    """Catch stale or mismatched SVGs by comparing exact bond graphs."""
    module = _module()
    traces = module._MECH_TRACES
    for spec in module._MECHANISMS:
        if spec.get("status") == "needs candidate":
            continue
        selection = module._mechanism_selection(spec)
        result = module._trace(selection["kind"], selection)
        prefix = selection["prefix"]

        if selection["kind"] in {"addition", "romp"}:
            for index, smiles in enumerate(result["intermediates"]):
                path = traces / f"{prefix}_{index:02d}.svg"
                expected = Chem.MolFromSmiles(smiles[0])
                assert _svg_edges(path) == _mol_edges(expected), (
                    spec["id"],
                    path.name,
                )
            final = Chem.MolFromSmiles(result["final"][0])
            assert _svg_edges(
                traces / f"{prefix}_final.svg"
            ) == _mol_edges(final)
            assert _svg_edges(traces / f"{prefix}.svg") == _mol_edges(final)
        else:
            for index, smiles_list in enumerate(result["intermediates"]):
                for offset, smile in enumerate(smiles_list):
                    suffix = (
                        ""
                        if len(smiles_list) == 1
                        else chr(ord("a") + offset)
                    )
                    path = (
                        traces
                        / f"{prefix}_{index + 1:02d}{suffix}.svg"
                    )
                    expected = Chem.MolFromSmiles(smile)
                    assert _svg_edges(path) == _mol_edges(expected), (
                        spec["id"],
                        path.name,
                    )
            final = Chem.MolFromSmiles(result["final"][0])
            assert _svg_edges(traces / f"{prefix}.svg") == _mol_edges(final)
