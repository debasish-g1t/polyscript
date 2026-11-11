"""Tests for the BRICS sub-package: decomposition, building, and batch utility."""

import os

import pandas as pd
import pytest
from rdkit import Chem

from polyscript.brics import BRICSBuild, BRICSDecompose, apply_BRICS


# ── Core BRICS functions ─────────────────────────────────────────────────────


def test_BRICSDecompose():
    """Decomposition returns the expected fragment sets."""
    src_mols = ["CCC(=O)CC(Cl)", "C1CCCCC1CC(N)CCO"]
    src_mols = [Chem.MolFromSmiles(sm) for sm in src_mols]

    res_mols = [BRICSDecompose(mol) for mol in src_mols]
    expected = [{"CCC(=O)CCCl"}, {"[15*]C1CCCCC1", "[8*]CC(N)CCO"}]
    assert res_mols == expected, f"Unexpected decomposition: {res_mols}"


def test_BRICSBuild():
    """Recombination produces a product with the expected SMILES and recipe fields."""
    src_mols = ["CCC(=O)CC(Cl)", "C1CCCCC1CC(N)CCO"]
    src_mols = [Chem.MolFromSmiles(sm) for sm in src_mols]

    res_mols = [BRICSDecompose(mol) for mol in src_mols]
    s_res_mols = [i for x in res_mols for i in x]

    build_mols = BRICSBuild(
        [Chem.MolFromSmiles(sm) for sm in s_res_mols], scrambleReagents=False
    )

    first_result = list(build_mols)[0]
    text, mol = first_result

    # The product SMILES should be the recombined molecule
    assert isinstance(mol, Chem.rdchem.Mol), "Second tuple element must be an RDKit Mol"
    assert Chem.MolToSmiles(mol) == "NC(CCO)CC1CCCCC1", (
        f"Unexpected product: {Chem.MolToSmiles(mol)}"
    )

    # The recipe text encodes fragment, seed, reaction — order may vary
    for marker in ("|fragment|", "|seed|", "|rxn_smiles|", "|rxn_smarts|"):
        assert marker in text, f"Missing recipe marker {marker!r}"

    assert "[15*]C1CCCCC1" in text, "Recipe should mention cyclohexyl fragment"
    assert "[8*]CC(N)CCO" in text, "Recipe should mention amino-alcohol fragment"


# ── apply_BRICS: return-value path (no output_dir) ───────────────────────────


SRC_SMILES = ["CCC(=O)CC(Cl)", "C1CCCCC1CC(N)CCO"]


def test_apply_BRICS_return_value():
    """When output_dir is None, apply_BRICS returns a list of result dicts."""
    results = apply_BRICS(SRC_SMILES, log_level=0)

    assert isinstance(results, list), "Should return a list"
    assert len(results) == 1, f"Expected 1 build result, got {len(results)}"

    r = results[0]
    assert isinstance(r, dict), "Each result must be a dict"
    assert "product_smiles" in r
    assert "recipe" in r
    assert "mol" in r

    # Product should be the recombined molecule
    assert r["product_smiles"] == "NC(CCO)CC1CCCCC1"

    # Recipe should contain both fragments (order-independent)
    assert "[15*]C1CCCCC1" in r["recipe"]
    assert "[8*]CC(N)CCO" in r["recipe"]

    # mol should be an RDKit Mol
    assert isinstance(r["mol"], Chem.rdchem.Mol)


def test_apply_BRICS_return_no_files_polluted(tmp_path):
    """When output_dir is None, no files are written to CWD."""
    cwd_before = set(os.listdir(os.getcwd()))
    apply_BRICS(SRC_SMILES, log_level=0)
    cwd_after = set(os.listdir(os.getcwd()))

    new_files = cwd_after - cwd_before
    offending = {
        f for f in new_files
        if f in ("decomposed_unique.txt", "brics_results.csv", "brics_results.parquet")
    }
    assert not offending, f"apply_BRICS wrote files to CWD: {offending}"


# ── apply_BRICS: is_decomposed flag ──────────────────────────────────────────

# These are the fragments that BRICSDecompose would produce from SRC_SMILES
PRE_DECOMPOSED = ["[15*]C1CCCCC1", "[8*]CC(N)CCO", "CCC(=O)CCCl"]


def test_apply_BRICS_is_decomposed_true():
    """When is_decomposed=True, pre-decomposed fragments go straight to build."""
    # The decomposed path (default) and is_decomposed path should produce
    # the same build results from the same fragments.
    results_fresh = apply_BRICS(SRC_SMILES, log_level=0)
    results_pre = apply_BRICS(PRE_DECOMPOSED, is_decomposed=True, log_level=0)

    assert len(results_fresh) == len(results_pre), (
        f"is_decomposed=True should match fresh decomposition: "
        f"{len(results_fresh)} vs {len(results_pre)}"
    )

    # Product SMILES should be identical
    fresh_products = [r["product_smiles"] for r in results_fresh]
    pre_products = [r["product_smiles"] for r in results_pre]
    assert sorted(fresh_products) == sorted(pre_products), (
        f"Mismatched products:\n  fresh={sorted(fresh_products)}\n  pre=  {sorted(pre_products)}"
    )


def test_apply_BRICS_is_decomposed_skips_invalid():
    """is_decomposed=True with non-BRICS SMILES still tries to build."""
    # Passing a non-fragment SMILES with is_decomposed=True won't crash,
    # but it also won't produce build results since it has no wildcards.
    results = apply_BRICS(["CCO"], is_decomposed=True, log_level=0)
    assert isinstance(results, list), "Should return a list, not crash"
    # CCO has no attachment points, so BRICSBuild won't match anything


def test_apply_BRICS_is_decomposed_file_output(tmp_path):
    """is_decomposed=True with output_dir writes the correct files."""
    apply_BRICS(
        PRE_DECOMPOSED, is_decomposed=True,
        output_dir=str(tmp_path), log_level=0,
    )

    # decomposed_unique.txt should contain the input fragments (deduped)
    content = (tmp_path / "decomposed_unique.txt").read_text()
    lines = set(content.strip().splitlines())
    expected = set(PRE_DECOMPOSED)
    assert lines == expected, f"Fragment file mismatch: {lines}"

    # CSV should exist with the same product
    df = pd.read_csv(str(tmp_path / "brics_results.csv"))
    assert len(df) == 1
    assert df.iloc[0]["product_smiles"] == "NC(CCO)CC1CCCCC1"


def test_apply_BRICS_is_decomposed_false_is_default():
    """The default (is_decomposed=False) behaves identically to before."""
    results_default = apply_BRICS(SRC_SMILES, log_level=0)
    results_explicit = apply_BRICS(SRC_SMILES, is_decomposed=False, log_level=0)
    assert len(results_default) == len(results_explicit)
    assert results_default[0]["product_smiles"] == results_explicit[0]["product_smiles"]


# ── apply_BRICS: file-output path (with output_dir) ──────────────────────────


def _run_to_dir(tmp_path, fmt="csv"):
    """Run apply_BRICS with output_dir=str(tmp_path), return file paths."""
    apply_BRICS(SRC_SMILES, output_dir=str(tmp_path), fmt=fmt, log_level=0)
    table_name = "brics_results.csv" if fmt == "csv" else "brics_results.parquet"
    return (
        tmp_path / "decomposed_unique.txt",
        tmp_path / table_name,
    )


def test_output_files_exist(tmp_path):
    """apply_BRICS writes both expected files to output_dir."""
    deco_path, table_path = _run_to_dir(tmp_path)
    assert deco_path.exists(), "Missing decomposed_unique.txt"
    assert table_path.exists(), "Missing brics_results.csv"


def test_output_files_not_in_cwd(tmp_path):
    """apply_BRICS writes to output_dir, NOT the current working directory."""
    cwd_before = set(os.listdir(os.getcwd()))
    _run_to_dir(tmp_path)
    cwd_after = set(os.listdir(os.getcwd()))
    assert cwd_before == cwd_after, "apply_BRICS polluted the original CWD"


def test_decomposed_unique_content(tmp_path):
    """decomposed_unique.txt contains the correct unique fragments."""
    deco_path, _ = _run_to_dir(tmp_path)
    lines = set(deco_path.read_text().strip().splitlines())
    expected = {"[15*]C1CCCCC1", "[8*]CC(N)CCO", "CCC(=O)CCCl"}
    assert lines == expected, f"Fragment mismatch: {lines}"


def test_results_csv_content(tmp_path):
    """brics_results.csv has product_smiles and recipe columns."""
    _, table_path = _run_to_dir(tmp_path)
    df = pd.read_csv(str(table_path))

    assert list(df.columns) == ["product_smiles", "recipe"]
    assert len(df) == 1, f"Expected 1 row, got {len(df)}"

    row = df.iloc[0]
    assert row["product_smiles"] == "NC(CCO)CC1CCCCC1"
    for marker in ("|fragment|", "|seed|", "|rxn_smiles|", "|rxn_smarts|"):
        assert marker in row["recipe"], f"Missing marker {marker!r} in recipe"
    assert "[15*]C1CCCCC1" in row["recipe"]
    assert "[8*]CC(N)CCO" in row["recipe"]


def test_results_parquet_content(tmp_path):
    """brics_results.parquet has the same content as CSV."""
    _, table_path = _run_to_dir(tmp_path, fmt="parquet")
    df = pd.read_parquet(str(table_path))

    assert list(df.columns) == ["product_smiles", "recipe"]
    assert len(df) == 1
    assert df.iloc[0]["product_smiles"] == "NC(CCO)CC1CCCCC1"


def test_output_dir_creates_subdirs(tmp_path):
    """output_dir is created automatically if it doesn't exist."""
    deep = tmp_path / "sub" / "nested"
    apply_BRICS(SRC_SMILES, output_dir=str(deep), log_level=0)
    assert (deep / "decomposed_unique.txt").exists(), "Nested output_dir not created"
    assert (deep / "brics_results.csv").exists(), "Nested CSV not created"


# ── Edge cases ───────────────────────────────────────────────────────────────


def test_apply_BRICS_invalid_smiles(tmp_path):
    """Invalid SMILES are silently skipped; valid ones still processed."""
    # CCCOCC (ether) has BRICS-cleavable bonds; CCO (ethanol) does not.
    results = apply_BRICS(["not_a_smiles", "CCCOCC"], log_level=0)
    assert len(results) > 0, "Valid SMILES should produce results"
    all_text = " ".join(r["recipe"] for r in results) + " ".join(
        r["product_smiles"] for r in results
    )
    assert "CCCOCC" in all_text, "Valid SMILES should be processed"


def test_apply_BRICS_empty_list():
    """Empty input returns an empty list, no crash."""
    results = apply_BRICS([], log_level=0)
    assert results == [], "Empty input should yield empty results"


def test_apply_BRICS_empty_csv_output(tmp_path):
    """Empty input with output_dir writes an empty CSV (columns only)."""
    apply_BRICS([], output_dir=str(tmp_path), log_level=0)
    df = pd.read_csv(str(tmp_path / "brics_results.csv"))
    assert list(df.columns) == ["product_smiles", "recipe"]
    assert len(df) == 0


def test_apply_BRICS_return_is_none_with_output_dir(tmp_path):
    """When output_dir is given, the return value is None."""
    result = apply_BRICS(SRC_SMILES, output_dir=str(tmp_path), log_level=0)
    assert result is None, "Should return None when output_dir is provided"


def test_fmt_default_is_csv(tmp_path):
    """Default format is CSV."""
    apply_BRICS(SRC_SMILES, output_dir=str(tmp_path), log_level=0)
    assert (tmp_path / "brics_results.csv").exists(), "Default should be CSV"
    assert not (tmp_path / "brics_results.parquet").exists(), "Parquet not requested"


def test_fmt_parquet(tmp_path):
    """Explicit parquet format works."""
    apply_BRICS(SRC_SMILES, output_dir=str(tmp_path), fmt="parquet", log_level=0)
    assert (tmp_path / "brics_results.parquet").exists(), "Parquet should exist"
    assert not (tmp_path / "brics_results.csv").exists(), "CSV not requested"
