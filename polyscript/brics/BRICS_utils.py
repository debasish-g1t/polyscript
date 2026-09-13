from __future__ import annotations

import os
from typing import Optional

import pandas as pd
from polyscript.brics.BRICSMod import BRICSBuild, BRICSDecompose
from polyscript.utils.logger import INFO, PolyLogger
from rdkit import Chem


def read_df(df_path: str, column_name: str):
    """Read a column of data from a CSV file.

    Args:
        df_path: Path to a CSV file.
        column_name: Name of the column to extract.

    Returns:
        numpy.ndarray: Values from the specified column.

    Raises:
        ValueError: If the file is not a CSV or the column is missing.
    """
    if not df_path.split(".")[-1] == "csv":
        raise ValueError("The given file path is not targeting a csv fil")
    df = pd.read_csv(df_path)
    if not column_name in df.columns:
        raise ValueError("The given column is not in the df !")
    data = df[column_name].values
    return data


def write_data(data, file_path: str, new_line: bool = True):
    """Write data to a text file.

    Args:
        data: A list of strings (written one per line) or a single
            string (written as-is).
        file_path: Destination file path.
        new_line: If True (default) and *data* is a list, write each
            element followed by a newline.  Ignored for non-list data.

    Returns:
        bool: ``True`` on success.
    """
    with open(file_path, "w+") as f:
        if type(data) == list:
            if new_line:
                for d in data:
                    f.write(d + "\n")
            else:
                f.write(str(data))
        else:
            f.write(data)
    return True


def apply_BRICS(
    data: list[str],
    max_depth: int = 1,
    output_dir: Optional[str] = None,
    fmt: str = "csv",
    is_decomposed: bool = False,
    log_level: int | str = 2,
    logger: Optional[PolyLogger] = None,
) -> Optional[list[dict]]:
    """Decompose and rebuild molecules via BRICS.

    Parameters
    ----------
    data : list[str]
        SMILES strings to process.  May be full molecules (the default)
        or already-decomposed BRICS fragments when *is_decomposed* is
        ``True``.
    max_depth : int
        Maximum BRICS build depth (default 1).
    output_dir : str, optional
        Directory to save output files to.  When provided, two files are
        written: ``decomposed_unique.txt`` (fragments, one per line) and
        ``brics_results.{csv,parquet}`` (table with ``product_smiles`` and
        ``recipe`` columns).  When ``None`` (default), no files are written
        and the results are returned as a list of dicts instead.
    fmt : str
        Output table format: ``"csv"`` (default) or ``"parquet"``.
        Ignored when *output_dir* is ``None``.
    is_decomposed : bool
        When ``False`` (default), *data* contains full molecules that will
        be decomposed into BRICS fragments before building.  When ``True``,
        *data* already contains BRICS fragments (with wildcards like
        ``[15*]C1CCCCC1``) — the decomposition step is skipped and the
        fragments go directly to ``BRICSBuild``.
    log_level : int or str
        Verbosity threshold used when *logger* is not provided.
        As an **int**: ``0`` = silent, ``1`` = warnings,
        ``2`` = info (default), ``3`` = debug.
        As a **str**: ``"SILENT"``, ``"WARNING"``, ``"INFO"``, ``"DEBUG"``.
    logger : PolyLogger, optional
        A pre-configured ``PolyLogger`` instance.  When given, its
        existing level and handlers are used as-is; *log_level* is
        ignored.  When ``None``, a console logger at *log_level* is
        created automatically.

    Returns
    -------
    list[dict] or None
        When *output_dir* is ``None``, returns a list of dicts keyed by
        ``"product_smiles"``, ``"recipe"``, and ``"mol"``.  When
        *output_dir* is given, returns ``None`` (results are written to
        disk instead).
    """
    if logger is None:
        logger = PolyLogger(name="brics", level=log_level)
        logger.add_console(level=log_level)

    # ── 1. Decompose (skip if already decomposed) ───────────────────────
    if is_decomposed:
        # Input is already BRICS fragments — use directly, deduplicate
        uniques = list(set(data))
        logger.info(
            "Using %d pre-decomposed fragments (skipping decomposition)",
            len(uniques),
        )
    else:
        Decompositions = []
        for smiles in data:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                logger.warning(f"Invalid SMILES >> {smiles}")
            if mol is not None:
                res = BRICSDecompose(mol)
                Decompositions += res

        uniques = list(set(Decompositions))
        logger.info("Unique blocks found -> %d", len(uniques))

    # ── 2. Build ─────────────────────────────────────────────────────────
    build_mols = [Chem.MolFromSmiles(smiles) for smiles in uniques]
    brics_iter = list(BRICSBuild(build_mols, maxDepth=max_depth))

    # ── 3. Collect results ───────────────────────────────────────────────
    results: list[dict] = []
    for b in brics_iter:
        if isinstance(b, tuple):
            recipe, mol = b
            product_smiles = Chem.MolToSmiles(mol)
            results.append({
                "product_smiles": product_smiles,
                "recipe": recipe,
                "mol": mol,
            })

    logger.info("New blocks found -> %d", len(brics_iter))

    # ── 4. Save or return ────────────────────────────────────────────────
    if output_dir is not None:
        os.makedirs(output_dir, exist_ok=True)

        # Fragments — human-readable, one SMILES per line
        write_data(uniques, os.path.join(output_dir, "decomposed_unique.txt"))

        # Results table — dataframe-compatible (CSV or Parquet)
        if results:
            df_out = pd.DataFrame([
                {"product_smiles": r["product_smiles"], "recipe": r["recipe"]}
                for r in results
            ])
        else:
            df_out = pd.DataFrame(columns=["product_smiles", "recipe"])

        if fmt == "parquet":
            out_path = os.path.join(output_dir, "brics_results.parquet")
            df_out.to_parquet(out_path, index=False)
        else:
            out_path = os.path.join(output_dir, "brics_results.csv")
            df_out.to_csv(out_path, index=False)

        logger.info(
            "Saved %d results to %s (decomposed_unique.txt, %s)",
            len(results), output_dir, os.path.basename(out_path),
        )
        return None

    return results


if __name__ == "__main__":
    smiles = read_df("./source.csv", "SMILES")
    apply_BRICS(smiles)
