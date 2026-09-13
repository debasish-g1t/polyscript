#!/usr/bin/env python3
"""
MCP Server for PolyScript — polymer chemistry sequence generation, parsing,
classification, validation, and depolymerization.

Exposes tools via FastMCP for AI agents to work with PolyScript sequences:
  - Parse PolyScript into monomers/reaction/polymer
  - Validate chemical correctness of a PolyScript sequence

Note: WE DO NOT PROVIDE EXACT EXAMPLES IN THE DOCUMENTATION, BECAUSE THE LLM REPLICATES THE EXAMPLES
OR HALLUCINATES BASED ON THE PROPER POLYSCRIPT EXAMPLE RATHER GENERATING NOVEL SEQUENCES.
"""


  # - Classify monomers by functional-group type
  # - Polymerize: generate polymer sequences from monomer pairs
  # - Depolymerize: recover monomers from polymer SMILES via reverse reaction

from __future__ import annotations

import json
import os
import tempfile
from typing import Optional

import pandas as pd
from mcp.server.fastmcp import FastMCP

# ── PolyScript imports ───────────────────────────────────────────────────────
from polyscript.classifier import PolyScriptClassifier
from polyscript.depolymerizer import PolyScriptDepolymerizer, DepolymerizeResult
from polyscript.polymerizer import Polymerizer
from polyscript.utils.parsers import PolyScriptParser, NestedPolyScriptParser
from polyscript.utils.validators import SeqValidator
from polyscript.utils.executor import BaseExecutor

# ---------------------------------------------------------------------------
# Global server
# ---------------------------------------------------------------------------
mcp = FastMCP("polyscript")

# ── Lazy / cached singletons ────────────────────────────────────────────────
_classifier: PolyScriptClassifier | None = None
_depolymerizer: PolyScriptDepolymerizer | None = None
_parser: PolyScriptParser | None = None
_nested_parser: NestedPolyScriptParser | None = None
_validator: SeqValidator | None = None
_polymerizer: Polymerizer | None = None


def _quiet() -> None:
    """Mute all constructors so they don't print to stderr by default."""
    BaseExecutor.set_log_level("SILENT")


def _get_parser() -> PolyScriptParser:
    global _parser
    if _parser is None:
        _quiet()
        _parser = PolyScriptParser()
    return _parser


def _get_nested_parser() -> NestedPolyScriptParser:
    global _nested_parser
    if _nested_parser is None:
        _quiet()
        _nested_parser = NestedPolyScriptParser()
    return _nested_parser


def _get_validator() -> SeqValidator:
    global _validator
    if _validator is None:
        _quiet()
        _validator = SeqValidator()
    return _validator


def _get_classifier() -> PolyScriptClassifier:
    global _classifier
    if _classifier is None:
        _quiet()
        _classifier = PolyScriptClassifier(include_carbonate=False)
    return _classifier


def _get_depolymerizer() -> PolyScriptDepolymerizer:
    global _depolymerizer
    if _depolymerizer is None:
        _quiet()
        _depolymerizer = PolyScriptDepolymerizer()
    return _depolymerizer


def _get_polymerizer() -> Polymerizer:
    global _polymerizer
    if _polymerizer is None:
        _quiet()
        _polymerizer = Polymerizer()
    return _polymerizer


def _depol_result_to_dict(r: DepolymerizeResult) -> dict:
    return {
        "polymer": r.polymer,
        "reaction": r.reaction,
        "monomers": r.monomers,
        "input_monomers": r.input_monomers,
        "polymer_type": r.polymer_type,
        "polymer_pattern": r.polymer_pattern,
        "seq_index": r.seq_index,
        "is_exact_match": r.is_exact_match,
    }


# ═══════════════════════════════════════════════════════════════════════════════
# Tools
# ═══════════════════════════════════════════════════════════════════════════════


# ── Parse ────────────────────────────────────────────────────────────────────


@mcp.tool()
def parse_single_sequence(sequence: str) -> dict:
    """
    PolyScript is a library that handles a specific form of polymer
    representation expressed in a character sequence.

    This function parses a **single** PolyScript sequence into its
    components.  It does **not** handle multiple sequences (use
    ``parse_nested_sequence`` or ``batch_validate_sequences`` for
    ``|||``-separated inputs).

    .. note::

       ``{...}`` is literal PolyScript syntax (required).  ``<...>``
       tags are placeholders — replace them with actual SMILES/SMARTS.

    .. important::

       The polymer **must** be a PSMILES containing ``*`` as wildcard
       atoms (polymer attachment points).  The ``*`` positions must be
       consistent with the atom mappings in the reaction SMARTS.

    Format
    ------
    ``{<monomer1>}+{<monomer2>}=>{<reaction>}=>{<polymer>}``

    Example
    -------
    Valid
    input: "{<monomer1>}+{<monomer2>}=>{<reaction>}=>{<polymer>}"
    output: {"monomers": ["<monomer1>", "<monomer2>"], "reaction": "<reaction>", "polymer": "<polymer>", "errors": [], "valid": true}

    Invalid
    input: "{<monomer1>}+{<monomer2>}<bad_delim>{<reaction>}=>{<polymer>}"
    output: {"monomers": null, "reaction": null, "polymer": null, "errors": ["<E-parse-|invalid-number-of-parts|>"], "valid": false}
    """
    parser = _get_parser()
    errors = parser.parse(sequence)
    return {
        "monomers": parser.monomers,
        "reaction": parser.reaction,
        "polymer": parser.polymer,
        "errors": errors,
        "valid": len(errors) == 0,
    }


@mcp.tool()
def parse_nested_sequence(sequence: str) -> dict:
    """Parse multiple PolyScript sub-sequences separated by '|||'.

    Each sub-sequence is parsed independently.  Returns a list of
    (monomers, reaction, polymer) tuples and per-sub-sequence errors.

    .. note::

       In the examples below ``{...}`` is literal PolyScript syntax
       (required).  ``<...>`` tags are placeholders — they are **not**
       part of the format.  Replace them with actual SMILES/SMARTS.

    .. important::

       The polymer **must** be a PSMILES containing ``*`` wildcards
       consistent with the reaction SMARTS atom mappings.

    Format
    ------
    ``{<monomer1>}+{<monomer2>}=>{<reaction>}=>{<polymer>}|||...``

    Example
    -------
    Valid
    input: "{<monomer1>}+{<monomer2>}=>{<rxn1>}=>{<poly1>}|||{<poly1>}+{none}=>{<rxn2>}=>{<poly2>}"
    output: {"sequences": [{"monomers": ["<monomer1>", "<monomer2>"], "reaction": "<rxn1>", "polymer": "<poly1>"}, {"monomers": ["<poly1>", "none"], "reaction": "<rxn2>", "polymer": "<poly2>"}], "errors": [[], []], "valid": true}

    Invalid
    input: "{<monomer1>}+{<monomer2>}=>{<rxn1>}=>{<poly1>}|||{<poly1>}+{none}=>{<rxn2>}=>{<poly2>}"
    output: {"sequences": [{"monomers": null, "reaction": null, "polymer": null}], "errors": [["<E-parse-|too-many-parts|>"]], "valid": false}
    """
    parser = _get_nested_parser()
    errors = parser.parse(sequence)
    return {
        "sequences": [
            {
                "monomers": mons,
                "reaction": rxn,
                "polymer": poly,
            }
            for mons, rxn, poly in parser.sequences
        ],
        "errors": errors,
        "valid": all(len(e) == 0 for e in errors) if errors else False,
    }


# ── Classify ─────────────────────────────────────────────────────────────────


# @mcp.tool()
# def classify_monomers(smiles_list: list[str]) -> dict:
#     """Classify each monomer SMILES by its functional-group type(s).

#     Returns a table (JSON-serializable) with boolean columns for each
#     monomer functional-group class.

#     Example
#     -------
#     input: ["<monomer_smiles_1>", "<monomer_smiles_2>"]
#     output: {"columns": ["smiles", "<type1>", "<type2>", ...], "rows": [{"smiles": "<monomer_smiles_1>", "<type1>": true, "<type2>": false, ...}, {"smiles": "<monomer_smiles_2>", "<type1>": false, "<type2>": true, ...}], "num_rows": 2}
#     """
#     classifier = _get_classifier()
#     df = pd.DataFrame({"smiles": smiles_list})
#     result = classifier.classify(df, col_name="smiles")

#     # Drop ROMol if present (non-serializable)
#     if "ROMol" in result.columns:
#         result = result.drop(columns=["ROMol"])

#     return {
#         "columns": result.columns.tolist(),
#         "rows": result.to_dict(orient="records"),
#         "num_rows": len(result),
#     }


# ── Validate ─────────────────────────────────────────────────────────────────


@mcp.tool()
def validate_sequence(sequence: str) -> dict:
    """Chemically validate a PolyScript sequence.

    Checks:
      - Monomer SMILES validity
      - Reaction SMARTS validity
      - Reaction execution (whether the reaction runs on the monomers)
      - Polymer product match

    Handles multi-sequences (``|||``).  Returns errors (empty = valid).

    .. note::

       In the examples below ``{...}`` is literal PolyScript syntax
       (required).  ``<...>`` tags are placeholders — they are **not**
       part of the format.  Replace them with actual SMILES/SMARTS.

    .. important::

       The polymer **must** be a PSMILES containing ``*`` wildcards
       consistent with the reaction SMARTS atom mappings.

    Format
    ------
    ``{<monomer1>}+{<monomer2>}=>{<reaction>}=>{<polymer>}``

    Example
    -------
    Valid
    input: "{<monomer1>}+{<monomer2>}=>{<reaction>}=>{<polymer>}"
    output: {"errors": [], "valid": true}

    Invalid
    input: "{<bad_monomer>}+{<monomer2>}=>{<reaction>}=>{<polymer>}"
    output: {"errors": ["<E-chem-|invalid-monomer-smiles|>", "<E-chem-|reaction-execution-error|>"], "valid": false}

    """
    validator = _get_validator()
    errors = validator.validate(sequence)
    return {
        "errors": errors,
        "valid": len(errors) == 0,
    }


@mcp.tool()
def is_valid_sequence(sequence: str) -> bool:
    """
    Return True when every sub-sequence passes all chemical validation.

    .. note::

       In the examples below ``{...}`` is literal PolyScript syntax
       (required).  ``<...>`` tags are placeholders — they are **not**
       part of the format.  Replace them with actual SMILES/SMARTS.

    .. important::

       The polymer **must** be a PSMILES containing ``*`` wildcards
       consistent with the reaction SMARTS atom mappings.

    Format
    ------
    ``{<monomer1>}+{<monomer2>}=>{<reaction>}=>{<polymer>}``

    Example
    -------
    Valid
    input: "{<monomer1>}+{<monomer2>}=>{<reaction>}=>{<polymer>}"
    output: true

    Invalid
    input: "{<bad_monomer>}+{<monomer2>}=>{<reaction>}=>{<polymer>}"
    output: false
    """

    return _get_validator().is_valid(sequence)


# ── Batch validate (includes parse) ──────────────────────────────────────


@mcp.tool()
def batch_validate_sequences(sequences: list[str]) -> dict:
    """Parse and chemically validate multiple PolyScript sequences.

    Accepts a list of PolyScript sequences and returns both parse *and*
    chemical validation results for every entry in a single round-trip.
    Handles single and nested (``|||``-separated) sequences automatically.

    Much more efficient than calling the single-sequence tools separately.

    .. note::

       ``{...}`` is literal PolyScript syntax (required).  ``<...>``
       tags are placeholders — replace them with actual SMILES/SMARTS.

    .. important::

       The polymer **must** be a PSMILES containing ``*`` wildcards
       consistent with the reaction SMARTS atom mappings.

    Format
    ------
    ``{<monomer1>}+{<monomer2>}=>{<reaction>}=>{<polymer>}``
    (or ``|||``-separated chains of the above)

    Each result entry contains:
      - ``index`` — position in the input list
      - ``monomers``, ``reaction``, ``polymer`` — top-level parsed
        components (first valid sub-sequence for nested inputs)
      - ``sub_sequences`` — for nested inputs, each ``|||``-separated
        sub-sequence as ``{sub_index, monomers, reaction, polymer}``
      - ``parse_errors`` — list of parse error tags
      - ``validation_errors`` — list of chemical validation error tags
      - ``valid`` — true when both parse and validation pass

    Example
    -------
    Input: [
      "{<monomer1>}+{<monomer2>}=>{<reaction>}=>{<polymer>}",
      "not-a-valid-sequence"
    ]
    output: {"results": [{"index": 0, "monomers": ["<monomer1>", "<monomer2>"], "reaction": "<reaction>", "polymer": "<polymer>", "sub_sequences": [], "parse_errors": [], "validation_errors": [], "valid": true}, {"index": 1, "monomers": null, "reaction": null, "polymer": null, "sub_sequences": [], "parse_errors": ["<E-parse-|invalid-number-of-parts|>"], "validation_errors": ["<E-parse-|invalid-number-of-parts|>"], "valid": false}], "total": 2, "valid_count": 1}
    """
    parser = _get_parser()
    nested_parser = _get_nested_parser()
    validator = _get_validator()

    results = []
    for i, seq in enumerate(sequences):
        entry: dict = {
            "index": i,
            "monomers": None,
            "reaction": None,
            "polymer": None,
            "sub_sequences": [],
            "parse_errors": [],
            "validation_errors": [],
            "valid": False,
        }

        # ── Parse ─────────────────────────────────────────────────────
        if "|||" in seq:
            sub_errors = nested_parser.parse(seq)
            entry["parse_errors"] = [e for sub in sub_errors for e in sub]
            first_valid = None
            for j, (mons, rxn, poly) in enumerate(nested_parser.sequences):
                entry["sub_sequences"].append({
                    "sub_index": j,
                    "monomers": mons,
                    "reaction": rxn,
                    "polymer": poly,
                })
                if first_valid is None and mons is not None:
                    first_valid = (mons, rxn, poly)
            if first_valid is not None:
                entry["monomers"] = first_valid[0]
                entry["reaction"] = first_valid[1]
                entry["polymer"] = first_valid[2]
        else:
            entry["parse_errors"] = parser.parse(seq)
            entry["monomers"] = parser.monomers
            entry["reaction"] = parser.reaction
            entry["polymer"] = parser.polymer

        # ── Validate ───────────────────────────────────────────────────
        try:
            entry["validation_errors"] = validator.validate(seq)
        except Exception:
            entry["validation_errors"] = ["<E-batch-|validation-crash|>"]

        entry["valid"] = (
            len(entry["parse_errors"]) == 0
            and len(entry["validation_errors"]) == 0
        )
        results.append(entry)

    return {
        "results": results,
        "total": len(results),
        "valid_count": sum(1 for r in results if r["valid"]),
    }


# ── Post-process candidates ───────────────────────────────────────────────


@mcp.tool()
def post_process_oneshot(candidates: list[dict]) -> dict:
    """
    Parse and chemically validate a batch of polymer candidates
    in a single call.  Accepts a list of candidate objects each
    with a ``polyscript_repr`` field.

    **This is for one-shot post-processing, not iterative calling.**
    Use ``validate_sequence`` or ``batch_validate_sequences`` for
    single or iterative validation during generation.

    .. note::

       ``{...}`` is literal PolyScript syntax (required).  ``<...>``
       tags are placeholders — replace them with actual SMILES/SMARTS.

    .. important::

       The polymer **must** be a PSMILES containing ``*`` wildcards
       consistent with the reaction SMARTS atom mappings.

    Format
    ------
    ``{<monomer1>}+{<monomer2>}=>{<reaction>}=>{<polymer>}``
    (or ``|||``-separated chains of the above)

    Each result entry:
      - ``index`` — position in the input list
      - ``monomers``, ``reaction``, ``polymer`` — parsed components
      - ``sub_sequences`` — nested sub-sequence details
      - ``parse_errors`` — parse error tags
      - ``validation_errors`` — chemical validation error tags
      - ``valid`` — true when both parse and validation pass

    Example
    -------
    Input: {"candidates": [
      {"polyscript_repr": "{<monomer1>}+{<monomer2>}=>{<reaction>}=>{<polymer>}"}
    ]}
    output: {"results": [{"index": 0, "monomers": ["<monomer1>", "<monomer2>"], "reaction": "<reaction>", "polymer": "<polymer>", "sub_sequences": [], "parse_errors": [], "validation_errors": [], "valid": true}], "total": 1, "valid_count": 1}
    """
    parser = _get_parser()
    nested_parser = _get_nested_parser()
    validator = _get_validator()

    results = []
    for i, cand in enumerate(candidates):
        seq = cand.get("polyscript_repr", "")

        entry: dict = {
            "index": i,
            "monomers": None,
            "reaction": None,
            "polymer": None,
            "sub_sequences": [],
            "parse_errors": [],
            "validation_errors": [],
            "valid": False,
        }

        # ── Parse ─────────────────────────────────────────────────────
        if "|||" in seq:
            sub_errors = nested_parser.parse(seq)
            entry["parse_errors"] = [e for sub in sub_errors for e in sub]
            first_valid = None
            for j, (mons, rxn, poly) in enumerate(nested_parser.sequences):
                entry["sub_sequences"].append({
                    "sub_index": j,
                    "monomers": mons,
                    "reaction": rxn,
                    "polymer": poly,
                })
                if first_valid is None and mons is not None:
                    first_valid = (mons, rxn, poly)
            if first_valid is not None:
                entry["monomers"] = first_valid[0]
                entry["reaction"] = first_valid[1]
                entry["polymer"] = first_valid[2]
        else:
            entry["parse_errors"] = parser.parse(seq)
            entry["monomers"] = parser.monomers
            entry["reaction"] = parser.reaction
            entry["polymer"] = parser.polymer

        # ── Validate ───────────────────────────────────────────────────
        try:
            entry["validation_errors"] = validator.validate(seq)
        except Exception:
            entry["validation_errors"] = ["<E-batch-|validation-crash|>"]

        entry["valid"] = (
            len(entry["parse_errors"]) == 0
            and len(entry["validation_errors"]) == 0
        )
        results.append(entry)

    return {
        "results": results,
        "total": len(results),
        "valid_count": sum(1 for r in results if r["valid"]),
    }


# # ── Polymerize ───────────────────────────────────────────────────────────────


# @mcp.tool()
# def polymerize(
#     monomer_smiles_1: list[str],
#     monomer_type_1: str,
#     polymer_class: str,
#     monomer_smiles_2: Optional[list[str]] = None,
#     monomer_type_2: Optional[str] = None,
# ) -> dict:
#     """Generate polymer sequences from one or two monomer lists.

#     When only one monomer list is provided, homopolymerization is performed.
#     When two monomer lists are provided, bipolymerization (copolymer) is
#     performed using the matching reaction rule.

#     Args:
#         monomer_smiles_1: SMILES strings for the first set of monomers.
#         monomer_type_1: Functional-group type tag (e.g. "vinyl", "acid").
#         polymer_class: Polymer class (e.g. "polyolefin", "polyester").
#         monomer_smiles_2: Optional second monomer list for copolymerization.
#         monomer_type_2: Functional-group type for second monomers (required
#             when monomer_smiles_2 is provided).
#     """
#     poly = _get_polymerizer()

#     mon_df1 = pd.DataFrame({"smiles": monomer_smiles_1})

#     mon_df2 = None
#     if monomer_smiles_2 is not None:
#         mon_df2 = pd.DataFrame({"smiles": monomer_smiles_2})

#     result = poly.bipolymerize(
#         mon_df1=mon_df1,
#         mon_type1=monomer_type_1,
#         P_class=polymer_class,
#         mon_df2=mon_df2,
#         mon_type2=monomer_type_2,
#         candidate_col_name="smiles",
#     )

#     if result is None or result.empty:
#         return {
#             "num_results": 0,
#             "results": [],
#             "message": "No polymer sequences were generated for the given inputs.",
#         }

#     return {
#         "num_results": len(result),
#         "columns": result.columns.tolist(),
#         "results": result.to_dict(orient="records"),
#     }


# @mcp.tool()
# def execute_polymerization(
#     mon1_path: str,
#     mon1_type: str,
#     polymer_class: str,
#     col_name: str = "psmiles",
#     mon2_path: Optional[str] = None,
#     mon2_type: Optional[str] = None,
#     output_path: Optional[str] = None,
# ) -> dict:
#     """Run polymerizer.execute() with a path-based config.

#     Reads monomer parquet files from disk and writes results.

#     Args:
#         mon1_path: Path to first monomer parquet file.
#         mon1_type: Functional-group type for first monomers.
#         polymer_class: Polymer class (e.g. "polyolefin").
#         col_name: Column name containing SMILES in the parquet files.
#         mon2_path: Optional second monomer parquet file.
#         mon2_type: Functional-group type for second monomers.
#         output_path: Where to write the output parquet.  If omitted, a
#             temporary file is used.
#     """
#     if output_path is None:
#         with tempfile.NamedTemporaryFile(suffix=".parquet", delete=False) as f:
#             output_path = f.name

#     config = {
#         "mon1_path": mon1_path,
#         "mon2_path": mon2_path,
#         "output_path": output_path,
#         "col_name": col_name,
#         "P_class": polymer_class,
#         "mon_type1": mon1_type,
#         "mon_type2": mon2_type,
#     }

#     poly = _get_polymerizer()
#     success, error = poly.execute(config)

#     return {
#         "success": success,
#         "error": error,
#         "output_path": output_path if success else None,
#     }


# # ── Depolymerize ─────────────────────────────────────────────────────────────


# @mcp.tool()
# def depolymerize(polymer_smiles: str) -> dict:
#     """Recover monomers from a polymer SMILES via reverse reaction.

#     Finds all reaction pathways that can produce the given polymer,
#     then reverses them to recover candidate monomers.  Returns the
#     list of successful recovery results.
#     """
#     depol = _get_depolymerizer()
#     results = depol.convert_psmiles(polymer_smiles)
#     return {
#         "polymer": polymer_smiles,
#         "num_pathways": len(results),
#         "pathways": [_depol_result_to_dict(r) for r in results],
#     }


# @mcp.tool()
# def depolymerize_batch(polymer_smiles_list: list[str]) -> dict:
#     """Batch depolymerization of multiple polymer SMILES.

#     Finds reverse -> forward reaction round-trips for each polymer.
#     """
#     depol = _get_depolymerizer()
#     results, stats = depol.convert_all(polymer_smiles_list)
#     return {
#         "stats": stats,
#         "num_pathways": len(results),
#         "pathways": [
#             _depol_result_to_dict(r) for r in results
#         ],
#     }


# @mcp.tool()
# def validate_depolymerization(polym_sequence: str) -> dict:
#     """Validate monomer recovery via depolymerization for a PolyScript sequence.

#     Parses the PolyScript sequence, then attempts a reverse -> forward
#     reaction round-trip to check whether the input monomers can be
#     recovered from the polymer.
#     """
#     depol = _get_depolymerizer()
#     results = depol.validate(polym_sequence)
#     return {
#         "sequence": polym_sequence,
#         "num_matches": len(results),
#         "matches": [
#             _depol_result_to_dict(r) for r in results
#         ],
#     }


# @mcp.tool()
# def validate_depolymerization_batch(polym_sequences: list[str]) -> dict:
#     """Batch depolymerization validation for multiple PolyScript sequences."""
#     depol = _get_depolymerizer()
#     results, stats = depol.validate_all(polym_sequences)
#     return {
#         "stats": stats,
#         "num_matches": len(results),
#         "matches": [
#             _depol_result_to_dict(r) for r in results
#         ],
#     }


# ═══════════════════════════════════════════════════════════════════════════════
# Entry point
# ═══════════════════════════════════════════════════════════════════════════════


def main():
    """Run the MCP server.

    Transport is chosen via CLI arguments or the ``MCP_TRANSPORT`` env var:

    - ``stdio``  (default) — standard input/output, no port needed.
    - ``sse``    — Server-Sent Events over HTTP (default port 8000).
    - ``streamable-http`` — Streamable HTTP (default port 8000).

    Examples::

        python server.py                           # stdio
        python server.py --transport sse --port 9000
        MCP_TRANSPORT=sse python server.py         # port 8000
    """
    import argparse

    parser = argparse.ArgumentParser(description="PolyScript MCP Server")
    parser.add_argument(
        "--transport",
        choices=["stdio", "sse", "streamable-http"],
        default=os.environ.get("MCP_TRANSPORT", "stdio"),
        help="Transport protocol (default: stdio, or $MCP_TRANSPORT)",
    )
    parser.add_argument(
        "--host",
        default=os.environ.get("MCP_HOST", "127.0.0.1"),
        help="Host to bind when using sse / streamable-http (default: 127.0.0.1)",
    )
    parser.add_argument(
        "--port",
        type=int,
        default=int(os.environ.get("MCP_PORT", "8000")),
        help="Port to bind when using sse / streamable-http (default: 8000)",
    )
    args = parser.parse_args()

    if args.transport == "stdio":
        mcp.run(transport="stdio")
    elif args.transport == "sse":
        mcp.run(transport="sse", host=args.host, port=args.port)
    elif args.transport == "streamable-http":
        mcp.run(transport="streamable-http", host=args.host, port=args.port)


if __name__ == "__main__":
    main()
