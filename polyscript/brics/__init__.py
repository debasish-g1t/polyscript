"""
BRICS (Breaking of Retrosynthetically Interesting Chemical Substructures)
sub-package for PolyScript.

Provides fragment-based decomposition and recombination of molecules:

- :func:`FindBRICSBonds`    — locate cleavable bonds in a molecule.
- :func:`BreakBRICSBonds`   — break those bonds, returning labelled fragments.
- :func:`BRICSDecompose`    — full recursive decomposition into BRICS fragments.
- :func:`BRICSBuild`        — enumerate new molecules by recombining fragments.
- :func:`BRICSBuildDistribWrite` — build + persist results to disk.

High-level convenience (from ``BRICS_utils``):
- :func:`apply_BRICS`       — end-to-end: decompose → deduplicate → rebuild.
- :func:`read_df`           — load a SMILES column from a CSV.
- :func:`write_data`        — write text / list data to a file.
"""

from __future__ import annotations

# ── Core BRICS operations ────────────────────────────────────────────────────
from polyscript.brics.BRICSMod import (
    BRICSBuild,
    BRICSBuildDistribWrite,
    BRICSDecompose,
    BreakBRICSBonds,
    FindBRICSBonds,
)

# ── High-level utilities ─────────────────────────────────────────────────────
from polyscript.brics.BRICS_utils import (
    apply_BRICS,
    read_df,
    write_data,
)

__all__ = [
    # BRICSMod
    "FindBRICSBonds",
    "BreakBRICSBonds",
    "BRICSDecompose",
    "BRICSBuild",
    "BRICSBuildDistribWrite",
    # BRICS_utils
    "apply_BRICS",
    "read_df",
    "write_data",
]
