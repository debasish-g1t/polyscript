"""Full-set smoke tests for every configured extender driver.

Generate the input with::

    python tests/extract_extender_test_data.py

Every row for a configured reaction is run through the driver selected by
``rxn_map``.  Failures are collected per reaction and reported together with
the offending monomers.  The former ring-opening gaps (a thioimidate lactam and
a cyclic carbamate/boronate lactone) are now covered by ``ring_opening``, so no
chemistry is excused; any failure fails the test.
"""

from pathlib import Path

import pandas as pd
import pytest
from rdkit import RDLogger

from polyscript.extender import (
    ChainGrowthReaction,
    RingOpeningMetathesisReaction,
    RingOpeningReaction,
    StepGrowthReaction,
)
from polyscript.extender.map import reaction_config, rxn_map

RDLogger.DisableLog("rdApp.*")

_CSV_PATH = (
    Path(__file__).resolve().parent
    / "test_data"
    / "extender_test_data_all.csv"
)

pytestmark = pytest.mark.skipif(
    not _CSV_PATH.exists(),
    reason=(
        "full extender set not generated; run "
        "`python tests/extract_extender_test_data.py`"
    ),
)

_DF = pd.read_csv(_CSV_PATH) if _CSV_PATH.exists() else pd.DataFrame()
_REACTIONS = (
    sorted(
        reaction
        for reaction in _DF["reaction_smarts"].unique()
        if rxn_map.get(reaction) is not None
    )
    if not _DF.empty
    else []
)

# The two former gaps (a thioimidate lactam and a cyclic carbamate/boronate
# lactone) are now handled by the ring-opening propagation templates.  Keep the
# set empty so any failure is reported as a new failure rather than excused.
KNOWN_COUPLING_GAPS = {}


def _monomers(row) -> list[str]:
    return [
        value
        for value in (row["mon1"], row["mon2"])
        if isinstance(value, str) and value.strip()
    ]


def _driver(reaction: str, config: dict):
    rxn_type = config["type"]
    if rxn_type == "addition":
        return ChainGrowthReaction(**config["params"], max_products=1)
    if rxn_type == "step_growth":
        return StepGrowthReaction(
            **reaction_config(reaction), max_products=1
        )
    if rxn_type == "ring_opening":
        return RingOpeningReaction(
            **reaction_config(reaction), max_products=1
        )
    if rxn_type == "romp":
        return RingOpeningMetathesisReaction(
            **reaction_config(reaction), max_products=1
        )
    raise AssertionError(f"Unknown reaction type: {rxn_type!r}")


@pytest.mark.parametrize(
    "reaction",
    _REACTIONS,
    ids=[f"{index}-{key[:20]}" for index, key in enumerate(_REACTIONS)],
)
def test_large_set_reaction(reaction):
    config = rxn_map.get(reaction)
    assert config is not None, f"No rxn_map entry for {reaction!r}"

    driver = _driver(reaction, config)
    rows = _DF[_DF["reaction_smarts"] == reaction]
    failures: list[tuple] = []
    known_failures: list[tuple] = []

    for _, row in rows.iterrows():
        try:
            products = driver.extend(_monomers(row), 2)
            assert products, "driver returned no products"
            assert all("*" not in smiles for smiles in products), (
                "driver returned a wildcard-terminated product"
            )
        except Exception as exc:  # noqa: BLE001 - report every candidate
            monomer = str(row["mon1"])
            if monomer in KNOWN_COUPLING_GAPS:
                known_failures.append(
                    (row.get("row_index"), monomer,
                     KNOWN_COUPLING_GAPS[monomer])
                )
                continue
            failures.append(
                (
                    row.get("row_index"),
                    row["mon1"],
                    row["mon2"],
                    f"{type(exc).__name__}: {exc}",
                )
            )

    assert not failures, (
        f"{len(failures)}/{len(rows)} candidates failed for "
        f"{reaction[:60]}...\n"
        + "\n".join(
            f"  row={row_index} mon1={mon1} mon2={mon2} error={error}"
            for row_index, mon1, mon2, error in failures[:5]
        )
    )
