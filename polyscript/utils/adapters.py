# This is still experimental and requires to be tested and modified later

import re
from typing import Optional

import pandas as pd
import selfies as sf
from polyscript.utils.parsers import PolyScriptParser
from polyscript.utils.logger import INFO, PolyLogger


def extract_bigsmiles(psmiles_string: str, directed: bool = False) -> str:
    """
    Converts a Polymer SMILES (PSMILES) string to a formally valid BigSMILES representation.

    Args:
        psmiles_string (str): The PSMILES string (e.g., "[*]CCO[*]", "*CC*", or "[1*]CC[2*]").
        directed (bool): If True, assigns head-to-tail directionality using [>] and [<].
                         If False, uses the standard undirected descriptor [$].

    Returns:
        str: A structurally valid BigSMILES string.
    """
    # 1. Validate the input contains at least one attachment point
    if not psmiles_string or not re.search(r"\[\d*\*\]|\*", psmiles_string):
        raise ValueError("Invalid PSMILES: No attachment points ('*' or '[*]') found.")

    # 2. Normalize all PSMILES attachment points to a temporary placeholder.
    # The regex \[\d*\*\] captures [*], [1*], [2*], etc. The \* captures a raw *.
    core_string = re.sub(r"\[\d*\*\]|\*", "{TMP}", psmiles_string)

    # 3. Apply the correct BigSMILES bonding descriptors
    if directed:
        # Head-to-tail mapping requires exactly 2 attachment points in the repeat unit
        if core_string.count("{TMP}") != 2:
            raise ValueError(
                "Directed conversion requires exactly 2 attachment points in the repeat unit."
            )

        # Replace the first attachment with the "left" descriptor [>]
        core_string = core_string.replace("{TMP}", "[>]", 1)
        # Replace the second attachment with the "right" descriptor [<]
        core_string = core_string.replace("{TMP}", "[<]", 1)
    else:
        # Undirected uses [$] for all connection points
        core_string = core_string.replace("{TMP}", "[$]")

    # 4. Wrap in the BigSMILES stochastic object syntax {...}
    # Note: '[]' denotes an empty/unspecified terminal connection for the polymer ends.
    bigsmiles = f"{{[]{core_string}[]}}"

    return bigsmiles


class Adapter:
    """Abstract base class for format adapters.

    Subclasses implement ``convert`` to transform a string from one
    chemical notation into another (e.g. PolyScript → PSMILES).
    """

    def __init__(self, type):
        """Initialise the adapter with a format type label.

        Args:
            type: A string identifier for the target format
                (e.g. ``"psmiles"``, ``"smiles"``).
        """
        self.type = type

    def convert(self, custom_string: str) -> str:
        """Convert *custom_string* to the target format.

        Args:
            custom_string: Input string in a supported source format.

        Returns:
            The converted string.

        Raises:
            NotImplementedError: Must be overridden by subclasses.
        """
        raise NotImplementedError


class PSMILESAdapter(Adapter):
    """Adapter that converts various notations into PSMILES (Polymer SMILES)."""

    def __init__(self):
        super().__init__("psmiles")

    def convert(self, custom_string: str, _from="polyscript") -> str:
        """Convert *custom_string* to PSMILES.

        Args:
            custom_string: The input string.
            _from: Source format — ``"polyscript"`` or ``"selfies"``.

        Returns:
            A PSMILES string.

        Raises:
            ValueError: If the input cannot be parsed or the source
                format is unknown.
        """
        if _from == "polyscript":
            parser = PolyScriptParser()
            err = parser.parse(custom_string)
            if len(err) > 0 or parser.polymer is None:
                raise ValueError(err)
            else:
                return parser.polymer
        elif _from == "selfies":
            decoded = sf.decoder(custom_string)
            if decoded is None:
                raise ValueError("sf.decoder returned None")
            return decoded
        else:
            raise ValueError(f"Unknown _from: {_from}")


class SMILESAdapter(Adapter):
    """Adapter that converts various notations into dot-separated SMILES."""

    def __init__(self):
        super().__init__("smiles")

    def convert(self, custom_string: str, _from="polyscript") -> str:
        """Convert *custom_string* to dot-separated monomer SMILES.

        Args:
            custom_string: The input string.
            _from: Source format — ``"polyscript"`` or ``"selfies"``.

        Returns:
            A dot-separated SMILES string of the monomers.

        Raises:
            ValueError: If the input cannot be parsed or the source
                format is unknown.
        """
        if _from == "polyscript":
            parser = PolyScriptParser()
            err = parser.parse(custom_string)
            if len(err) > 0 or parser.monomers is None:
                raise ValueError(err)
            else:
                return ".".join(parser.monomers)
        elif _from == "selfies":
            decoded = sf.decoder(custom_string)
            if decoded is None:
                raise ValueError("sf.decoder returned None")
            return decoded
        else:
            raise ValueError(f"Unknown _from: {_from}")


class SELFIESAdapter(Adapter):
    """Adapter that converts various notations into SELFIES."""

    def __init__(self):
        super().__init__("selfies")

    def convert(self, custom_string: str, _from="polyscript") -> str:
        """Convert *custom_string* to a SELFIES string.

        Args:
            custom_string: The input string.
            _from: Source format — ``"polyscript"``, ``"smiles"``, or
                ``"psmiles"``.

        Returns:
            A SELFIES-encoded string.

        Raises:
            ValueError: If the input cannot be parsed or the source
                format is unsupported.
        """
        if _from == "polyscript":
            parser = PolyScriptParser()
            err = parser.parse(custom_string)
            if len(err) > 0 or parser.polymer is None:
                raise ValueError(err)
            else:
                smiles_for_selfies = parser.polymer.replace("[*]", "[Y]").replace(
                    "*", "[Y]"
                )
                encoded = sf.encoder(smiles_for_selfies)
                if encoded is None:
                    raise ValueError("sf.encoder returned None")
                return encoded
        elif _from == "smiles":
            encoded = sf.encoder(custom_string)
            if encoded is None:
                raise ValueError("sf.encoder returned None")
            return encoded
        elif _from == "psmiles":
            smiles_for_selfies = custom_string.replace("[*]", "[Y]").replace("*", "[Y]")
            encoded = sf.encoder(smiles_for_selfies)
            if encoded is None:
                raise ValueError("sf.encoder returned None")
            return encoded
        else:
            raise ValueError(f"Unsupported conversion from {_from}")


class BigSMILESAdapter(Adapter):
    """Adapter that converts various notations into BigSMILES."""

    def __init__(self):
        super().__init__("bigsmiles")

    def convert(self, custom_string: str, _from="polyscript") -> str:
        """Convert *custom_string* to a BigSMILES string.

        Args:
            custom_string: The input string.
            _from: Source format — ``"polyscript"`` or ``"psmiles"``.

        Returns:
            A structurally valid BigSMILES string.

        Raises:
            ValueError: If the input cannot be parsed or the source
                format is unknown.
        """
        if _from == "polyscript":
            parser = PolyScriptParser()
            err = parser.parse(custom_string)
            if len(err) > 0 or parser.polymer is None:
                raise ValueError(err)
            else:
                return extract_bigsmiles(parser.polymer)
        elif _from == "psmiles":
            return extract_bigsmiles(custom_string)
        else:
            raise ValueError(f"Unsupported conversion from {_from}")


def extract_psmiles(polym: str) -> str:
    """Convert a PolyScript string to PSMILES."""
    return PSMILESAdapter().convert(polym, _from="polyscript")


def extract_smiles(polym: str) -> str:
    """Convert a PolyScript string to dot-separated monomer SMILES."""
    return SMILESAdapter().convert(polym, _from="polyscript")


def extract_selfies(polym: str) -> str:
    """Convert a PolyScript string to SELFIES."""
    return SELFIESAdapter().convert(polym, _from="polyscript")


def enrich(
    df: pd.DataFrame,
    progress_batch: int = 5_000,
    logger: Optional[PolyLogger] = None,
) -> pd.DataFrame:
    """Add ``psmiles``, ``smiles``, ``selfies``, ``bigsmiles`` columns."""
    if logger is None:
        logger = PolyLogger(name="enrich", level=INFO)
        logger.add_console(level=INFO)

    total = len(df)
    logger.info("Enriching %d rows …", total)

    # Step 1: polym → psmiles + smiles (both share PolyScriptParser.parse)
    psmiles_adapter = PSMILESAdapter()
    smiles_adapter = SMILESAdapter()

    psmiles_col: list[str | None] = []
    smiles_col: list[str | None] = []
    err_psmiles = 0
    err_smiles = 0

    for idx, polym in enumerate(df["polym"]):
        try:
            psmiles_col.append(psmiles_adapter.convert(polym, _from="polyscript"))
        except Exception:
            psmiles_col.append(None)
            err_psmiles += 1

        try:
            smiles_col.append(smiles_adapter.convert(polym, _from="polyscript"))
        except Exception:
            smiles_col.append(None)
            err_smiles += 1

        if (idx + 1) % progress_batch == 0:
            logger.info(
                "[psmiles/smiles]  %d/%d  (%.1f%%)  err: p=%d s=%d",
                idx + 1, total, (idx + 1) / total * 100, err_psmiles, err_smiles,
            )

    df["psmiles"] = psmiles_col
    # Use "monomer_smiles" to avoid clobbering any existing "smiles" column
    # (e.g. reference CSVs already have a "smiles" column with the PSMILES)
    df["monomer_smiles"] = smiles_col

    if err_psmiles:
        logger.warning("psmiles parse failures: %d", err_psmiles)
    if err_smiles:
        logger.warning("monomer_smiles parse failures: %d", err_smiles)

    # Step 2: psmiles → selfies (reuses cached psmiles, avoids re-parsing)
    selfies_adapter = SELFIESAdapter()
    selfies_col: list[str | None] = []
    err_selfies = 0

    for idx, psm in enumerate(psmiles_col):
        if psm is None:
            selfies_col.append(None)
            continue
        try:
            selfies_col.append(selfies_adapter.convert(psm, _from="psmiles"))
        except Exception:
            selfies_col.append(None)
            err_selfies += 1

        if (idx + 1) % progress_batch == 0:
            logger.info(
                "[selfies]         %d/%d  (%.1f%%)  err: %d",
                idx + 1, total, (idx + 1) / total * 100, err_selfies,
            )

    df["selfies"] = selfies_col
    if err_selfies:
        logger.warning("selfies encode failures: %d", err_selfies)

    # Step 3: psmiles → bigsmiles
    bigsmiles_adapter = BigSMILESAdapter()
    bigsmiles_col: list[str | None] = []
    err_bigsmiles = 0

    for idx, psm in enumerate(psmiles_col):
        if psm is None:
            bigsmiles_col.append(None)
            continue
        try:
            bigsmiles_col.append(bigsmiles_adapter.convert(psm, _from="psmiles"))
        except Exception:
            bigsmiles_col.append(None)
            err_bigsmiles += 1

        if (idx + 1) % progress_batch == 0:
            logger.info(
                "[bigsmiles]       %d/%d  (%.1f%%)  err: %d",
                idx + 1, total, (idx + 1) / total * 100, err_bigsmiles,
            )

    df["bigsmiles"] = bigsmiles_col
    if err_bigsmiles:
        logger.warning("bigsmiles failures: %d", err_bigsmiles)

    logger.info("Enrichment complete (%d rows, %d columns)", len(df), len(df.columns))
    return df
