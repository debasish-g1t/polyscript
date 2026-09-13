class PolyScriptParser:
    """Parser for a single PolyScript sequence.

    Parses the monomers, reaction SMARTS, and polymer SMILES from a
    PolyScript string.  Results are stored in instance attributes.
    """

    def __init__(self):
        """Initialise an empty PolyScriptParser instance.

        All parse-target attributes are set to ``None`` and the
        sequences list is empty.
        """
        self.monomers = None
        self.reaction = None
        self.polymer = None
        self.errors = []
        # List of (monomers, reaction, polymer) tuples.  For a single-
        # sequence parse this holds one entry; use NestedPolyScriptParser
        # for multi-sequence inputs separated by "|||".
        self.sequences = []

    def get_parts(self, sequence):
        """Split a PolyScript sequence into its three constituent parts.

        Args:
            sequence: A raw PolyScript string, optionally wrapped in ``$``.

        Returns:
            A list of three strings: ``[monomers, reaction, polymer]``.
        """
        if sequence[0] == "$" and sequence[-1] == "$":
            return sequence[1:-1].split("=>")
        else:
            # remove $ from start if exists, and from end if exists
            sequence = sequence.strip("$")
            return sequence.split("=>")

    def parse_monomers(self, monomer_parts):
        """Parse the monomers portion of a PolyScript sequence.

        Stores the extracted monomer SMILES strings in ``self.monomers``.

        Args:
            monomer_parts: The monomers substring, e.g.
                ``{CCO}+{CN}`` (with or without outer ``$``).
        """
        monomers = monomer_parts.split("}+{")
        mon_smiles = []
        for monomer in monomers:
            # Strip leading {
            if monomer.startswith("{"):
                monomer = monomer[1:]
            # Strip trailing }  (separate 'if' — not 'elif' — so single monomers
            # like "{C=C=O}" get BOTH braces removed)
            if monomer.endswith("}"):
                monomer = monomer[:-1]
            mon_smiles.append(monomer)
        self.monomers = mon_smiles

    def _parse_single(self, sequence):
        """Parse a single PolyScript sequence (no ||| splitting).

        Returns a tuple of (monomers, reaction, polymer) or None on failure.
        Errors are accumulated in self.errors.
        """
        parts = []
        try:
            parts = self.get_parts(sequence)
        except Exception:
            self.errors.append("<E-parse-|invalid-parts|>")
            return None
        if len(parts) < 3:
            self.errors.append("<E-parse-|invalid-number-of-parts|>")
            return None
        elif len(parts) > 3:
            self.errors.append("<E-parse-|too-many-parts|>")
            return None

        try:
            self.parse_monomers(parts[0])
        except Exception:
            self.errors.append("<E-parse-|invalid-monomer-parts|>")
            return None
        return (list(self.monomers), parts[1], parts[2])

    def parse(self, sequence):
        """Parse a **single** PolyScript sequence.

        The sequence may be wrapped in ``$…$`` delimiters (they are
        stripped automatically).  Multi-sequence inputs separated by
        ``|||`` are also accepted -- only the **first** sub-sequence is
        parsed (use ``NestedPolyScriptParser`` when you need all of them).

        After a successful parse:
        - ``self.sequences`` is a single-element list of
          (monomers, reaction, polymer).
        - ``self.monomers``, ``self.reaction``, ``self.polymer`` hold the
          parsed values.
        - ``self.errors`` contains parse errors (empty on success).
        """
        self.errors = []
        self.monomers = None
        self.reaction = None
        self.polymer = None

        # Strip optional $…$ wrapper (handles both "$…$" and lone "$" edges)
        inner = sequence
        if sequence.startswith("$") and sequence.endswith("$"):
            inner = sequence[1:-1]
        elif sequence.startswith("$"):
            inner = sequence[1:]
        elif sequence.endswith("$"):
            inner = sequence[:-1]

        # When the input contains "|||", take only the first sub-sequence
        # (NestedPolyScriptParser provides full multi-sequence support).
        if "|||" in inner:
            inner = inner.split("|||")[0].strip()

        result = self._parse_single(inner)
        if result is None:
            self.sequences = []
            return self.errors

        self.sequences = [result]
        self.monomers = list(result[0])
        self.reaction = result[1]
        self.polymer = result[2]
        return self.errors


class NestedPolyScriptParser:
    """Parser for PolyScript inputs that may contain multiple sub-sequences
    separated by ``|||``.

    Each sub-sequence is parsed independently by a ``PolyScriptParser``
    instance.
    """

    def __init__(self):
        """Initialise a NestedPolyScriptParser with empty state."""
        self.errors = []
        self.sequences = []

    def parse(self, sequence):
        self.errors = []
        self.sequences = []

        # Strip the outer $…$ wrapper when it wraps the whole multi-sequence
        inner = sequence
        if sequence.startswith("$") and sequence.endswith("$"):
            inner = sequence[1:-1]

        if "|||" not in inner:
            # Single sequence — delegate to PolyScriptParser
            parser = PolyScriptParser()
            parser.parse(inner)
            self.errors = [parser.errors]
            self.sequences = [(parser.monomers, parser.reaction, parser.polymer)]
            return self.errors

        # Multi-sequence: split on "|||" and parse each independently
        splits = inner.split("|||")
        for split in splits:
            split = split.strip()
            if not split:
                continue
            parser = PolyScriptParser()
            parser.parse(split)
            self.errors.append(parser.errors)
            self.sequences.append((parser.monomers, parser.reaction, parser.polymer))
        return self.errors


# Alias for the reaction-validation use case
SeqParser = PolyScriptParser
