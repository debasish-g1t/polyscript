"""Polymer-type classifier based on SMiPoly functional-group logic.

This module provides ``PolyScriptClassifier``, which classifies monomer
SMILES strings into polymer-relevant functional-group categories (mono-
and poly-functional, as well as olefinic) using SMARTS patterns and
exclusion rules derived from SMiPoly, with additional classes from OMG.
"""

from __future__ import annotations

import pandas as pd
from polyscript.classifier.class_reference import exc_lst, mon_index_map, mon_lst, mon_type_index
from polyscript.utils.executor import BaseExecutor
from polyscript.utils.logger import INFO
from rdkit import Chem


class PolyScriptClassifier(BaseExecutor):
    _log_level = INFO  # class-level default, change via set_log_level()

    def __init__(
        self,
        MinFG: int = 2,
        MaxFG: int = 4,
        include_carbonate: bool = True,
        **kwargs,
    ):
        """Initialize the PolyScript classifier.

        Args:
            MinFG: Minimum number of functional groups for
                poly-functional classification (default 2).
            MaxFG: Maximum number of functional groups for
                poly-functional classification (default 4).
            include_carbonate: If True, append CO and HCHO rows
                to enable carbonate monomer classification.
            **kwargs: Passed to ``BaseExecutor.__init__``.
        """
        super().__init__(**kwargs)
        self.MinFG = MinFG
        self.MaxFG = MaxFG
        self.mon_lst = mon_lst
        self.exc_lst = exc_lst
        self.mon_index_map = mon_index_map
        self.mon_type_index = mon_type_index
        self.include_carbonate = include_carbonate

        self.logger.debug(
            "PolyScriptClassifier initialised (MinFG=%d, MaxFG=%d)", MinFG, MaxFG
        )

    def _extract_FG(self, mol, patt):
        """Count non-overlapping functional-group matches for a pattern.

        Filters out consecutive duplicate matches that differ by a single
        atom shift (symmetric difference of size 2), which typically
        represent the same functional group seen from a different
        orientation rather than distinct occurrences.

        Args:
            mol: An RDKit ``Mol`` object.
            patt: An RDKit ``Mol`` object representing a SMARTS pattern.

        Returns:
            int: Number of distinct functional-group matches.
        """
        numFG = 0
        matchs = mol.GetSubstructMatches(
            patt
        )  # finding the substructures based on the pattern
        if len(matchs) >= 2:  # if higher number of substructures are found
            not_match = []
            for i in range(0, len(matchs) - 1):
                # Checks if the current and next match are NOT extensively overlapping by a
                # single atom shift (i.e., symmetric difference of atom indices is not 2).
                if len(set(matchs[i]) ^ set(matchs[i + 1])) != 2:
                    pass
                else:
                    not_match.append(i + 1)
                numFG = len(
                    [matchs[i] for i in range(0, len(matchs)) if i not in not_match]
                )
        else:
            # if there are only 2 number of matches just add the length as number of functional group
            numFG = len(matchs)
        return numFG

    def _get_mono_FG(self, mol, type_index):
        """Check for mono-functional (single-monomer) classification.

        Tests whether *mol* contains at least one substructure match
        from the monomer SMARTS list for *type_index* (ignoring
        duplicate overlaps).  Also verifies that none of the exclusion
        patterns for that type match.

        Args:
            mol: An RDKit ``Mol`` object, or ``None``.
            type_index: String key into ``self.mon_lst`` and
                ``self.exc_lst`` (e.g. ``"1"`` for vinyl).

        Returns:
            list: ``[bool, int]`` — the bool is ``True`` if the
            monomer type is present with no exclusion violations;
            the int is the total number of raw substructure matches.
        """
        if mol is None:
            return [False, 0]
        else:
            chk_c = 0
            fchk_c = 0
            chk = []
            if len(self.mon_lst[str(type_index)]) != 0:
                for mon in self.mon_lst[str(type_index)]:
                    patt = Chem.MolFromSmarts(mon)
                    if mol.HasSubstructMatch(patt):
                        chk_c = len(mol.GetSubstructMatches(patt))
                        fchk_c = fchk_c + chk_c
                        chk_excl = []
                        for excl in self.exc_lst[str(type_index)]:
                            excl_patt = Chem.MolFromSmarts(excl)
                            if mol.HasSubstructMatch(excl_patt):
                                chk_excl.append(False)
                            else:
                                chk_excl.append(True)
                        if False in chk_excl:
                            chk.append(False)
                        else:
                            chk.append(True)
                    else:
                        chk.append(False)
                if True in chk:
                    fchk = True
                else:
                    fchk = False
            else:
                fchk = False
        return [fchk, fchk_c]

    def _get_poly_FG(self, mol, type_index):
        """Check for poly-functional (multi-monomer) classification.

        Counts all functional-group matches (using ``_extract_FG``) and
        requires the total to fall within ``[MinFG, MaxFG]``.  Also
        validates that none of the exclusion patterns match.

        Args:
            mol: An RDKit ``Mol`` object, or ``None``.
            type_index: String key into ``self.mon_lst`` and
                ``self.exc_lst``.

        Returns:
            list: ``[bool, int]`` — the bool is ``True`` if the
            monomer type is present within the FG count range and
            has no exclusion violations; the int is the count of
            functional-group matches.
        """
        if mol is None:
            return [False, 0]
        chk_c = 0  # check count
        fchk_c = 0  # functional check count
        if len(self.mon_lst[str(type_index)]) != 0:
            for mon in self.mon_lst[str(type_index)]:
                patt = Chem.MolFromSmarts(mon)
                chk_c = self._extract_FG(mol, patt)
                fchk_c = fchk_c + chk_c
            if self.MinFG <= fchk_c <= self.MaxFG:
                chk = []
                for excl in self.exc_lst[str(type_index)]:
                    excl_patt = Chem.MolFromSmarts(excl)
                    if mol.HasSubstructMatch(excl_patt):
                        chk.append(False)
                    else:
                        chk.append(True)
                if False in chk:
                    fchk = False
                else:
                    fchk = True
            else:
                fchk = False
        else:
            fchk = False
        return [fchk, fchk_c]

    def _classify(self, df, ids, mapper_func):
        """Apply a classification mapper to a group of monomer types.

        Iterates over the type indices listed in ``self.mon_index_map[ids]``
        and creates a new DataFrame column for each type name (from
        ``self.mon_type_index``), populated by applying *mapper_func*
        to the ``ROMol`` column.

        Args:
            df: A pandas DataFrame with a ``ROMol`` column.
            ids: Index into ``self.mon_index_map`` (0 for mono, 1 for poly).
            mapper_func: Callable ``(mol, type_index) -> [bool, int]``.
        """
        for i in self.mon_index_map[ids]:
            df[self.mon_type_index[str(i)]] = [
                e[0] for e in df["ROMol"].apply(mapper_func, type_index=i)
            ]

    def classify(self, df, col_name):
        """Run full mono- and poly-functional-group classification.

        Parses SMILES in *df[col_name]*, optionally appends CO / HCHO
        rows for carbonate classification, runs mono-functional
        (``_get_mono_FG``) and poly-functional (``_get_poly_FG``)
        classifiers, and returns the enriched DataFrame.

        Args:
            df: Input DataFrame containing SMILES strings.
            col_name: Column name holding SMILES strings.

        Returns:
            pandas.DataFrame: The input DataFrame with additional
            classification columns (one per monomer type) and the
            ``ROMol`` column dropped.  Returns the original DataFrame
            unchanged if *col_name* is not found.
        """
        # read source file
        DF01 = df
        n_rows = len(DF01)
        self.logger.info(
            "Starting classification (%d rows, column='%s')", n_rows, col_name
        )

        # optionally append CO and HCHO for carbonate classification
        DF02 = DF01.copy()
        if col_name in DF02.columns.to_list():
            if self.include_carbonate:
                # adding two rows for CO and HCHO that will later be used
                # to form carbonate with other monomers
                DFadd = pd.DataFrame(
                    [
                        ["[C-]#[O+]"],  # CO
                        ["C=O"],  # HCHO
                    ],
                    columns=[col_name],
                )
                DF02 = pd.concat([DF02, DFadd], ignore_index=True)
                self.logger.debug(
                    "Appended CO / HCHO rows for carbonate classification"
                )
        else:
            self.logger.error(
                "Column '%s' not found in dataframe (available: %s)",
                col_name,
                DF02.columns.to_list(),
            )
            return df

        # drop NA of smiles, and add chemical structure
        DF02["ROMol"] = DF02[col_name].apply(lambda x: Chem.MolFromSmiles(x))
        DF02[col_name] = DF02["ROMol"].apply(lambda x: Chem.MolToSmiles(x))

        null_count = DF02["ROMol"].isna().sum()
        if null_count:
            self.logger.warning("%d SMILES could not be parsed", null_count)

        # implementing mono func
        self.logger.debug("Running mono-functional-group classification")
        self._classify(DF02, 0, self._get_mono_FG)

        # implementing poly func
        self.logger.debug("Running poly-functional-group classification")
        self._classify(DF02, 1, self._get_poly_FG)

        DF02 = DF02.drop("ROMol", axis=1)
        self.logger.info("Classification complete (%d rows)", len(DF02))
        return DF02
