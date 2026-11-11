"""Shared helpers and the reaction-driver interface for chain extension."""

import re

from rdkit import Chem
from rdkit.Chem import AllChem


def get_rxn(smarts: str | list):
    """Compile one or more reaction SMARTS strings.

    Returns a list of ``rdkit.Chem.rdChemReactions.ChemicalReaction`` objects,
    one per input SMARTS.
    """
    rxn_list = []
    if isinstance(smarts, str):
        smarts_list = [smarts]
    else:
        smarts_list = smarts

    if not smarts_list:
        raise ValueError("No SMARTS provided for verification")

    for smarts_str in smarts_list:
        rxn = AllChem.ReactionFromSmarts(str(smarts_str))
        if rxn is None:
            raise ValueError(f"Invalid reaction SMARTS: {smarts_str!r}")
        else:
            rxn_list.append(rxn)

    return rxn_list


def normalize_reactant_wildcards(smarts: str) -> str:
    """Rewrite growing-end wildcards on the reactant side as ``[#0]``.

    In SMARTS, a bare ``*`` is the *universal* atom query and matches any
    atom; against a wildcard-ended chain it would also match ordinary
    backbone bonds.  The extension convention, however, uses ``*`` to mean
    "the dummy atom at the reactive end", which is ``[#0]``.  Only the
    reactant side is rewritten; ``*`` in the product remains a new dummy
    atom.

    Wildcards inside a bracketed atom expression (including recursive SMARTS
    such as ``!$([O,S]C=*)``) are not touched, because they are part of the
    query rather than attachment points.
    """
    parts = str(smarts).split(">")
    if len(parts) < 2:
        return str(smarts)

    reactants = parts[0]
    rewritten = []
    index = 0
    while index < len(reactants):
        char = reactants[index]
        if char == "[":
            end = index
            depth = 0
            while end < len(reactants):
                if reactants[end] == "[":
                    depth += 1
                elif reactants[end] == "]":
                    depth -= 1
                    if depth == 0:
                        break
                end += 1
            bracket = reactants[index : end + 1]
            mapped = re.fullmatch(r"\[\*:(\d+)\]", bracket)
            if bracket == "[*]":
                rewritten.append("[#0]")
            elif mapped:
                rewritten.append(f"[#0:{mapped.group(1)}]")
            else:
                rewritten.append(bracket)
            index = end + 1
        elif char == "*":
            rewritten.append("[#0]")
            index += 1
        else:
            rewritten.append(char)
            index += 1

    parts[0] = "".join(rewritten)
    return ">".join(parts)


def prepare_monomers(monomers) -> list:
    """Normalise monomers to Mol objects and remove duplicates."""
    if isinstance(monomers, (str, Chem.Mol)):
        monomers = [monomers]

    prepared = []
    seen = set()
    for index, monomer in enumerate(monomers):
        if isinstance(monomer, Chem.Mol):
            mol = Chem.Mol(monomer)
        else:
            mol = Chem.MolFromSmiles(str(monomer))
            if mol is None:
                raise ValueError(
                    f"Invalid monomer SMILES at index {index}: {monomer!r}"
                )
        smiles = Chem.MolToSmiles(mol)
        if smiles in seen:
            continue
        seen.add(smiles)
        prepared.append(mol)
    return prepared


class BaseReaction:
    def __init__(self, rxn_type):
        self.rxn_type = rxn_type

    def extend(self, monomers: list, length: int):
        raise NotImplementedError("Child class must implement it")
