# Add reference to SMiPoly with deatils

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

# ---------------------------------------------------------------------------
# Module-level SMARTS cache (1.11)
# ---------------------------------------------------------------------------
# SMARTS patterns from monL / mons / excls are compiled to RDKit Mol objects
# once and reused.  Without this, seq_successive alone recompiles 7 patterns
# for every product of every reaction — millions of redundant calls.

_smarts_cache: dict = {}


def _get_smarts(pattern: str):
    """Return a cached RDKit Mol from a SMARTS string.

    Args:
        pattern (str): A SMARTS pattern string.

    Returns:
        rdkit.Chem.Mol: The compiled RDKit molecule object.
    """
    mol = _smarts_cache.get(pattern)
    if mol is None:
        mol = Chem.MolFromSmarts(pattern)
        _smarts_cache[pattern] = mol
    return mol


# Cache for AllChem.ReactionToSmarts (1.13) — keyed by id(reaction)
_rxn_smarts_cache: dict = {}


def _get_rxn_smarts(reaction):
    """Return a cached SMARTS string for an RDKit reaction.

    Args:
        reaction (rdkit.Chem.rdChemReactions.ChemicalReaction): The
            reaction object to serialise.

    Returns:
        str: The SMARTS string representation of the reaction.
    """
    key = id(reaction)
    smarts = _rxn_smarts_cache.get(key)
    if smarts is None:
        smarts = AllChem.ReactionToSmarts(reaction)
        _rxn_smarts_cache[key] = smarts
    return smarts


def genmol(s):
    """
    Generates a molecular object from a SMILES string.

    Args:
        s (str): A SMILES (Simplified Molecular Input Line Entry System) string
            representing the molecular structure.

    Returns:
        rdkit.Chem.Mol or numpy.nan: A molecular object if the SMILES string
        is valid, otherwise returns numpy.nan.

    """
    try:
        m = Chem.MolFromSmiles(s)
    except:
        m = np.nan
    return m


def genc_smi(m):
    """
    Generates a RDkit canonical SMILES string from a molecule object.

    Args:
        m (rdkit.Chem.Mol): A molecule object, from the RDKit library.

    Returns:
        str or np.nan: The SMILES string representation of the molecule if successful,
        otherwise returns np.nan.

    """
    try:
        cS = Chem.MolToSmiles(m)
    except:
        cS = np.nan
    return cS


# count the number of the targetted functional group


def count_fg(m, patt):
    """
    Counts the number of functional groups (FG) in a molecule
    based on a given pattern.

    Args:
        m (rdkit.Chem.Mol): The molecule object to search
            for substructure matches.
        patt (rdkit.Chem.Mol): The pattern molecule used
            to identify substructure matches.

    Returns:
        int: The number of functional groups identified in the molecule.

    """
    numFG = 0
    matchs = m.GetSubstructMatches(patt)
    if len(matchs) >= 2:
        not_match = []
        for i in range(0, len(matchs) - 1):
            if len(set(matchs[i]) ^ set(matchs[i + 1])) != 2:
                pass
            else:
                not_match.append(i + 1)
            numFG = len(
                [matchs[i] for i in range(0, len(matchs)) if i not in not_match]
            )
    else:
        numFG = len(matchs)
    return numFG


# classify candidate compounds for mono-FG monomer


def monomer_sel_mfg(m, mons, excls):
    """
    Determining whether the given small molecule compound
    qualifies as a self-polymerizable monomer and
    categolize it into a monomer class.

    Args:
        m (rdkit.Chem.Mol): The molecule to be analyzed.
            If None or NaN, the function returns default values.
        mons (list of str): A list of SMARTS strings representing
            monomer patterns to match against the molecule.
        excls (list of str): A list of SMARTS strings representing
            exclusion patterns to check against the molecule.

    Returns:
        list: A list containing:

            - fchk (bool): True if the molecule matches any monomer pattern
            and does not match any exclusion pattern, otherwise False.

            - fchk_c (int): The total count of substructure matches
            for all monomer patterns.

    """
    if pd.notna(m):
        chk_c = 0
        fchk_c = 0
        chk = []
        if len(mons) != 0:
            for mon in mons:
                patt = _get_smarts(mon)
                if m.HasSubstructMatch(patt):
                    chk_c = len(m.GetSubstructMatches(patt))
                    fchk_c = fchk_c + chk_c
                    chk_excl = []
                    for excl in excls:
                        excl_patt = _get_smarts(excl)
                        if m.HasSubstructMatch(excl_patt):
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
    else:
        fchk = False
    return [fchk, fchk_c]


def monomer_sel_pfg(m, mons, excls, minFG, maxFG):
    """
    Determining whether the given small molecule compound qualifies as
    a monomer or not. If so, count a number of polymerizeble
    functional group and categolize it into a monomer class.

    Args:
        m (rdkit.Chem.Mol): The monomer molecule to evaluate.
        mons (list of str): A list of SMARTS patterns representing the
            functional groups to count in the monomer.
        excls (list of str): A list of SMARTS patterns representing the
            exclusion patterns to check against the monomer.
        minFG (int): The minimum number of functional groups required.
        maxFG (int): The maximum number of functional groups allowed.

    Returns:
        list: A list containing:

            - fchk (bool): True if the monomer satisfies the conditions,
            False otherwise.

            - fchk_c (int): The total count of functional groups found in
            the monomer.

    """
    if pd.notna(m):
        chk_c = 0
        fchk_c = 0
        if len(mons) != 0:
            for mon in mons:
                patt = _get_smarts(mon)
                chk_c = count_fg(m, patt)
                fchk_c = fchk_c + chk_c
            if minFG <= fchk_c <= maxFG:
                chk = []
                for excl in excls:
                    excl_patt = _get_smarts(excl)
                    if m.HasSubstructMatch(excl_patt):
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
    else:
        fchk = False
    return [fchk, fchk_c]


def run_react_update_pre(prod_P, reaction, pre, p2="none"):
    """Run a single reaction step and update the sequence string.

    Args:
        prod_P (rdkit.Chem.Mol): The reactant molecule.
        reaction (rdkit.Chem.rdChemReactions.ChemicalReaction): The
            reaction to apply.
        pre (str): The existing sequence formation string to extend.
        p2 (str): Identifier for the second reactant (default ``"none"``).

    Returns:
        tuple: ``(new_product, updated_sequence)`` where *new_product* is
        the sanitised product molecule and *updated_sequence* is the
        extended sequence string.
    """
    prods = reaction.RunReactants([prod_P])
    rxn_text = _get_rxn_smarts(reaction)
    new_prod = prods[0][0]
    Chem.SanitizeMol(new_prod)
    pre = handle_sequence_formation(
        pre, genc_smi(prod_P), p2, rxn_text, genc_smi(new_prod)
    )
    return new_prod, pre


# define sequential polymerization for chain polymerization except polyolefine
def seq_chain(prod_P, targ_mon1, Ps_rxnL, mon_dic, monL, pre, max_iter=6):
    """
    This function applied to multifunctional monomers
    for chain polymerization except polyolefine.
    Processes a molecular structure by applying a sequential reactions
    based on specific substructure matches.

    Args:
        prod_P (rdkit.Chem.Mol): The input molecule to be processed.
        targ_mon1 (str): Target monomer type, used to
            determine processing logic.
        Ps_rxnL (dict): A dictionary of polymerization reaction objects
            indexed by integers.
        mon_dic (dict): A dictionary containing monomer class
            (not used in this function).
        monL (list): A list of monomer SMARTS patterns indexed by integers.
        pre (str): The sequence introduction string.

    Returns:
        rdkit.Chem.Mol: The processed molecule after applying the reactions.

    """
    run_count = 0
    if prod_P is not None and prod_P.GetNumAtoms() > 0:
        if targ_mon1 not in ["vinyl", "cOle"]:
            seqFG2 = _get_smarts(monL[[202][0]])
            seqFG3 = _get_smarts(monL[[203][0]])
            seqFG4 = _get_smarts(monL[[204][0]])
            while prod_P.HasSubstructMatch(seqFG2) and run_count < max_iter:
                prod_P, pre = run_react_update_pre(prod_P, Ps_rxnL[202], pre)
                run_count += 1
            while prod_P.HasSubstructMatch(seqFG3) and run_count < max_iter:
                prod_P, pre = run_react_update_pre(prod_P, Ps_rxnL[203], pre)
                run_count += 1
            while prod_P.HasSubstructMatch(seqFG4) and run_count < max_iter:
                prod_P, pre = run_react_update_pre(prod_P, Ps_rxnL[204], pre)
                run_count += 1
        else:
            prod_P = prod_P
    return prod_P, pre


# define sequential polymerization for chain polymerization except polyolefine
def add_fg_mf_mon_chain(prod_P, targ_mon1, monL):
    """Append functional-group marker suffixes to a chain polymer SMILES.

    Detects remaining functional-group substructures (FG2, FG3, FG4)
    on the product and appends them as markers after a ``|||`` separator.
    Skipped for ``"vinyl"`` and ``"cOle"`` monomer types.

    Args:
        prod_P (rdkit.Chem.Mol): The product molecule to inspect.
        targ_mon1 (str): Target monomer type.
        monL (dict): Monomer list indexed by integer keys, containing
            SMARTS patterns.

    Returns:
        str: The SMILES string with appended functional-group markers.
    """
    outro = "|||"
    if prod_P is not None and prod_P.GetNumAtoms() > 0:
        if targ_mon1 not in ["vinyl", "cOle"]:
            seqFG2 = _get_smarts(monL[[202][0]])
            seqFG3 = _get_smarts(monL[[203][0]])
            seqFG4 = _get_smarts(monL[[204][0]])
            if prod_P.HasSubstructMatch(seqFG2):
                outro += "FG2|"
            if prod_P.HasSubstructMatch(seqFG3):
                outro += "FG3|"
            if prod_P.HasSubstructMatch(seqFG4):
                outro += "FG4|"
        else:
            prod_P = prod_P
    prod_p_smi = genc_smi(prod_P)
    return prod_p_smi + outro


# define sequential polymerization for successive polymerization
def seq_successive(prod_P, targ_rxn, monL, Ps_rxnL, P_class, pre, max_iter=10):
    """
    This function applied to multifunctional monomers for
    successive polymerization.
    Processes a molecular structure by applying a sequential reactions
    based on specific substructure matches.

    Args:
        prod_P (rdkit.Chem.Mol): The product molecule to be processed.
        targ_rxn (Any): Target reaction
            (not used in the current implementation).
        monL (list): A list containing SMARTS patterns for functional groups.
        Ps_rxnL (list): A list of polymerization reaction objects
            to be applied to the product molecule.
        P_class (str): The polymer class of the product molecule,
            which determines the reaction sequence.
        pre (str): Pre-text / string that denotes the main formation of product before
            nested functional execution
        max_iter (str): Maximum iteration to go through

    Returns:
        rdkit.Chem.Mol: The processed product molecule
        after applying the reaction sequence.

    Notes:
        - The function uses substructure matching to determine
          which reactions to apply.
        - The behavior of the function depends on the `P_class`
          of the molecule.
        - Specific reaction sequences are applied for classes
          such as 'polyolefin', 'polyoxazolidone', 'polyimide',
          and 'polyester'.
        - If the `P_class` is not recognized,
          the product molecule is returned unchanged.

    """
    counter = 0
    if prod_P is not None and prod_P.GetNumAtoms() > 0:
        seqFG0 = _get_smarts(monL[[200][0]])
        seqFG1 = _get_smarts(monL[[201][0]])
        seqFG2 = _get_smarts(monL[[202][0]])
        seqFG3 = _get_smarts(monL[[203][0]])
        seqFG4 = _get_smarts(monL[[204][0]])
        seqFG5 = _get_smarts(monL[[205][0]])
        seqFG6 = _get_smarts(monL[[206][0]])
        if P_class not in [
            "polyolefin",
            "polyoxazolidone",
        ]:
            # print("[+] Section for not polyolefic ,polyoxazolidone is activated")
            while prod_P.HasSubstructMatch(seqFG1) and counter < max_iter:
                # print(
                #     f"[seq_succesive] [FG1] Processing monomer {Chem.MolToSmiles(prod_P)}"
                # )
                prod_P, pre = run_react_update_pre(prod_P, Ps_rxnL[201], pre)
                counter += 1
            while prod_P.HasSubstructMatch(seqFG2) and counter < max_iter:
                # print(
                #     f"[seq_succesive] [FG2] Processing monomer {Chem.MolToSmiles(prod_P)}"
                # )
                prod_P, pre = run_react_update_pre(prod_P, Ps_rxnL[202], pre)
                counter += 1
            while prod_P.HasSubstructMatch(seqFG3) and counter < max_iter:
                # print(
                #     f"[seq_succesive] [FG3] Processing monomer {Chem.MolToSmiles(prod_P)}"
                # )
                prod_P, pre = run_react_update_pre(prod_P, Ps_rxnL[203], pre)
                counter += 1
            while prod_P.HasSubstructMatch(seqFG4) and counter < max_iter:
                # print(
                #     f"[seq_succesive] [FG4] Processing monomer {Chem.MolToSmiles(prod_P)}"
                # )
                prod_P, pre = run_react_update_pre(prod_P, Ps_rxnL[204], pre)
                counter += 1
            while prod_P.HasSubstructMatch(seqFG5) and counter < max_iter:
                # print(
                #     f"[seq_succesive] [FG5] Processing monomer {Chem.MolToSmiles(prod_P)}"
                # )
                prod_P, pre = run_react_update_pre(prod_P, Ps_rxnL[205], pre)
                counter += 1
            while prod_P.HasSubstructMatch(seqFG6) and counter < max_iter:
                # print(
                #     f"[seq_succesive] [FG6] Processing monomer {Chem.MolToSmiles(prod_P)}"
                # )
                if P_class == "polyimide":
                    prod_P, pre = run_react_update_pre(prod_P, Ps_rxnL[207], pre)
                    counter += 1
                elif P_class == "polyester":
                    prod_P, pre = run_react_update_pre(prod_P, Ps_rxnL[206], pre)
                    counter += 1
                else:
                    counter += 1
        elif P_class in [
            "polyoxazolidone",
        ]:
            # print("[+] Section for polyoxazolidone is activated")
            while prod_P.HasSubstructMatch(seqFG1) and counter < max_iter:
                # print(
                #     f"[seq_succesive] [FG1] [POx] Processing monomer {Chem.MolToSmiles(prod_P)}"
                # )
                prod_P, pre = run_react_update_pre(prod_P, Ps_rxnL[201], pre)
                counter += 1
            while prod_P.HasSubstructMatch(seqFG5) and counter < max_iter:
                # print(
                #     f"[seq_succesive] [FG5] [POx] Processing monomer {Chem.MolToSmiles(prod_P)}"
                # )
                prod_P, pre = run_react_update_pre(prod_P, Ps_rxnL[208], pre)
                counter += 1
        elif P_class in [
            "polyolefin",
        ]:
            # print("[+] Section for polyolefin is activated")
            while prod_P.HasSubstructMatch(seqFG0) and counter < max_iter:
                # print(
                #     f"[seq_succesive] [FG0] [P0] Processing monomer {Chem.MolToSmiles(prod_P)}"
                # )
                prod_P, pre = run_react_update_pre(prod_P, Ps_rxnL[200], pre)
                counter += 1
        else:
            prod_P = prod_P
    return prod_P, pre


# define sequential polymerization for successive polymerization
def add_fg_mf_mon_suc(prod_P, monL, P_class):
    """Append functional-group marker suffixes to a successive polymer SMILES.

    Detects remaining functional-group substructures (FG0–FG6) based on
    the polymer class and appends them as markers after a ``|||`` separator.
    Different polymer classes check different FG subsets.

    Args:
        prod_P (rdkit.Chem.Mol): The product molecule to inspect.
        monL (dict): Monomer list indexed by integer keys, containing
            SMARTS patterns.
        P_class (str): Polymer class label (e.g. ``"polyolefin"``,
            ``"polyoxazolidone"``).

    Returns:
        str: The SMILES string with appended functional-group markers.
    """
    outro = "|||"
    if prod_P is not None and prod_P.GetNumAtoms() > 0:
        seqFG0 = _get_smarts(monL[[200][0]])
        seqFG1 = _get_smarts(monL[[201][0]])
        seqFG2 = _get_smarts(monL[[202][0]])
        seqFG3 = _get_smarts(monL[[203][0]])
        seqFG4 = _get_smarts(monL[[204][0]])
        seqFG5 = _get_smarts(monL[[205][0]])
        seqFG6 = _get_smarts(monL[[206][0]])
        if P_class not in [
            "polyolefin",
            "polyoxazolidone",
        ]:
            # print("[+] Section for not polyolefic ,polyoxazolidone is activated")
            if prod_P.HasSubstructMatch(seqFG1):
                outro += "FG1"
            if prod_P.HasSubstructMatch(seqFG2):
                outro += "FG2"
            if prod_P.HasSubstructMatch(seqFG3):
                outro += "FG3"
            if prod_P.HasSubstructMatch(seqFG4):
                outro += "FG4"
            if prod_P.HasSubstructMatch(seqFG5):
                outro += "FG5"
            if prod_P.HasSubstructMatch(seqFG6):
                outro += "FG6"
        elif P_class in [
            "polyoxazolidone",
        ]:
            # print("[+] Section for polyoxazolidone is activated")
            if prod_P.HasSubstructMatch(seqFG1):
                outro += "FG1"
            if prod_P.HasSubstructMatch(seqFG5):
                outro += "FG5"
        elif P_class in [
            "polyolefin",
        ]:
            # print("[+] Section for polyolefin is activated")
            if prod_P.HasSubstructMatch(seqFG0):
                outro += "FG0"
        else:
            prod_P = prod_P
    smi = genc_smi(prod_P)
    return smi + outro


def handle_sequence_formation(intro, p1, p2, rxn_text, prod_p_smi):
    """Build a PolyScript sequence string for a reaction step.

    Encodes the reaction pathway in the format::

        {reactant1}+{reactant2}=>reaction_smarts=>product_smiles

    If *intro* (previous sequence) is non-empty, the new step is appended
    after a ``|||`` separator.

    Args:
        intro (str): The existing sequence string (may be empty).
        p1 (str): SMILES of the first reactant.
        p2 (str): SMILES of the second reactant.
        rxn_text (str): SMARTS string of the reaction.
        prod_p_smi (str): SMILES of the reaction product.

    Returns:
        str: The complete or extended PolyScript sequence string.
    """
    if intro != "":
        seq = (
            intro
            + "|||"
            + "{"
            + f"{p1}"
            + "}"
            + "+"
            + "{"
            + f"{p2}"
            + "}"
            + "=>"
            + rxn_text
            + "=>"
            + prod_p_smi
        )
    else:
        seq = (
            "{"
            + f"{p1}"
            + "}"
            + "+"
            + "{"
            + f"{p2}"
            + "}"
            + "=>"
            + rxn_text
            + "=>"
            + prod_p_smi
        )
    return seq


# homopolymerization
def homopolymA(mon1, mons, excls, targ_mon1, Ps_rxnL, mon_dic, monL):
    """
    Generates a polymer CRU formed from a single monomer by
    iteratively reacting a monomer until no further reactions are possible.

    Args:
        mon1 (rdkit.Chem.Mol): The initial monomer
            to start the polymerization process.
        mons (list): A list of SMARTS strings representing monomer patterns to
            match against the molecule.
        excls (list): A list of SMARTS strings representing exclusion patterns
            to check against the molecule.
        targ_mon1 (object): The target monomer class for
            the polymerization process.
        Ps_rxnL (list): A dictionary of polymerization reaction objects
            indexed by integers.
        mon_dic (dict):  A dictionary containing monomer class.
        monL (list): A list of monomer SMARTS patterns indexed by integers.

    Returns:
        list: A list of SMILES strings representing the generated homopolymers.

    """
    prod_P = mon1
    prod_Ps = []
    counter = 0
    while monomer_sel_mfg(prod_P, mons, excls)[0] == True:
        prods = Ps_rxnL[mon_dic[targ_mon1]].RunReactants([prod_P])
        for prod_P in prods:
            try:
                Chem.SanitizeMol(prod_P[0])
                prod_P = prod_P[0]
                # prod_P = seq_chain(
                #     prod_P,
                #     targ_mon1=targ_mon1,
                #     Ps_rxnL=Ps_rxnL,
                #     mon_dic=mon_dic,
                #     monL=monL,
                # )
                prod_p_smi = add_fg_mf_mon_chain(
                    prod_P,
                    targ_mon1=targ_mon1,
                    monL=monL,
                )
                prod_Ps.append(prod_p_smi)
            except:
                pass
        counter += 1
    return prod_Ps


# sequence generation from homopolymerization
def seq_homopolymA(mon1, mons, excls, targ_mon1, Ps_rxnL, mon_dic, monL, max_iter=10):
    """
    Generates a polymer CRU formed from a single monomer by
    iteratively reacting a monomer until no further reactions are possible.

    this function returns the exact sequence of reaction from monomer to product unlike `homopolymA` function

    Args:
        mon1 (rdkit.Chem.Mol): The initial monomer
            to start the polymerization process.
        mons (list): A list of SMARTS strings representing monomer patterns to
            match against the molecule.
        excls (list): A list of SMARTS strings representing exclusion patterns
            to check against the molecule.
        targ_mon1 (object): The target monomer class for
            the polymerization process.
        Ps_rxnL (list): A dictionary of polymerization reaction objects
            indexed by integers.
        mon_dic (dict):  A dictionary containing monomer class.
        monL (list): A list of monomer SMARTS patterns indexed by integers.
        max_iter (int): The maximum number of iterations to perform to avoid infinite loops.

    Returns:
        list: A list of sequence strings representing the entire pathway as specific syntax.

    """
    prod_P = mon1
    returnable = {}
    counter = 0
    rxn_key = mon_dic[targ_mon1]
    if rxn_key not in Ps_rxnL:
        return []
    rxn = Ps_rxnL[rxn_key]
    rxn_text = _get_rxn_smarts(rxn)

    while monomer_sel_mfg(prod_P, mons, excls)[0] is True and counter < max_iter:
        intro = returnable.get(prod_P, "")
        prev_smi = genc_smi(prod_P)
        prods = rxn.RunReactants([prod_P])
        for result in prods:
            try:
                Chem.SanitizeMol(result[0])
                new_prod = result[0]
                seq = handle_sequence_formation(
                    intro, prev_smi, "none", rxn_text, genc_smi(new_prod)
                )
                new_prod, seq = seq_chain(
                    new_prod,
                    targ_mon1=targ_mon1,
                    Ps_rxnL=Ps_rxnL,
                    mon_dic=mon_dic,
                    monL=monL,
                    pre=seq,
                )
                returnable[new_prod] = seq
                prod_P = new_prod
            except Exception:
                pass
        counter += 1
    return [value for value in returnable.values() if value is not None]


# binarypolymerization
def bipolymA(reactant, targ_rxn, monL, Ps_rxnL, P_class):
    """
    Generates a polymer CRU formed from two monomers by iteratively reacting
    a monomer until no further reactions are possible.

    Args:
        reactant (tuple): A tuple of reactant molecules
            to be used in the reaction.
        targ_rxn (rdkit.Chem.rdChemReactions.ChemicalReaction):
            The target chemical reaction to apply.
        monL (list):  A list of monomer SMARTS patterns indexed by integers.
        Ps_rxnL (dict): A list of monomer SMARTS patterns indexed by integers.
        P_class (type): A class type used for polymer processing.

    Returns:
        list: A list of SMILES strings representing
        the generated polymer products.

    """
    prods = targ_rxn.RunReactants(reactant)
    prod_Ps = []
    for prod_P in prods:
        try:
            Chem.SanitizeMol(prod_P[0])
            prod_P = prod_P[0]
            # print("[+] copolymer sequence generation started")
            # prod_P = seq_successive(
            #     prod_P,
            #     targ_rxn=targ_rxn,
            #     monL=monL,
            #     Ps_rxnL=Ps_rxnL,
            #     P_class=P_class,
            # )
            # print("[+] copolymer sequence generation completed")
            prod_p_smi = genc_smi(prod_P)
            prod_Ps.append(prod_p_smi)
        except:
            # print("[+] Skipping the sequence generation due to timeout")
            pass
    return prod_Ps


# sequence generation from binarypolymerization
def seq_bipolymA(reactant, targ_rxn, monL, Ps_rxnL, P_class):
    """
    Generates a polymer CRU formed from two monomers by iteratively reacting
    a monomer until no further reactions are possible.

    Args:
        reactant (tuple): A tuple of reactant molecules
            to be used in the reaction.
        targ_rxn (rdkit.Chem.rdChemReactions.ChemicalReaction):
            The target chemical reaction to apply.
        monL (list):  A list of monomer SMARTS patterns indexed by integers.
        Ps_rxnL (dict): A list of monomer SMARTS patterns indexed by integers.
        P_class (type): A class type used for polymer processing.

    Returns:
        list: A list of SMILES strings representing
        the generated polymer products.

    """
    prods = targ_rxn.RunReactants(reactant)
    rxn_text = _get_rxn_smarts(targ_rxn)
    prod_Ps = []
    for prod_P in prods:
        try:
            Chem.SanitizeMol(prod_P[0])
            prod_P = prod_P[0]
            prod_p_smi = handle_sequence_formation(
                "",
                genc_smi(reactant[0]),
                genc_smi(reactant[1]),
                rxn_text,
                genc_smi(prod_P),
            )
            prod_Ps.append(prod_p_smi)
            prod_P, pre = seq_successive(
                prod_P,
                targ_rxn=targ_rxn,
                monL=monL,
                Ps_rxnL=Ps_rxnL,
                P_class=P_class,
                pre=prod_p_smi,
            )
            if pre not in prod_Ps:
                prod_Ps.append(pre)
        except:
            pass

    return prod_Ps
