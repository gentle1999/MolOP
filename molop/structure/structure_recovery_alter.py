import itertools
import time
from random import shuffle
from typing import List, Tuple, Union

from openbabel import openbabel as ob
from openbabel import pybel
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds

from molop.config import molopconfig
from molop.logger.logger import moloplogger
from molop.structure.structure import (
    bond_list,
    get_bond_pairs,
    get_radical_number,
    get_under_bonded_number,
    omol_to_rdmol_by_graph,
    rdmol_to_omol,
    reset_atom_index,
)
from molop.utils.consts import (
    get_possible_metal_radicals,
    metal_valence_avialable_minor,
    metal_valence_avialable_prior,
)
from molop.utils.functions import is_metal

pt = Chem.GetPeriodicTable()
HETEROATOM = (9, 8, 17, 7, 35, 54, 16, 34, 15)
DEBUG_TAG = "[STRUCTURE RECOVERY]"


def has_bridge(omol: pybel.Molecule, atom_idx_1: int, atom_idx_2: int) -> bool:
    atom_1 = omol.OBMol.GetAtom(atom_idx_1)
    atom_2 = omol.OBMol.GetAtom(atom_idx_2)
    neighbors_1 = [atom.GetIdx() for atom in ob.OBAtomAtomIter(atom_1)]
    neighbors_2 = [atom.GetIdx() for atom in ob.OBAtomAtomIter(atom_2)]
    for neighbor_1 in neighbors_1:
        if neighbor_1 in neighbors_2:
            return True
    return False


def make_connections(omol: pybel.Molecule, factor: float = 1.45) -> pybel.Molecule:
    """
    Make connections between atoms in the pybel molecule.
    Parameters:
        omol (pybel.Molecule): The pybel molecule to be connected.
    Returns:
        pybel.Molecule: The connected pybel molecule.
    """
    donate_smarts = pybel.Smarts(
        "[Nv0,Cv1,Nv3,Clv1,Clv2,Clv3,Brv1,Brv2,Brv3,Iv1,Iv2,Iv3]"
    )
    accept_smarts = pybel.Smarts(
        "[Hv0,Bv2,Bv3,Cv0,Cv1,Cv2,Cv3,Nv1,Nv2,Ov0,Ov1,Clv0,Siv3,Pv2,Sv0,Sv1,Brv0,Iv0]"
    )
    donate_atoms = list(itertools.chain(*donate_smarts.findall(omol)))
    accept_atoms = list(itertools.chain(*accept_smarts.findall(omol)))

    for atom in donate_atoms:
        pairs = [(atom, accept_atom) for accept_atom in accept_atoms]
        pairs = sorted(
            pairs,
            key=lambda x: omol.OBMol.GetAtom(x[0]).GetDistance(
                omol.OBMol.GetAtom(x[1])
            ),
        )
        for pair in pairs:
            donate_atom = omol.OBMol.GetAtom(pair[0])
            accept_atom = omol.OBMol.GetAtom(pair[1])
            distance = donate_atom.GetDistance(accept_atom)
            if (
                distance
                < (
                    pt.GetRcovalent(donate_atom.GetAtomicNum())
                    + pt.GetRcovalent(accept_atom.GetAtomicNum())
                )
                * factor
                and pair[0] in donate_atoms
                and pair[1] in accept_atoms
                and pair[0] != pair[1]
                # and not has_bridge(omol, pair[0], pair[1])
            ):
                if omol.OBMol.GetBond(pair[0], pair[1]) is None:
                    omol.OBMol.AddBond(pair[0], pair[1], 1)
                    moloplogger.debug(
                        f"{DEBUG_TAG} | Add bond {pair[0]} - {pair[1]}, distance {distance}"
                    )
                    continue
                if omol.OBMol.GetBond(pair[0], pair[1]).GetBondOrder() == 0:
                    omol.OBMol.GetBond(pair[0], pair[1]).SetBondOrder(1)
                    moloplogger.debug(
                        f"{DEBUG_TAG} | Set bond order {pair[0]} - {pair[1]} to 1"
                    )

            donate_atoms = list(itertools.chain(*donate_smarts.findall(omol)))
            accept_atoms = list(itertools.chain(*accept_smarts.findall(omol)))

    return omol


def pre_clean(omol: pybel.Molecule) -> pybel.Molecule:
    """
    Clean the pybel molecule before structure recovery.
    This function suppose no formal charge in the molecule.

    Parameters:
        omol (pybel.Molecule): The pybel molecule to be cleaned.
    Returns:
        pybel.Molecule: The cleaned pybel molecule.
    """
    smarts = pybel.Smarts("[Cv5,Nv5,Pv5]=,#[*]")
    while res := smarts.findall(omol):
        idxs = res.pop(0)
        omol.OBMol.GetBond(idxs[0], idxs[1]).SetBondOrder(
            omol.OBMol.GetBond(idxs[0], idxs[1]).GetBondOrder() - 1
        )

    smarts = pybel.Smarts("[#6]1([#6]2)([#6]3)[#7]23[#6]1")
    while res := smarts.findall(omol):
        idxs = res.pop(0)
        for idx in idxs:
            indexs = set(idxs) - set([idx])
            if all(omol.OBMol.GetBond(idx, idx_2) for idx_2 in indexs):
                if omol.OBMol.GetAtom(idx).GetAtomicNum() == 7:
                    bcp_n = idx
                if omol.OBMol.GetAtom(idx).GetAtomicNum() == 6:
                    bcp_c = idx
        omol.OBMol.DeleteBond(omol.OBMol.GetBond(bcp_n, bcp_c))
        moloplogger.debug(f"{DEBUG_TAG} | Fix N-BCP: {bcp_n} - {bcp_c}")

    smarts = pybel.Smarts("[#6]1([#6]2)[#7]2[#6]1")
    while res := smarts.findall(omol):
        idxs = res.pop(0)
        for idx in idxs:
            indexs = set(idxs) - set([idx])
            if all(omol.OBMol.GetBond(idx, idx_2) for idx_2 in indexs):
                if omol.OBMol.GetAtom(idx).GetAtomicNum() == 7:
                    amine_n = idx
                if omol.OBMol.GetAtom(idx).GetAtomicNum() == 6:
                    butyl_c = idx
        omol.OBMol.DeleteBond(omol.OBMol.GetBond(amine_n, butyl_c))
        moloplogger.debug(f"{DEBUG_TAG} | Fix N-BCP: {amine_n} - {butyl_c}")

    smarts = pybel.Smarts("[Siv5]-[O,F]")
    while res := smarts.findall(omol):
        idxs = res.pop(0)
        omol.OBMol.DeleteBond(omol.OBMol.GetBond(idxs[0], idxs[1]))
        moloplogger.debug(
            f"{DEBUG_TAG} | Remove Over bonding Si: {idxs[0]} - {idxs[1]}"
        )

    return omol


def final_check_omol(
    omol: pybel.Molecule, total_charge: int, total_radical_electrons: int
) -> bool:
    """
    Check if the final structure is valid.
    Parameters:
        omol (pybel.Molecule): The pybel molecule to be checked.
        total_charge (int): The total charge of the molecule.
        total_radical_electrons (int): The total radical of the molecule.
    Returns:
        bool: True if the final structure is valid, False otherwise.
    """
    if sum(atom.OBAtom.GetFormalCharge() for atom in omol.atoms) != total_charge:
        moloplogger.debug(
            f"{DEBUG_TAG} | Charge check failed, total charge: "
            f"{total_charge}, actual charge: {sum(atom.OBAtom.GetFormalCharge() for atom in omol.atoms)}"
        )
        return False
    if (
        sum(get_radical_number(atom.OBAtom) for atom in omol.atoms)
        != total_radical_electrons
    ):
        moloplogger.debug(
            f"{DEBUG_TAG} | Radical check failed, total radical: "
            f"{total_radical_electrons}, actual radical: {sum(get_radical_number(atom.OBAtom) for atom in omol.atoms)}"
        )
        return False
    return True


def final_check(
    rwmol: Chem.rdchem.RWMol, total_charge: int, total_radical_electrons: int
) -> bool:
    """
    Check if the final structure is valid.
    Parameters:
        rwmol (Chem.rdchem.RWMol): The editable rdkit molecule.
        total_charge (int): The total charge of the molecule.
        total_radical_electrons (int): The total radical of the molecule.
    Returns:
        bool: True if the final structure is valid, False otherwise.
    """
    if sum(atom.GetFormalCharge() for atom in rwmol.GetAtoms()) != total_charge:
        moloplogger.debug(f"{DEBUG_TAG} | Charge check failed")
        return False
    if (
        sum(atom.GetNumRadicalElectrons() for atom in rwmol.GetAtoms())
        != total_radical_electrons
    ):
        moloplogger.debug(f"{DEBUG_TAG} | Radical check failed")
        return False
    moloplogger.debug(f"{DEBUG_TAG} | Final check passed")
    return True


def omol_score(omol_tuple: Tuple[pybel.Molecule, int, int]) -> int:
    """
    Calculate the structural recovery score of a molecule, the lower the score the better
    the structural recovery.
    The criteria are 2 points for each unbonded electron and 1 point for each absolute
    value of charge.

    This function can only be used for comparison between isomers recovered from the same
    set of atomic coordinates.

    Parameters:
        omol_tuple (Tuple[pybel.Molecule, int, int]): A tuple containing a pybel.Molecule
            object and its charge and radical number.

    Returns:
        int: The structural recovery score of a molecule.
    """
    score = 0
    score += 2 * abs(
        sum(abs(get_radical_number(atom.OBAtom)) for atom in omol_tuple[0].atoms)
        - omol_tuple[2]
    )
    score += sum(abs(atom.OBAtom.GetFormalCharge()) for atom in omol_tuple[0].atoms)
    return score


def structure_score(rwmol: Chem.rdchem.RWMol) -> float:
    """
    Calculate the score of the given rdkit molecule.

    The score function is designed to rank the possible structures of the given rdkit molecule.
    The score is calculated by the sum of the total valence and the number of formal charge.
    The higher the score, the better the structure.

    Parameters:
        rwmol (Chem.rdchem.RWMol): The rdkit molecule to be scored.
    Returns:
        float: The score of the given rdkit molecule.
    """
    # the more total valence the better
    Chem.SanitizeMol(rwmol)
    total_valence = sum(atom.GetTotalValence() for atom in rwmol.GetAtoms())
    # the less formal charge the better
    total_formal_charge = sum(abs(atom.GetFormalCharge()) for atom in rwmol.GetAtoms())
    total_num_radical = sum(atom.GetNumRadicalElectrons() for atom in rwmol.GetAtoms())
    num_fragments = len(Chem.GetMolFrags(rwmol))
    metal_atoms = [atom for atom in rwmol.GetAtoms() if is_metal(atom.GetAtomicNum())]
    metal_score = 0
    for atom in metal_atoms:
        if atom.GetFormalCharge() in metal_valence_avialable_prior[atom.GetSymbol()]:
            metal_score += 50
        if atom.GetFormalCharge() < 0:
            metal_score += atom.GetFormalCharge() * 100

    score = (
        10 * total_valence
        + 10 / (total_formal_charge + 1)
        + 20 / (total_num_radical + 1)
        + 50 / (num_fragments + 1)
        + metal_score
    )
    """moloplogger.info(
        f"{Chem.CanonSmiles(Chem.MolToSmiles(rwmol))} score: {score}, "
        f"total valence: {total_valence}, total formal charge: {total_formal_charge},"
        f"total radical: {total_num_radical}, num fragments: {num_fragments}, metal score: {metal_score}"
    )"""
    return score


def clean_carbine_neighbor_unsaturated(omol: pybel.Molecule) -> pybel.Molecule:
    """
    clean the carbine neighbor unsaturated.
    Parameters:
        omol (pybel.Molecule): The molecule to be transformed.
    """
    smarts = pybel.Smarts("[*]-[*]=[*]")
    res = list(smarts.findall(omol))
    while len(res):
        idxs = res.pop(0)
        atom1 = omol.OBMol.GetAtom(idxs[0])
        atom2 = omol.OBMol.GetAtom(idxs[1])
        atom3 = omol.OBMol.GetAtom(idxs[2])
        if get_radical_number(atom1) == 2 and get_radical_number(atom3) == 0:
            moloplogger.debug(
                f"{DEBUG_TAG} | fixing carbine neighbor: {atom1.GetIdx()} "
                f"and {atom2.GetIdx()} {atom3.GetIdx()}"
            )
            omol.OBMol.GetBond(idxs[1], idxs[2]).SetBondOrder(
                int(omol.OBMol.GetBond(idxs[1], idxs[2]).GetBondOrder() - 1)
            )
            omol.OBMol.GetBond(idxs[0], idxs[1]).SetBondOrder(
                int(omol.OBMol.GetBond(idxs[0], idxs[1]).GetBondOrder() + 1)
            )
    return omol


def eliminate_high_positive_charge_atoms(
    omol: pybel.Molecule, given_charge: int
) -> Tuple[pybel.Molecule, int]:
    smarts = pybel.Smarts("[*+1,*+2,*+3]-[Ov1+0,Nv2+0,Sv1+0]")
    while res := smarts.findall(omol):
        idxs = res.pop(0)
        omol.OBMol.GetAtom(idxs[1]).SetFormalCharge(-1)
        moloplogger.debug(f"{DEBUG_TAG} | Eliminate atom {idxs[1]} with -1 charge")
        given_charge += 1
    return omol, given_charge


def eliminate_CN_in_doubt(
    omol: pybel.Molecule, given_charge: int
) -> Tuple[pybel.Molecule, int]:
    smarts = pybel.Smarts("[#6v4+0]=,#[#7v4+1,#15v4+1]")
    doubt_pair = smarts.findall(omol)
    CN_in_doubt = len(doubt_pair)
    if CN_in_doubt > 0:
        moloplogger.debug(
            f"{DEBUG_TAG} | Fixing CN in doubt, number of pairs: {CN_in_doubt}, "
            f"charge to be allocated: {given_charge}"
        )
    if CN_in_doubt % 2 == 0 and CN_in_doubt > 0:
        for atom_1, atom_2 in doubt_pair[: CN_in_doubt // 2]:
            omol.OBMol.GetAtom(atom_1).SetFormalCharge(-1)
            omol.OBMol.GetBond(atom_1, atom_2).SetBondOrder(
                omol.OBMol.GetBond(atom_1, atom_2).GetBondOrder() - 1
            )
            omol.OBMol.GetAtom(atom_2).SetFormalCharge(0)
            given_charge += 2
            moloplogger.debug(
                f"{DEBUG_TAG} | Fix CN in doubt: {atom_1} - {atom_2}, "
                f"charge to be allocated: {given_charge}"
            )
    return omol, given_charge


def eliminate_carboxyl(
    omol: pybel.Molecule, given_charge: int
) -> Tuple[pybel.Molecule, int]:
    smarts = pybel.Smarts("[Ov1+0]-C=O")
    res = list(smarts.findall(omol))
    while len(res):
        idxs = res.pop(0)
        omol.OBMol.GetAtom(idxs[0]).SetFormalCharge(
            int(omol.OBMol.GetAtom(idxs[0]).GetFormalCharge() - 1)
        )
        given_charge += 1
        moloplogger.debug(f"{DEBUG_TAG} | Eliminate atom {idxs[0]} with -1 charge")
    return omol, given_charge


def eliminate_carbine_neighbor_heteroatom(
    omol: pybel.Molecule, given_charge: int
) -> Tuple[pybel.Molecule, int]:
    """
    Fix the carbine neighbor heteroatom.
    Parameters:
        omol (pybel.Molecule): The molecule to be transformed.
        given_charge (int): The given charge.
    Returns:
        int: The remaining charge to be allocated.
    """
    for atom in omol.atoms:
        if get_radical_number(atom.OBAtom) == 2:
            for neibhor in ob.OBAtomAtomIter(atom.OBAtom):
                if get_radical_number(neibhor):
                    return omol, given_charge
            for neighbor in ob.OBAtomAtomIter(atom.OBAtom):
                if (
                    neighbor.GetAtomicNum() in HETEROATOM
                    and neighbor.GetFormalCharge() == 0
                    and get_radical_number(neighbor) == 0
                ):
                    atom.OBAtom.GetBond(neighbor).SetBondOrder(
                        atom.OBAtom.GetBond(neighbor).GetBondOrder() + 1
                    )
                    atom.OBAtom.SetFormalCharge(atom.OBAtom.GetFormalCharge() - 1)
                    neighbor.SetFormalCharge(1)
                    moloplogger.debug(
                        f"{DEBUG_TAG} | eliminating carbine neighbor: {atom.OBAtom.GetIdx()} and "
                        f"{neighbor.GetIdx()}, charge to be allocated: {given_charge}"
                    )
                    break
    return omol, given_charge


def eliminate_NNN(
    omol: pybel.Molecule, given_charge: int
) -> Tuple[pybel.Molecule, int]:
    smarts = pybel.Smarts("[#7v1+0]-[#7v2+0]-[#7v1+0]")
    while res := smarts.findall(omol):
        idxs = res.pop(0)
        atom1 = omol.OBMol.GetAtom(idxs[0])
        atom2 = omol.OBMol.GetAtom(idxs[1])
        atom3 = omol.OBMol.GetAtom(idxs[2])
        omol.OBMol.GetBond(idxs[0], idxs[1]).SetBondOrder(
            int(omol.OBMol.GetBond(idxs[0], idxs[1]).GetBondOrder() + 1)
        )
        omol.OBMol.GetBond(idxs[1], idxs[2]).SetBondOrder(
            int(omol.OBMol.GetBond(idxs[1], idxs[2]).GetBondOrder() + 1)
        )
        atom1.SetFormalCharge(atom1.GetFormalCharge() - 1)
        atom2.SetFormalCharge(atom2.GetFormalCharge() + 1)
        atom3.SetFormalCharge(atom3.GetFormalCharge() - 1)
        given_charge += 1
        moloplogger.debug(f"{DEBUG_TAG} | Fix NNN-: {idxs[0]}, {idxs[1]}, {idxs[2]}")
    return omol, given_charge


def break_one_bond(
    omol: pybel.Molecule, given_charge: int = 0, given_radical: int = 0
) -> pybel.Molecule:
    """
    Breaks one chemical bond under specific conditions.

    This function checks if there are any unsatisfied valence states in the molecule,
    and based on the charge and radical conditions, finds and breaks an appropriate bond.

    Parameters:
        omol (pybel.Molecule): The molecule to break bonds in.
        given_charge (int, optional): The allowed charge. Defaults to 0.
        given_radical (int, optional): The given radical. Defaults to 0.
    """
    # Check if there are any unsatisfied valence states and if there is an allowed charge
    # or given radical
    smarts = pybel.Smarts("[*+0]#,=[*+0]")
    # Loop to find suitable bonds
    while res := smarts.findall(omol):
        if (
            sum(get_radical_number(atom) for atom in ob.OBMolAtomIter(omol.OBMol))
            >= abs(given_charge) + given_radical
        ):
            return omol
        # Get the indices of the first suitable bond's atoms
        idxs = res.pop(0)
        # Log debug information indicating the bond to be broken
        moloplogger.debug(f"{DEBUG_TAG} | break bond {idxs[0]} and {idxs[1]}")
        # Get and reduce the bond order of the found bond
        bond = omol.OBMol.GetBond(idxs[0], idxs[1])
        bond.SetBondOrder(bond.GetBondOrder() - 1)

    if all(bond.GetBondOrder() == 1 for bond in ob.OBMolBondIter(omol.OBMol)):
        for bond in ob.OBMolBondIter(omol.OBMol):
            if (
                sum(get_radical_number(atom) for atom in ob.OBMolAtomIter(omol.OBMol))
                >= abs(given_charge) + given_radical
            ):
                return omol
            omol.OBMol.DeleteBond(bond)
    return omol


def clean_neighbor_radicals(omol: pybel.Molecule) -> pybel.Molecule:
    """
    Clean the neighbor radicals.
    Parameters:
        omol (pybel.Molecule): The molecule to be transformed.
    """
    for bond in ob.OBMolBondIter(omol.OBMol):
        if get_radical_number(bond.GetBeginAtom()) and get_radical_number(
            bond.GetEndAtom()
        ):
            bond.SetBondOrder(bond.GetBondOrder() + 1)
            moloplogger.debug(
                f"{DEBUG_TAG} | Cleaning bond {bond.GetBeginAtom().GetIdx()} "
                f"and {bond.GetEndAtom().GetIdx()}"
            )
    return omol


def eliminate_charge_spliting(
    omol: pybel.Molecule, given_charge: int
) -> Tuple[pybel.Molecule, int]:
    if (
        all(atom.GetFormalCharge() == 0 for atom in ob.OBMolAtomIter(omol.OBMol))
        and sum(get_radical_number(atom) for atom in ob.OBMolAtomIter(omol.OBMol)) >= 2
    ):
        radical_atoms = [
            atom for atom in ob.OBMolAtomIter(omol.OBMol) if get_radical_number(atom)
        ]
        moloplogger.debug(
            f"{DEBUG_TAG} | Eliminate charge splitting, number of radical atoms: {len(radical_atoms)}"
        )
        while len(radical_atoms) > abs(given_charge) + 1:
            for atom in radical_atoms:
                if atom.GetAtomicNum() == 8:
                    atom.SetFormalCharge(atom.GetFormalCharge() - 1)
                    given_charge += 1
                    radical_atoms.remove(atom)
                    break
            else:
                break
        while len(radical_atoms) > abs(given_charge) + 1:
            for atom in radical_atoms:
                if atom.GetAtomicNum() == 7:
                    atom.SetFormalCharge(atom.GetFormalCharge() - 1)
                    given_charge += 1
                    radical_atoms.remove(atom)
                    break
            else:
                break
        while len(radical_atoms) > abs(given_charge) + 1:
            for atom in radical_atoms:
                if atom.GetAtomicNum() == 6 and not any(
                    _atom
                    for _atom in ob.OBAtomAtomIter(atom)
                    if _atom.GetAtomicNum() in HETEROATOM
                ):
                    atom.SetFormalCharge(atom.GetFormalCharge() - 1)
                    given_charge += 1
                    radical_atoms.remove(atom)
                    break
            else:
                break
        while len(radical_atoms) > abs(given_charge) + 1:
            for atom in radical_atoms:
                if atom.GetAtomicNum() == 6:
                    atom.SetFormalCharge(atom.GetFormalCharge() - 1)
                    given_charge += 1
                    radical_atoms.remove(atom)
                    break
            else:
                break
    return omol, given_charge


def get_one_step_resonance(omol: pybel.Molecule) -> List[pybel.Molecule]:
    """
    Get one step resonance structures from a molecule.

    Parameters:
        omol (pybel.Molecule): The molecule to be processed.

    Returns:
        List[pybel.Molecule]: A list of one step resonance structures.
    """
    # Get a list of atoms with unpaired electrons (radicals)
    smarts = pybel.Smarts("[*]-,=[*]=,#[*]")
    res = list(smarts.findall(omol))
    result = []
    for idxs in res:
        if get_radical_number(omol.OBMol.GetAtom(idxs[0])) >= 1:
            new_omol = rdmol_to_omol(omol_to_rdmol_by_graph(omol))
            new_omol.OBMol.GetBond(idxs[1], idxs[2]).SetBondOrder(
                new_omol.OBMol.GetBond(idxs[1], idxs[2]).GetBondOrder() - 1
            )
            new_omol.OBMol.GetBond(idxs[0], idxs[1]).SetBondOrder(
                new_omol.OBMol.GetBond(idxs[0], idxs[1]).GetBondOrder() + 1
            )

            result.append(new_omol)
    return result


def get_radical_resonances(omol: pybel.Molecule) -> List[pybel.Molecule]:
    """
    Retrieves a list of molecular structures that exhibit radical resonance.

    Parameters:
        omol (pybel.Molecule): The input molecule object

    Returns:
        List[pybel.Molecule]: A list of molecules with radical resonance structures
    """

    # Initialize the resonance molecules list with the input molecule
    resonances = {omol.write("smi"): omol}

    new_resonances = get_one_step_resonance(omol)
    for new_resonance in new_resonances:
        resonances[new_resonance.write("smi")] = new_resonance

    for temp_omol in new_resonances:
        new_new_resonances = get_one_step_resonance(temp_omol)
        for new_resonance in new_new_resonances:
            resonances[new_resonance.write("smi")] = new_resonance

    # Return the list of resonance molecules
    return list(resonances.values())


def eliminate_1_3_dipole(
    omol: pybel.Molecule, given_charge: int
) -> Tuple[pybel.Molecule, int]:
    smarts = pybel.Smarts("[*-1]-,=[N+0,O+0]-,=[*]")
    res = list(smarts.findall(omol))
    while len(res):
        idxs = res.pop(0)
        atom1 = omol.OBMol.GetAtom(idxs[0])
        atom2 = omol.OBMol.GetAtom(idxs[1])
        atom3 = omol.OBMol.GetAtom(idxs[2])
        if (
            get_radical_number(atom3)
            and pt.GetNOuterElecs(atom2.GetAtomicNum()) + atom2.GetTotalValence() == 8
        ):
            moloplogger.debug(
                f"{DEBUG_TAG} | Eliminate 1,3 dipole: {atom1.GetIdx()} "
                f"{atom2.GetIdx()} {atom3.GetIdx()}"
            )
            atom2.SetFormalCharge(atom2.GetFormalCharge() + 1)
            omol.OBMol.GetBond(idxs[1], idxs[2]).SetBondOrder(
                int(omol.OBMol.GetBond(idxs[1], idxs[2]).GetBondOrder() + 1)
            )
            given_charge -= 1
    return omol, given_charge


def eliminate_positive_charges(
    omol: pybel.Molecule, given_charge: int
) -> Tuple[pybel.Molecule, int]:
    smarts = pybel.Smarts("[#6v3+0,#6v2+0,#1v0+0]")
    while given_charge > 0 and (res := smarts.findall(omol)):
        idxs = res.pop(0)
        omol.OBMol.GetAtom(idxs[0]).SetFormalCharge(1)
        given_charge -= 1
    for atom in omol.atoms:
        if given_charge <= 0:
            break
        if get_radical_number(atom.OBAtom) >= 1 and atom.OBAtom.GetFormalCharge() == 0:
            to_add = min(get_radical_number(atom.OBAtom), given_charge)
            atom.OBAtom.SetFormalCharge(to_add)
            given_charge -= to_add
    return omol, given_charge


def eliminate_negative_charges(
    omol: pybel.Molecule, given_charge: int
) -> Tuple[pybel.Molecule, int]:
    for atom in omol.atoms:
        if given_charge >= 0:
            break
        if (
            atom.atomicnum in HETEROATOM
            and atom.OBAtom.GetFormalCharge() == 0
            and get_radical_number(atom.OBAtom)
        ):
            to_add = min(get_radical_number(atom.OBAtom), abs(given_charge))
            atom.OBAtom.SetFormalCharge(-to_add)
            given_charge += to_add
            moloplogger.debug(
                f"{DEBUG_TAG} | Eliminate negative charge: {atom.OBAtom.GetIdx()} "
                f"with charge {to_add}"
            )

    # Step 2.3.2: Try to find the heteroatom with negative charge first.
    smarts = pybel.Smarts("[#6v3+0]")
    while given_charge < 0 and (res := smarts.findall(omol)):
        idxs = res.pop(0)
        to_add = min(get_radical_number(omol.OBMol.GetAtom(idxs[0])), abs(given_charge))
        omol.OBMol.GetAtom(idxs[0]).SetFormalCharge(-to_add)
        given_charge += to_add
        moloplogger.debug(
            f"{DEBUG_TAG} | Eliminate negative charge: {atom.OBAtom.GetIdx()} "
            f"with charge {to_add}"
        )
    smarts = pybel.Smarts("[#1v0+0]")
    while given_charge < 0 and (res := smarts.findall(omol)):
        idxs = res.pop(0)
        to_add = min(get_radical_number(omol.OBMol.GetAtom(idxs[0])), abs(given_charge))
        omol.OBMol.GetAtom(idxs[0]).SetFormalCharge(-to_add)
        given_charge += to_add
        moloplogger.debug(
            f"{DEBUG_TAG} | Eliminate negative charge: {atom.OBAtom.GetIdx()} "
            f"with charge {to_add}"
        )
    smarts = pybel.Smarts("[#6v2+0,#6v1+0,#6v0+0]")
    while given_charge < 0 and (res := smarts.findall(omol)):
        idxs = res.pop(0)
        to_add = min(get_radical_number(omol.OBMol.GetAtom(idxs[0])), abs(given_charge))
        omol.OBMol.GetAtom(idxs[0]).SetFormalCharge(-to_add)
        given_charge += to_add
        moloplogger.debug(
            f"{DEBUG_TAG} | Eliminate negative charge: {atom.OBAtom.GetIdx()} "
            f"with charge {to_add}"
        )
    return omol, given_charge


def clean_resonances_0(omol: pybel.Molecule) -> pybel.Molecule:
    smarts = pybel.Smarts("[*-]-[*]=[*]~[*+]")
    res = list(smarts.findall(omol))
    while len(res):
        idxs = res.pop(0)
        if (
            pt.GetDefaultValence(omol.OBMol.GetAtom(idxs[-1]).GetAtomicNum())
            > omol.OBMol.GetAtom(idxs[-1]).GetTotalValence()
        ):
            moloplogger.debug(f"Cleaning resonance 0: {idxs}")
            omol.OBMol.GetBond(idxs[0], idxs[1]).SetBondOrder(
                omol.OBMol.GetBond(idxs[0], idxs[1]).GetBondOrder() + 1
            )
            omol.OBMol.GetBond(idxs[1], idxs[2]).SetBondOrder(
                omol.OBMol.GetBond(idxs[1], idxs[2]).GetBondOrder() - 1
            )
            omol.OBMol.GetAtom(idxs[0]).SetFormalCharge(
                omol.OBMol.GetAtom(idxs[0]).GetFormalCharge() + 1
            )
            omol.OBMol.GetAtom(idxs[2]).SetFormalCharge(
                omol.OBMol.GetAtom(idxs[2]).GetFormalCharge() - 1
            )
    return omol


def clean_resonances_1(omol: pybel.Molecule) -> pybel.Molecule:
    """
    `[*-]=[*+]=[*]>>[#6]#[*+]-[*-]`

    Parameters:
        omol (pybel.Molecule): The input molecule object

    Returns:
        pybel.Molecule: The cleaned molecule.
    """
    smarts = pybel.Smarts("[*-]=[*+]=[*+0]")
    res = list(smarts.findall(omol))
    while len(res):
        idxs = res.pop(0)
        moloplogger.debug(f"Cleaning resonance 1: {idxs}")
        omol.OBMol.GetBond(idxs[1], idxs[2]).SetBondOrder(
            omol.OBMol.GetBond(idxs[1], idxs[2]).GetBondOrder() - 1
        )
        omol.OBMol.GetBond(idxs[0], idxs[1]).SetBondOrder(
            omol.OBMol.GetBond(idxs[0], idxs[1]).GetBondOrder() + 1
        )
        omol.OBMol.GetAtom(idxs[0]).SetFormalCharge(
            omol.OBMol.GetAtom(idxs[0]).GetFormalCharge() + 1
        )
        omol.OBMol.GetAtom(idxs[-1]).SetFormalCharge(
            omol.OBMol.GetAtom(idxs[0]).GetFormalCharge() - 1
        )
    return omol


def clean_resonances_2(omol: pybel.Molecule) -> pybel.Molecule:
    """
    `[#8]=[#6](-[!-])-[*]=[*]-[#7-,#6-]>>[#8-]-[#6](-[!-])=[*]-[*]=[#7,#6]`

    Parameters:
        omol (pybel.Molecule): The input molecule object

    Returns:
        pybel.Molecule: The cleaned molecule.
    """
    smarts = pybel.Smarts("[#8]=[#6](-[!-])-[*]=[*]-[#7-,#6-]")
    while res := smarts.findall(omol):
        idxs = res.pop(0)
        moloplogger.debug(f"Cleaning resonance 2: {idxs}")
        omol.OBMol.GetBond(idxs[4], idxs[5]).SetBondOrder(2)
        omol.OBMol.GetBond(idxs[3], idxs[4]).SetBondOrder(1)
        omol.OBMol.GetBond(idxs[1], idxs[3]).SetBondOrder(2)
        omol.OBMol.GetBond(idxs[0], idxs[1]).SetBondOrder(1)
        omol.OBMol.GetAtom(idxs[0]).SetFormalCharge(-1)
        omol.OBMol.GetAtom(idxs[-1]).SetFormalCharge(0)
    return omol


def clean_resonances_3(omol: pybel.Molecule) -> pybel.Molecule:
    """
    `[#7+]=[*]-[*]=[*]-[#8-]>>[#7]-[*]=[*]-[*]=[#8]`

    Parameters:
        omol (pybel.Molecule): The input molecule object

    Returns:
        pybel.Molecule: The cleaned molecule.
    """
    smarts = pybel.Smarts("[#7+]=[*]-[*]=[*]-[#8-]")
    res = list(smarts.findall(omol))
    while len(res):
        idxs = res.pop(0)
        moloplogger.debug(f"Cleaning resonance 3: {idxs}")
        omol.OBMol.GetBond(idxs[3], idxs[4]).SetBondOrder(2)
        omol.OBMol.GetBond(idxs[2], idxs[3]).SetBondOrder(1)
        omol.OBMol.GetBond(idxs[1], idxs[2]).SetBondOrder(2)
        omol.OBMol.GetBond(idxs[0], idxs[1]).SetBondOrder(1)
        omol.OBMol.GetAtom(idxs[0]).SetFormalCharge(0)
        omol.OBMol.GetAtom(idxs[-1]).SetFormalCharge(0)
    return omol


def clean_resonances_4(omol: pybel.Molecule) -> pybel.Molecule:
    """
    `[#7+,#8+]=[*]-[#6-,#7-,#8-]>>[#7,#8]=[*]-[#6,#7,#8]`

    Parameters:
        omol (pybel.Molecule): The input molecule object

    Returns:
        pybel.Molecule: The cleaned molecule.
    """
    smarts = pybel.Smarts("[#7+,#8+]=[*]-[#6-,#7-,#8-]")
    res = list(smarts.findall(omol))
    while len(res):
        idxs = res.pop(0)
        moloplogger.debug(f"Cleaning resonance 4: {idxs}")
        omol.OBMol.GetBond(idxs[1], idxs[2]).SetBondOrder(
            omol.OBMol.GetBond(idxs[1], idxs[2]).GetBondOrder() + 1
        )
        omol.OBMol.GetBond(idxs[0], idxs[1]).SetBondOrder(
            omol.OBMol.GetBond(idxs[0], idxs[1]).GetBondOrder() - 1
        )
        omol.OBMol.GetAtom(idxs[0]).SetFormalCharge(0)
        omol.OBMol.GetAtom(idxs[-1]).SetFormalCharge(0)
    return omol


def clean_resonances_5(omol: pybel.Molecule) -> pybel.Molecule:
    """
    `[#8,#16]=[*]-[#6-,#7-]>>[#8-,#16-]-[*]=[#6,#7]`

    Parameters:
        omol (pybel.Molecule): The input molecule object

    Returns:
        pybel.Molecule: The cleaned molecule.
    """
    smarts = pybel.Smarts("[#7+0,#8+0,#16+0]=[*+0]-[#6-,#7-]")
    res = list(smarts.findall(omol))
    while len(res):
        idxs = res.pop(0)
        moloplogger.debug(f"Cleaning resonance 5: {idxs}")
        omol.OBMol.GetBond(idxs[1], idxs[2]).SetBondOrder(
            omol.OBMol.GetBond(idxs[1], idxs[2]).GetBondOrder() + 1
        )
        omol.OBMol.GetBond(idxs[0], idxs[1]).SetBondOrder(
            omol.OBMol.GetBond(idxs[0], idxs[1]).GetBondOrder() - 1
        )
        omol.OBMol.GetAtom(idxs[0]).SetFormalCharge(-1)
        omol.OBMol.GetAtom(idxs[-1]).SetFormalCharge(0)
    return omol


def clean_resonances_6(omol: pybel.Molecule) -> pybel.Molecule:
    """
    `[#6]=[#6]=[#6-,#7-]>>[#6-]-[#6]#[#6,#7]`

    Parameters:
        omol (pybel.Molecule): The input molecule object

    Returns:
        pybel.Molecule: The cleaned molecule.
    """
    smarts = pybel.Smarts("[#6]=[#6]=[#6-,#7-]")
    res = list(smarts.findall(omol))
    while len(res):
        idxs = res.pop(0)
        moloplogger.debug(f"Cleaning resonance 6: {idxs}")
        omol.OBMol.GetBond(idxs[1], idxs[2]).SetBondOrder(3)
        omol.OBMol.GetBond(idxs[0], idxs[1]).SetBondOrder(1)
        omol.OBMol.GetAtom(idxs[0]).SetFormalCharge(-1)
        omol.OBMol.GetAtom(idxs[-1]).SetFormalCharge(0)
    return omol


def clean_resonances_7(omol: pybel.Molecule) -> pybel.Molecule:
    """
    `[*-]1-[*](=[*])-[*]=[*]-[*]=[*]1>>[*]1=[*](-[*-])-[*]=[*]-[*]=[*]1`

    Parameters:
        omol (pybel.Molecule): The input molecule object

    Returns:
        pybel.Molecule: The cleaned molecule.
    """
    smarts = pybel.Smarts("[*-]1-[*](=[*])-[*]=[*]-[*]=[*]1")
    res = list(smarts.findall(omol))
    while len(res):
        idxs = res.pop(0)
        moloplogger.debug(f"Cleaning resonance 7: {idxs}")
        omol.OBMol.GetBond(idxs[1], idxs[2]).SetBondOrder(1)
        omol.OBMol.GetBond(idxs[0], idxs[1]).SetBondOrder(2)
        omol.OBMol.GetAtom(idxs[0]).SetFormalCharge(0)
        omol.OBMol.GetAtom(idxs[-1]).SetFormalCharge(-1)
    return omol


def clean_resonances_8(omol: pybel.Molecule) -> pybel.Molecule:
    """
    `[*-]1-[*]=[*]-[*](=[*])-[*]=[*]1>>[*]1=[*]-[*]=[*](-[*-])-[*]=[*]1`

    Parameters:
        omol (pybel.Molecule): The input molecule object

    Returns:
        pybel.Molecule: The cleaned molecule.
    """
    smarts = pybel.Smarts("[*-]1-[*]=[*]-[*](=[*])-[*]=[*]1")
    res = list(smarts.findall(omol))
    while len(res):
        idxs = res.pop(0)
        moloplogger.debug(f"Cleaning resonance 8: {idxs}")
        omol.OBMol.GetBond(idxs[3], idxs[4]).SetBondOrder(1)
        omol.OBMol.GetBond(idxs[2], idxs[3]).SetBondOrder(2)
        omol.OBMol.GetBond(idxs[1], idxs[2]).SetBondOrder(1)
        omol.OBMol.GetBond(idxs[0], idxs[1]).SetBondOrder(2)
        omol.OBMol.GetAtom(idxs[0]).SetFormalCharge(0)
        omol.OBMol.GetAtom(idxs[-1]).SetFormalCharge(-1)
    return omol


def clean_resonances_9(omol: pybel.Molecule) -> pybel.Molecule:
    """
    Clean splited charges

    `[*+]-[*-]>>[*]=[*]`

    Parameters:
        omol (pybel.Molecule): The input molecule object

    Returns:
        pybel.Molecule: The cleaned molecule.
    """
    smarts = pybel.Smarts("[*+]-,=[*-]")
    res = list(smarts.findall(omol))
    while len(res):
        idxs = res.pop(0)
        pos_atom = omol.OBMol.GetAtom(idxs[0])
        neg_atom = omol.OBMol.GetAtom(idxs[1])
        if (
            pt.GetDefaultValence(pos_atom.GetAtomicNum()) - pos_atom.GetTotalValence()
            >= 1
            and pt.GetDefaultValence(neg_atom.GetAtomicNum())
            - neg_atom.GetTotalValence()
            >= 1
        ):
            moloplogger.debug(f"Cleaning splited charges {idxs}")
            omol.OBMol.GetBond(idxs[0], idxs[1]).SetBondOrder(
                int(omol.OBMol.GetBond(idxs[0], idxs[1]).GetBondOrder() + 1)
            )
            pos_atom.SetFormalCharge(0)
            neg_atom.SetFormalCharge(0)
    return omol


def clean_resonances_10(omol: pybel.Molecule) -> pybel.Molecule:
    """
    Clean double radical

    `[*]-[*]=,#[*]-[*]>>[*]-[*]-[*]`

    Parameters:
        omol (pybel.Molecule): The input molecule object

    Returns:
        pybel.Molecule: The cleaned molecule.
    """
    smarts = pybel.Smarts("[*]-[*]=,#[*]-[*]")
    res = list(smarts.findall(omol))
    while len(res):
        idxs = res.pop(0)
        atom1 = omol.OBMol.GetAtom(idxs[0])
        atom2 = omol.OBMol.GetAtom(idxs[-1])
        if get_radical_number(atom1) and get_radical_number(atom2):
            moloplogger.debug(f"Cleaning double radical")
            omol.OBMol.GetBond(idxs[1], idxs[2]).SetBondOrder(
                int(omol.OBMol.GetBond(idxs[1], idxs[2]).GetBondOrder() - 1)
            )
            omol.OBMol.GetBond(idxs[0], idxs[1]).SetBondOrder(
                int(omol.OBMol.GetBond(idxs[0], idxs[1]).GetBondOrder() + 1)
            )
            omol.OBMol.GetBond(idxs[2], idxs[3]).SetBondOrder(
                int(omol.OBMol.GetBond(idxs[2], idxs[3]).GetBondOrder() + 1)
            )
    return omol


def clean_resonances_11(omol: pybel.Molecule) -> pybel.Molecule:
    smarts = pybel.Smarts("[#7v3+0,#8v2+0,#16v2+0]-,=[*+1]")
    res = list(smarts.findall(omol))
    while len(res):
        idxs = res.pop(0)
        atom1 = omol.OBMol.GetAtom(idxs[1])
        if pt.GetDefaultValence(atom1.GetAtomicNum()) - atom1.GetTotalValence() >= 1:
            moloplogger.debug(f"Cleaning N, O, S with neighboring positive charge")
            omol.OBMol.GetBond(idxs[0], idxs[1]).SetBondOrder(
                int(omol.OBMol.GetBond(idxs[0], idxs[1]).GetBondOrder() + 1)
            )
            omol.OBMol.GetAtom(idxs[0]).SetFormalCharge(1)
            omol.OBMol.GetAtom(idxs[1]).SetFormalCharge(0)
    return omol


def clean_resonances_12(omol: pybel.Molecule) -> pybel.Molecule:
    smarts = pybel.Smarts("[#7v3+0,#8v2+0,#16v2+0]-[*]=[*]-[*+1]")
    res = list(smarts.findall(omol))
    while len(res):
        idxs = res.pop(0)
        atom1 = omol.OBMol.GetAtom(idxs[0])
        atom4 = omol.OBMol.GetAtom(idxs[3])
        if pt.GetDefaultValence(atom4.GetAtomicNum()) - atom4.GetTotalValence() >= 1:
            moloplogger.debug(f"Cleaning N, O, S with neighboring positive charge")
            omol.OBMol.GetBond(idxs[0], idxs[1]).SetBondOrder(
                int(omol.OBMol.GetBond(idxs[0], idxs[1]).GetBondOrder() + 1)
            )
            omol.OBMol.GetBond(idxs[1], idxs[2]).SetBondOrder(
                int(omol.OBMol.GetBond(idxs[1], idxs[2]).GetBondOrder() - 1)
            )
            omol.OBMol.GetBond(idxs[2], idxs[3]).SetBondOrder(
                int(omol.OBMol.GetBond(idxs[2], idxs[3]).GetBondOrder() + 1)
            )
            omol.OBMol.GetAtom(idxs[0]).SetFormalCharge(1)
            omol.OBMol.GetAtom(idxs[3]).SetFormalCharge(0)
    return omol


def clean_resonances(omol: pybel.Molecule) -> pybel.Molecule:
    omol = clean_resonances_0(omol)
    omol = clean_resonances_1(omol)
    omol = clean_resonances_2(omol)
    omol = clean_resonances_3(omol)
    omol = clean_resonances_4(omol)
    omol = clean_resonances_9(omol)
    omol = clean_resonances_5(omol)
    omol = clean_resonances_6(omol)
    omol = clean_resonances_7(omol)
    omol = clean_resonances_8(omol)
    omol = clean_resonances_9(omol)
    omol = clean_resonances_10(omol)
    omol = clean_resonances_11(omol)
    omol = clean_resonances_12(omol)
    return omol


def xyz2omol(
    xyz_block: str, total_charge: int = 0, total_radical_electrons: int = 0
) -> Union[pybel.Molecule, None]:
    moloplogger.debug(
        f"{DEBUG_TAG} | Structure recovery started"
        f" with total charge {total_charge} and total radical {total_radical_electrons}"
    )
    if total_radical_electrons < 0:
        moloplogger.debug(
            f"{DEBUG_TAG} | Cannot recover the structure with negative radical"
        )
        return None
    omol = pybel.readstring("xyz", xyz_block)
    moloplogger.debug(f"{DEBUG_TAG} | Input smiles: {omol.write('smi')}")
    for atom in omol.atoms:
        atom.OBAtom.SetFormalCharge(0)
    omol = make_connections(omol)
    omol = pre_clean(omol)

    for atom in omol.atoms:
        if get_under_bonded_number(atom.OBAtom) < 0:
            atom.OBAtom.SetFormalCharge(-get_under_bonded_number(atom.OBAtom))
        if atom.OBAtom.GetAtomicNum() == 5 and get_under_bonded_number(atom.OBAtom) > 0:
            atom.OBAtom.SetFormalCharge(-get_under_bonded_number(atom.OBAtom))
    given_charge = total_charge - sum(
        atom.OBAtom.GetFormalCharge() for atom in omol.atoms
    )

    omol, given_charge = eliminate_NNN(omol, given_charge)
    omol, given_charge = eliminate_high_positive_charge_atoms(omol, given_charge)
    omol, given_charge = eliminate_CN_in_doubt(omol, given_charge)
    omol, given_charge = eliminate_carboxyl(omol, given_charge)
    omol = clean_carbine_neighbor_unsaturated(omol)
    omol, given_charge = eliminate_carbine_neighbor_heteroatom(omol, given_charge)
    omol = clean_neighbor_radicals(omol)
    omol = clean_carbine_neighbor_unsaturated(omol)
    omol, given_charge = eliminate_charge_spliting(omol, given_charge)
    omol = break_one_bond(omol, given_charge, total_radical_electrons)

    moloplogger.debug(
        f"{DEBUG_TAG} | Given charge: {given_charge}, total charge: {total_charge}, smiles: {omol.write('smi')}"
    )
    possible_resonances = get_radical_resonances(omol)
    moloplogger.debug(
        f"{DEBUG_TAG} | Possible resonance structures number: {len(possible_resonances)}"
    )
    recovered_resonances = []
    for resonance in possible_resonances:
        charge = given_charge
        moloplogger.debug(
            f"{DEBUG_TAG} | charge to be allocated: {charge}, smiles: {resonance.write('smi')}"
        )
        resonance, charge = eliminate_1_3_dipole(resonance, charge)
        resonance, charge = eliminate_positive_charges(resonance, charge)
        resonance, charge = eliminate_negative_charges(resonance, charge)
        resonance = clean_neighbor_radicals(resonance)
        moloplogger.debug(
            f"{DEBUG_TAG} | charge to be allocated: {charge}, pre-cleaned resonance smiles: {resonance.write('smi')}"
        )
        resonance = clean_resonances(resonance)
        moloplogger.debug(
            f"{DEBUG_TAG} | charge to be allocated: {charge}, cleaned resonance smiles: {resonance.write('smi')}"
        )
        if final_check_omol(resonance, total_charge, total_radical_electrons):
            recovered_resonances.append((resonance, charge, total_radical_electrons))
    if len(recovered_resonances) == 0:
        moloplogger.debug(
            f"{DEBUG_TAG} | Cannot recover the structure with given charge and radical"
        )
        return None

    recovered_resonances.sort(key=omol_score)
    final_omol = recovered_resonances[0][0]

    moloplogger.debug(
        f"{DEBUG_TAG} | Final charge: {sum(atom.OBAtom.GetFormalCharge() for atom in final_omol.atoms)} "
        f"total charge: {total_charge}; Final radical: {get_radical_number(final_omol.OBMol.GetAtom(1))}"
        f" total radical: {total_radical_electrons}; Final smiles: {final_omol.write('smi')}"
    )
    return final_omol


def rdmol_structure_recovery(
    xyz_block: str, total_charge: int = 0, total_radical_electrons: int = 0
) -> Union[Chem.rdchem.RWMol, None]:
    omol = xyz2omol(xyz_block, total_charge, total_radical_electrons)
    if omol is None:
        return None
    rdmol = omol_to_rdmol_by_graph(omol)
    if rdmol is None:
        return None
    rwmol = Chem.RWMol(rdmol)
    return rwmol


def xyz2rwmol(
    xyz_block: str, total_charge: int = 0, total_radical_electrons: int = 0
) -> Union[Chem.rdchem.RWMol, None]:
    rwmol = Chem.RWMol(Chem.MolFromXYZBlock(xyz_block))

    removable_metal_atoms = [
        atom.GetIdx()
        for atom in rwmol.GetAtoms()
        if is_metal(atom.GetAtomicNum())
        and pt.GetElementSymbol(atom.GetAtomicNum()) in metal_valence_avialable_prior
    ]

    # Recursively remove metal atoms and guess the possible valence state of the metal,
    # and the remaining parts will be reconstructed with the new total charge
    for atom_idx in removable_metal_atoms:
        # Skip if the atom has non-zero formal charge, which means it has been processed
        if rwmol.GetAtomWithIdx(atom_idx).GetFormalCharge() != 0:
            continue
        # If the metal atom is the last atom, set the formal charge and radical number and return the molecule
        if rwmol.GetNumAtoms() == 1:
            rwmol.GetAtomWithIdx(atom_idx).SetFormalCharge(total_charge)
            rwmol.GetAtomWithIdx(atom_idx).SetNumRadicalElectrons(
                total_radical_electrons
            )
            return rwmol
        # Try to remove the metal atom with different valence states
        temp_rwmol = Chem.RWMol(rwmol)
        metal_coords = temp_rwmol.GetConformer().GetAtomPosition(atom_idx)
        temp_rwmol.RemoveAtom(atom_idx)
        possible_rwmols = []
        # Try different valence states of the metal atom
        for possible_metal_valence in metal_valence_avialable_prior[
            rwmol.GetAtomWithIdx(atom_idx).GetSymbol()
        ]:
            for radical_num in get_possible_metal_radicals(
                rwmol.GetAtomWithIdx(atom_idx).GetSymbol(), possible_metal_valence
            ):
                # Reconstruct the remaining parts with the new total charge
                moloplogger.debug(
                    f"{DEBUG_TAG} | Try metal atom {rwmol.GetAtomWithIdx(atom_idx).GetSymbol()} "
                    f"idx {atom_idx} with valence {possible_metal_valence} and radical {radical_num}"
                )
                sub_rwmol = xyz2rwmol(
                    Chem.MolToXYZBlock(temp_rwmol),
                    total_charge - possible_metal_valence,
                    total_radical_electrons - radical_num,
                )
                if sub_rwmol is not None:
                    metal_idx = sub_rwmol.AddAtom(rwmol.GetAtomWithIdx(atom_idx))
                    sub_rwmol.GetAtomWithIdx(metal_idx).SetFormalCharge(
                        possible_metal_valence
                    )
                    sub_rwmol.GetAtomWithIdx(metal_idx).SetNumRadicalElectrons(0)
                    sub_rwmol.GetConformer().SetAtomPosition(metal_idx, metal_coords)
                    possible_rwmols.append(
                        reset_atom_index(
                            sub_rwmol,
                            list(range(0, atom_idx))
                            + [metal_idx]
                            + list(range(atom_idx, metal_idx)),
                        )
                    )
        # Choose the best structure among the possible structures, ranked by the score function
        if possible_rwmols:
            possible_rwmols = sorted(possible_rwmols, key=structure_score, reverse=True)
            smiles_list = "\n".join(
                Chem.CanonSmiles(Chem.MolToSmiles(rwmol)) for rwmol in possible_rwmols
            )
            moloplogger.debug(
                f"{DEBUG_TAG} | Possible resonance structures: \n{smiles_list}"
            )
            return Chem.RWMol(possible_rwmols[0])

        for possible_metal_valence in metal_valence_avialable_minor[
            rwmol.GetAtomWithIdx(atom_idx).GetSymbol()
        ]:
            for radical_num in get_possible_metal_radicals(
                rwmol.GetAtomWithIdx(atom_idx).GetSymbol(), possible_metal_valence
            ):
                # Reconstruct the remaining parts with the new total charge
                moloplogger.debug(
                    f"{DEBUG_TAG} | Try metal atom {rwmol.GetAtomWithIdx(atom_idx).GetSymbol()} "
                    f"idx {atom_idx} with unusual valence {possible_metal_valence} and radical {radical_num}"
                )
                sub_rwmol = xyz2rwmol(
                    Chem.MolToXYZBlock(temp_rwmol),
                    total_charge - possible_metal_valence,
                    total_radical_electrons - radical_num,
                )
                if sub_rwmol is not None:
                    metal_idx = sub_rwmol.AddAtom(rwmol.GetAtomWithIdx(atom_idx))
                    sub_rwmol.GetAtomWithIdx(metal_idx).SetFormalCharge(
                        possible_metal_valence
                    )
                    sub_rwmol.GetAtomWithIdx(metal_idx).SetNumRadicalElectrons(0)
                    sub_rwmol.GetConformer().SetAtomPosition(metal_idx, metal_coords)
                    possible_rwmols.append(
                        reset_atom_index(
                            sub_rwmol,
                            list(range(0, atom_idx))
                            + [metal_idx]
                            + list(range(atom_idx, metal_idx)),
                        )
                    )
        if possible_rwmols:
            possible_rwmols = sorted(possible_rwmols, key=structure_score, reverse=True)
            smiles_list = "\n".join(
                Chem.CanonSmiles(Chem.MolToSmiles(rwmol)) for rwmol in possible_rwmols
            )
            moloplogger.debug(
                f"{DEBUG_TAG} | Possible resonance structures: \n{smiles_list}"
            )
            return Chem.RWMol(possible_rwmols[0])
        else:
            return None

    # If no metal atom is found, try to recover the structure
    rwmol = rdmol_structure_recovery(xyz_block, total_charge, total_radical_electrons)
    if rwmol is None:
        return None
    if final_check(rwmol, total_charge, total_radical_electrons):
        return rwmol
    return None


def make_dative_bonds(rwmol: Chem.rdchem.RWMol, ratio=1.3) -> Chem.rdchem.RWMol:
    """
    Make dative bonds between the metal atoms and the non-metal atoms.
    Parameters:
        rwmol (Chem.rdchem.RWMol): The editable rdkit molecule.

    Returns:
        Chem.rdchem.RWMol: The editable rdkit molecule with dative bonds.
    """
    Chem.Kekulize(rwmol, clearAromaticFlags=True)
    for atom in rwmol.GetAtoms():
        atom.SetNoImplicit(True)
    exp_ratio = {
        "N": 1.45,
        "O": ratio,
        "S": ratio,
        "P": ratio,
    }
    metal_atoms = [
        atom.GetIdx() for atom in rwmol.GetAtoms() if is_metal(atom.GetAtomicNum())
    ]
    datived_rwmol = []
    for shuffled_metal_atoms in itertools.permutations(metal_atoms, len(metal_atoms)):
        temp_rwmol = Chem.RWMol(rwmol)
        for metal_atom in shuffled_metal_atoms:
            negative_atoms = [
                atom.GetIdx()
                for atom in temp_rwmol.GetAtoms()
                if atom.GetIdx() != metal_atom
                and atom.GetFormalCharge() < 0
                and temp_rwmol.GetConformer()
                .GetAtomPosition(atom.GetIdx())
                .Distance(temp_rwmol.GetConformer().GetAtomPosition(metal_atom))
                <= ratio
                * (
                    pt.GetRcovalent(
                        temp_rwmol.GetAtomWithIdx(metal_atom).GetAtomicNum()
                    )
                    + pt.GetRcovalent(
                        temp_rwmol.GetAtomWithIdx(atom.GetIdx()).GetAtomicNum()
                    )
                )
            ]
            negative_atoms.sort(
                key=lambda x: temp_rwmol.GetConformer()
                .GetAtomPosition(x)
                .Distance(temp_rwmol.GetConformer().GetAtomPosition(metal_atom))
            )
            while (
                temp_rwmol.GetAtomWithIdx(metal_atom).GetFormalCharge() > 0
                and negative_atoms
            ):
                negative_atom = negative_atoms.pop(0)
                temp_rwmol.AddBond(
                    metal_atom,
                    negative_atom,
                    bond_list[
                        -temp_rwmol.GetAtomWithIdx(negative_atom).GetFormalCharge()
                    ],
                )
                temp_rwmol.GetAtomWithIdx(metal_atom).SetFormalCharge(
                    temp_rwmol.GetAtomWithIdx(metal_atom).GetFormalCharge()
                    + temp_rwmol.GetAtomWithIdx(negative_atom).GetFormalCharge()
                )
                temp_rwmol.GetAtomWithIdx(negative_atom).SetFormalCharge(0)
        for metal_atom in metal_atoms:
            for dative_atom in [
                idxs[0]
                for idxs in temp_rwmol.GetSubstructMatches(
                    Chem.MolFromSmarts("[Ov2+0,Sv2+0,Nv3+0,Pv3+0]")
                )
            ]:
                if not temp_rwmol.GetBondBetweenAtoms(metal_atom, dative_atom):
                    if temp_rwmol.GetConformer().GetAtomPosition(metal_atom).Distance(
                        temp_rwmol.GetConformer().GetAtomPosition(dative_atom)
                    ) < exp_ratio[
                        temp_rwmol.GetAtomWithIdx(dative_atom).GetSymbol()
                    ] * (
                        pt.GetRcovalent(
                            temp_rwmol.GetAtomWithIdx(metal_atom).GetAtomicNum()
                        )
                        + pt.GetRcovalent(
                            temp_rwmol.GetAtomWithIdx(dative_atom).GetAtomicNum()
                        )
                    ):
                        temp_rwmol.AddBond(
                            dative_atom, metal_atom, Chem.BondType.DATIVE
                        )
        remained_negative_atoms = [
            atom.GetIdx()
            for atom in temp_rwmol.GetAtoms()
            if atom.GetFormalCharge() < 0 and not is_metal(atom.GetAtomicNum())
        ]
        for remained_negative_atom in remained_negative_atoms:
            metal_atoms = [
                atom.GetIdx()
                for atom in temp_rwmol.GetAtoms()
                if is_metal(atom.GetAtomicNum())
            ]
            if not metal_atoms:
                continue
            metal_atoms.sort(
                key=lambda x: temp_rwmol.GetConformer()
                .GetAtomPosition(x)
                .Distance(
                    temp_rwmol.GetConformer().GetAtomPosition(remained_negative_atom)
                )
            )
            if temp_rwmol.GetConformer().GetAtomPosition(metal_atoms[0]).Distance(
                temp_rwmol.GetConformer().GetAtomPosition(remained_negative_atom)
            ) <= ratio * (
                pt.GetRcovalent(
                    temp_rwmol.GetAtomWithIdx(metal_atoms[0]).GetAtomicNum()
                )
                + pt.GetRcovalent(
                    temp_rwmol.GetAtomWithIdx(remained_negative_atom).GetAtomicNum()
                )
            ):
                temp_rwmol.AddBond(
                    remained_negative_atom,
                    metal_atoms[0],
                    bond_list[
                        -temp_rwmol.GetAtomWithIdx(
                            remained_negative_atom
                        ).GetFormalCharge()
                    ],
                )
                temp_rwmol.GetAtomWithIdx(metal_atoms[0]).SetFormalCharge(
                    temp_rwmol.GetAtomWithIdx(metal_atoms[0]).GetFormalCharge()
                    + temp_rwmol.GetAtomWithIdx(
                        remained_negative_atom
                    ).GetFormalCharge()
                )
                temp_rwmol.GetAtomWithIdx(remained_negative_atom).SetFormalCharge(0)
        datived_rwmol.append(temp_rwmol)
    datived_rwmol_smiles = [Chem.MolToSmiles(rwmol) for rwmol in datived_rwmol]
    moloplogger.debug(
        f"{DEBUG_TAG} | Possible resonance structures with dative bonds: \n{datived_rwmol_smiles}"
    )
    return sorted(
        datived_rwmol,
        key=lambda x: sum(abs(atom.GetFormalCharge()) for atom in x.GetAtoms()),
    )[0]


def xyz2rdmol(
    xyz_block: str,
    total_charge: int = 0,
    total_radical_electrons: int = 0,
    make_dative: bool = True,
) -> Union[Chem.rdchem.Mol, None]:
    rwmol = xyz2rwmol(xyz_block, total_charge, total_radical_electrons)
    if rwmol is None:
        return None
    if make_dative:
        rwmol = make_dative_bonds(rwmol)
    moloplogger.debug(
        f"{DEBUG_TAG} | Final charge: {sum(atom.GetFormalCharge() for atom in rwmol.GetAtoms())} "
        f"total charge: {total_charge}; Final radical: {sum(atom.GetNumRadicalElectrons() for atom in rwmol.GetAtoms())}"
        f" total radical: {total_radical_electrons}; Final smiles: {Chem.MolToSmiles(rwmol)}"
    )
    return rwmol.GetMol()


def rdmol_check(
    rdmol: Chem.rdchem.Mol, total_atom_num: int, total_charge=0, total_radical=0
):
    """
    Check whether the rdmol is the same as the target molecule.
    Parameters:
        rdmol (Chem.rdchem.Mol): The rdmol to be checked.
        total_atom_num (int): The total atom number of the target molecule.
        total_charge (int): The total charge of the target molecule.
        total_radical (int): The total radical number of the target molecule.
    Returns:
        bool: Whether the rdmol is the same as the target molecule.
    """
    charge = 0
    radical = 0
    has_metal = False
    if rdmol.GetNumAtoms() != total_atom_num:
        return False
    for atom in rdmol.GetAtoms():
        charge += atom.GetFormalCharge()
        if atom.GetNumRadicalElectrons() != 2:
            radical += atom.GetNumRadicalElectrons()
        if is_metal(atom.GetAtomicNum()):
            has_metal = True
    moloplogger.debug(
        f"{DEBUG_TAG} check rdmol, target charge: {total_charge} charge: {charge}, "
        f"target radical: {total_radical} radical: {radical}"
    )
    if has_metal:
        return charge == total_charge
    else:
        return charge == total_charge and radical == total_radical


def omol_to_rdmol(
    omol: pybel.Molecule, total_charge=0, total_radical=0
) -> Union[Chem.rdchem.Mol, None]:
    """
    Convert a pybel molecule to a rdkit molecule.

    Parameters:
        omol (pybel.Molecule): The pybel molecule to be converted.
        total_charge (int): The total charge of the target molecule.
        total_radical (int): The total radical number of the target molecule.
    Returns:
        Chem.rdchem.Mol: The rdkit molecule.
    """
    total_atom_num = omol.OBMol.NumAtoms()
    moloplogger.debug("try tranform by mol2")
    rdmol = Chem.MolFromMol2Block(omol.write("mol2"), removeHs=False)
    if rdmol is not None:
        if rdmol_check(rdmol, total_atom_num, total_charge, total_radical):
            return rdmol
    moloplogger.debug("try tranform by graph")
    rdmol = omol_to_rdmol_by_graph(omol)
    if rdmol is not None:
        if rdmol_check(rdmol, total_atom_num, total_charge, total_radical):
            return rdmol
    return None
