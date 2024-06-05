"""
This module contains functions for recovering the structure of the molecule.

Offer a molecular graph recovery algorithm from the simple coodinates of atoms based on the initial work by [OpenBabel](https://openbabel.org/index.html), which can be easily used in the file reading process. This algorithm is different from the rdDetermineBonds (Original code implemented by [Jensen group](https://github.com/jensengroup/xyz2mol) and integrated in RDKit from the 2022.09 release, which is not suitable for the free radicals and complex containing metal.
  
Although our algorithm overcome the free radicals and metal problem and tested, it is still not perfect. There is no denying that, rdDetermineBonds works well for normal organic molecules. Thus, we would give molecule structure recovered by rdDetermineBonds first, if error happens, we will use our algorithm to recover the molecule structure instead. We hope that this strategy can take advantage of both approaches.
"""

# TODO SMARTS
import time
import itertools
from typing import List, Tuple

from openbabel import openbabel as ob
from openbabel import pybel
from rdkit import Chem

from molop.utils import is_metal
from molop.config import molopconfig
from molop.logger.logger import logger

pt = Chem.GetPeriodicTable()
HETEROATOM = (9, 8, 17, 7, 35, 54, 16, 34, 15)
DEBUG_TAG = "[DEBUG]|[STRUCTURE RECOVERY]"


def get_under_bonded_number(atom: ob.OBAtom) -> int:
    """
    OpenBabel determines that there are a large number of UNDER BONDED atoms in the molecule, which may be actual free radicals, positively or negatively charged.

    This function will finds out the UNDER BONDED number of the atom.

    Parameters:
        atom (ob.OBAtom): An object containing atomic information.

    Returns:
        The UNDER BONDED number of the atom.
    """
    if atom.IsMetal():
        return 0
    if atom.GetAtomicNum() in (7, 15):
        if atom.GetExplicitValence() == 5:
            # For nitrogen and phosphorus atoms, the radical number is 0 when the explicit valence is 5
            return 0
    elif atom.GetAtomicNum() == 16:
        if atom.GetExplicitValence() in (4, 6):
            # For sulfur atoms, the radical number is 0 when the explicit valence is 4 or 6
            return 0
        if atom.GetExplicitValence() == 5:
            # For sulfur atoms, the radical number is 0 when the explicit valence is 4 or 6
            return 1
    # For other atoms, the radical number is the default valence electrons plus formal charge minus the explicit valence
    if (
        pt.GetDefaultValence(atom.GetAtomicNum()) > atom.GetExplicitValence()
        and atom.GetFormalCharge() > 0
    ):
        return (
            pt.GetDefaultValence(atom.GetAtomicNum())
            - atom.GetFormalCharge()
            - atom.GetExplicitValence()
        )
    return (
        pt.GetDefaultValence(atom.GetAtomicNum())
        + atom.GetFormalCharge()
        - atom.GetExplicitValence()
    )


def clean_neighbor_radicals(
    omol: pybel.Molecule, given_charge: int = 0, given_radical: int = 0
):
    """
    This function cleans up the neighbor radicals of a given molecule.

    Parameters:
        omol (pybel.Molecule): The molecule to clean.
    """
    radical_atoms: List[int] = [
        ratom.idx for ratom in omol.atoms if get_under_bonded_number(ratom.OBAtom) >= 1
    ]
    if has_metal(omol):
        deradical_num = sum(
            get_under_bonded_number(ratom.OBAtom) for ratom in omol.atoms
        ) - abs(given_charge)
    else:
        deradical_num = (
            sum(get_under_bonded_number(ratom.OBAtom) for ratom in omol.atoms)
            - abs(given_charge)
            - abs(given_radical)
        )
    if len(radical_atoms) < 2:
        return
    logger.debug(
        f"{DEBUG_TAG} cleaning radicals {radical_atoms}, given_charge={given_charge}, given_radical={given_radical}"
    )
    for radical_atom_1, radical_atom_2 in list(
        itertools.permutations(radical_atoms, 2)
    ):
        if deradical_num <= 1:
            return
        now_bonding = (
            False
            if omol.OBMol.GetBond(radical_atom_1, radical_atom_2) is None
            else True
        )
        distance = omol.OBMol.GetAtom(radical_atom_1).GetDistance(
            omol.OBMol.GetAtom(radical_atom_2)
        )
        default_bond_length = pt.GetRcovalent(
            omol.OBMol.GetAtom(radical_atom_1).GetAtomicNum()
        ) + pt.GetRcovalent(omol.OBMol.GetAtom(radical_atom_2).GetAtomicNum())
        logger.debug(
            f"{DEBUG_TAG} cleaning radicals {pt.GetElementSymbol(omol.OBMol.GetAtom(radical_atom_1).GetAtomicNum())} "
            + f"{radical_atom_1} and {pt.GetElementSymbol(omol.OBMol.GetAtom(radical_atom_1).GetAtomicNum())} "
            + f"{radical_atom_2}, now bonding: {now_bonding}, distance: {distance}, default bond length: {default_bond_length}"
        )
        if (
            now_bonding
            and get_under_bonded_number(omol.OBMol.GetAtom(radical_atom_1)) >= 1
            and get_under_bonded_number(omol.OBMol.GetAtom(radical_atom_2)) >= 1
        ):  # the two atoms should be bonded
            """
            This code block checks if two radical atoms are bonded and if they have a valence greater than their explicit valence.
            If they do, it increases the bond order between them by the minimum radical value of the two atoms.
            """
            omol.OBMol.GetBond(radical_atom_1, radical_atom_2).SetBondOrder(
                omol.OBMol.GetBond(radical_atom_1, radical_atom_2).GetBondOrder()
                + min(
                    get_under_bonded_number(omol.OBMol.GetAtom(radical_atom_1)),
                    get_under_bonded_number(omol.OBMol.GetAtom(radical_atom_2)),
                ),
            )
            deradical_num -= 2
            logger.debug(
                f"{DEBUG_TAG} bonding {radical_atom_1} and {radical_atom_2}, smiles now is {omol.write('smi')}"
            )
            continue
        if (
            distance <= 1.3 * default_bond_length
            and not now_bonding
            and get_under_bonded_number(omol.OBMol.GetAtom(radical_atom_1))
            and get_under_bonded_number(omol.OBMol.GetAtom(radical_atom_2))
        ):
            """
            This code block checks if two radical atoms are within a certain distance and if they are not already bonded.
            If they are, it adds a bond between them with order 1.
            """
            omol.OBMol.AddBond(radical_atom_1, radical_atom_2, 1)
            logger.debug(
                f"{DEBUG_TAG} Add bonding {radical_atom_1} and {radical_atom_2}, smiles now is {omol.write('smi')}"
            )
            logger.debug(f"{DEBUG_TAG} smiles now is {omol.OBMol.NumAtoms()}")
            deradical_num -= 2
            continue
        if (
            distance <= 1.45 * default_bond_length
            and not now_bonding
            and get_under_bonded_number(omol.OBMol.GetAtom(radical_atom_1))
            and get_under_bonded_number(omol.OBMol.GetAtom(radical_atom_2))
            and omol.OBMol.GetAtom(radical_atom_1).GetAtomicNum() in HETEROATOM
            and not omol.OBMol.GetAtom(radical_atom_2).IsMetal()
        ):
            omol.OBMol.AddBond(radical_atom_1, radical_atom_2, 1)
            logger.debug(
                f"{DEBUG_TAG} Add bonding {radical_atom_1} and {radical_atom_2}, smiles now is {omol.write('smi')}"
            )
            deradical_num -= 2
            continue


def break_one_bond(
    omol: pybel.Molecule, charge_to_be_allowed: int = 0, given_radical: int = 0
):
    if (
        sum(get_under_bonded_number(atom.OBAtom) for atom in omol.atoms) == 0
        and abs(charge_to_be_allowed) + given_radical > 0
    ):
        logger.debug(f"{DEBUG_TAG} breaking bond")
        for atom_1 in omol.atoms:
            for atom_2 in ob.OBAtomAtomIter(atom_1.OBAtom):
                if atom_1.OBAtom.GetBond(atom_2).GetBondOrder() > 1:
                    logger.debug(
                        f"{DEBUG_TAG} break bond {atom_1.idx} and {atom_2.GetIdx()}"
                    )
                    atom_1.OBAtom.GetBond(atom_2).SetBondOrder(
                        atom_1.OBAtom.GetBond(atom_2).GetBondOrder() - 1
                    )
                    return


def fix_under_bonded_dipole(omol: pybel.Molecule):
    """
    Sometimes, the neighbored atoms with posstive and negative charges are under bonded.

    e.g. `[H][C+]=[C-][H]>>[H][C]#[C][H]`

    This function fixes under-bonded dipole moments in a molecule.

    Parameters:
        omol (pybel.Molecule): The input molecule.
    """
    # smarts = pybel.Smarts("[*+]-[*-]")
    for atom_1 in omol.atoms:
        if atom_1.OBAtom.GetFormalCharge() == 1:
            for atom_2 in ob.OBAtomAtomIter(atom_1.OBAtom):
                if atom_2.GetFormalCharge() == -1:
                    if (
                        pt.GetDefaultValence(atom_1.OBAtom.GetAtomicNum())
                        > atom_1.OBAtom.GetExplicitValence()
                        and pt.GetDefaultValence(atom_2.GetAtomicNum())
                        > atom_2.GetExplicitValence()
                    ):
                        logger.debug(f"{DEBUG_TAG} fix under bonded dipole")
                        atom_1.OBAtom.GetBond(atom_2).SetBondOrder(
                            atom_1.OBAtom.GetBond(atom_2).GetBondOrder() + 1
                        )
                        atom_1.OBAtom.SetFormalCharge(0)
                        atom_2.SetFormalCharge(0)


def fix_dipole_type_a(mol: pybel.Molecule):
    """
    (Type A): Dipole like `[CH2]=[O+]-[NH-]`.

    OpenBabel will set the central atom to be neutral and the other two atoms neutral with radical 1 like `[CH2]-[O]-[NH]`.
    Thus, find the combination of the two atoms with radical 1 and one heteroatom between them with radical 0.
    Negative charges are mutually resonant at any position of the atoms on either side.
    So it is straightforward to set one of the atoms as neutral and the other as charge -1.

    Known issues: Molecule like `C=[O+][N-]N[CH-]C(=O)C` can not distinguish between the part CON and NNC.

    Parameters:
        mol (pybel.Molecule): The molecule to be fixed.
    """
    # Get all atoms with radical 1
    radical_atoms = [
        atom for atom in mol.atoms if get_under_bonded_number(atom.OBAtom) == 1
    ]

    # Iterate through all combinations of radical atoms
    for atom_1, atom_2 in itertools.combinations(radical_atoms, 2):
        # Find the intersections between the bonds of atom_1 and atom_2
        intersections = set(
            (atom.idx for atom in mol.atoms if atom_1.OBAtom.GetBond(atom.OBAtom))
        ) & set((atom.idx for atom in mol.atoms if atom_2.OBAtom.GetBond(atom.OBAtom)))

        # If there is only one intersection, proceed with fixing the dipole
        if len(intersections) == 1:
            # Convert the intersection index to 0-based
            center_idx = intersections.pop() - 1

            # Check if the center atom satisfies the conditions for dipole fixing
            if (
                mol.atoms[center_idx].atomicnum not in (7, 8)
                or get_under_bonded_number(atom_1.OBAtom) != 1
                or get_under_bonded_number(atom_2.OBAtom) != 1
                or get_under_bonded_number(mol.atoms[center_idx].OBAtom) != 0
            ):
                continue

            # Check the bond orders between the center atom and the radical atoms
            if (
                mol.OBMol.GetBond(atom_1.idx, center_idx + 1).GetBondOrder() == 1
                and mol.OBMol.GetBond(atom_2.idx, center_idx + 1).GetBondOrder() == 1
            ):
                mol.OBMol.GetBond(atom_1.idx, center_idx + 1).SetBondOrder(2)
                mol.atoms[center_idx].OBAtom.SetFormalCharge(1)
                atom_2.OBAtom.SetFormalCharge(-1)
                logger.debug(
                    f"{DEBUG_TAG} Fix dipole a {mol.atoms[center_idx].OBAtom.GetIdx()}"
                )
            elif (
                mol.OBMol.GetBond(atom_1.idx, center_idx + 1).GetBondOrder() == 2
                and mol.OBMol.GetBond(atom_2.idx, center_idx + 1).GetBondOrder() == 1
            ):
                mol.OBMol.GetBond(atom_1.idx, center_idx + 1).SetBondOrder(3)
                mol.atoms[center_idx].OBAtom.SetFormalCharge(1)
                atom_2.OBAtom.SetFormalCharge(-1)
                logger.debug(
                    f"{DEBUG_TAG} Fix dipole a {mol.atoms[center_idx].OBAtom.GetIdx()}"
                )
            elif (
                mol.OBMol.GetBond(atom_1.idx, center_idx + 1).GetBondOrder() == 1
                and mol.OBMol.GetBond(atom_2.idx, center_idx + 1).GetBondOrder() == 2
            ):
                mol.OBMol.GetBond(atom_2.idx, center_idx + 1).SetBondOrder(3)
                mol.atoms[center_idx].OBAtom.SetFormalCharge(1)
                atom_1.OBAtom.SetFormalCharge(-1)
                logger.debug(
                    f"{DEBUG_TAG} Fix dipole a {mol.atoms[center_idx].OBAtom.GetIdx()}"
                )


def fix_dipole_type_b(mol: pybel.Molecule):
    """
    (Type B): Dipole like `[CH]#[N+]-[O-]`. This type only allow N (maybe P) to be the center.

    OpenBabel will set the central atom N (maybe P) with radical 1, the neutral atom with radical 2, and the negative atom with radical 1 like `[CH]-[N]-[O]`
    Thus, find the combination of the three atoms follow rule above.
    Negative charges are mutually resonant at any position of the atoms on either side.
    Nevertheless, I pact that atoms with radical 1 carry a formal charge -1.

    Parameters:
        mol (pybel.Molecule): The molecule to be fixed.
    """
    # Get radical 1 and radical 2 atoms
    radical_1_atoms = [
        atom for atom in mol.atoms if get_under_bonded_number(atom.OBAtom) == 1
    ]
    radical_2_atoms = [
        atom for atom in mol.atoms if get_under_bonded_number(atom.OBAtom) == 2
    ]

    # Iterate over all possible combinations of radical 1 atoms and radical 2 atoms
    for (atom_1, atom_2), atom_3 in itertools.product(
        itertools.permutations(radical_1_atoms, 2), radical_2_atoms
    ):
        # Check if atom 2 is not a nitrogen or phosphorus atom
        if atom_2.atomicnum not in (7, 15):
            continue

        # Check if bonds exist between atom 1, atom 2, and atom 3
        if mol.OBMol.GetBond(atom_1.idx, atom_2.idx) and mol.OBMol.GetBond(
            atom_2.idx, atom_3.idx
        ):
            # Set formal charges and bond order
            atom_1.OBAtom.SetFormalCharge(-1)
            atom_2.OBAtom.SetFormalCharge(1)
            mol.OBMol.GetBond(atom_2.idx, atom_3.idx).SetBondOrder(3)
            logger.debug(f"{DEBUG_TAG} Fix dipole b {atom_2.idx}")


def get_one_step_resonance(omol: pybel.Molecule) -> List[pybel.Molecule]:
    # Get a list of atoms with unpaired electrons (radicals)
    radical_atoms = [
        satom for satom in omol.atoms if get_under_bonded_number(satom.OBAtom) >= 1
    ]
    result = []
    # Iterate through the list of radical atoms
    for radical_atom in radical_atoms:
        # Iterate over neighboring atom 1 of the radical atom
        for neighbour_1_atom in ob.OBAtomAtomIter(radical_atom.OBAtom):
            # Iterate over neighboring atom 2 of neighboring atom 1
            for neighbour_2_atom in ob.OBAtomAtomIter(neighbour_1_atom):
                # If the bond order between neighbor 2 and neighbor 1 is greater than or equal to 2
                if (
                    neighbour_2_atom.GetBond(neighbour_1_atom).GetBondOrder() >= 2
                    and neighbour_1_atom.GetFormalCharge() == 0
                    and neighbour_2_atom.GetFormalCharge() == 0
                ):
                    # Create a new molecule object based on the input molecule
                    new_omol = pybel.Molecule(omol)
                    new_omol.OBMol.GetBond(
                        neighbour_2_atom.GetIdx(), neighbour_1_atom.GetIdx()
                    ).SetBondOrder(
                        new_omol.OBMol.GetBond(
                            neighbour_2_atom.GetIdx(), neighbour_1_atom.GetIdx()
                        ).GetBondOrder()
                        - 1
                    )
                    new_omol.OBMol.GetBond(
                        radical_atom.OBAtom.GetIdx(), neighbour_1_atom.GetIdx()
                    ).SetBondOrder(
                        new_omol.OBMol.GetBond(
                            radical_atom.OBAtom.GetIdx(), neighbour_1_atom.GetIdx()
                        ).GetBondOrder()
                        + 1
                    )
                    # clean_neighbor_radicals(new_omol)
                    result.append(new_omol)
    return result


def get_radical_resonances(omol: pybel.Molecule) -> List[pybel.Molecule]:
    """
    Retrieves a list of molecular structures that exhibit radical resonance.

    Parameters:
        omol (pybel.Molecule): The input molecule object

    Returns:
        A list of molecules with radical resonance structures
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


def omol_score(omol_tuple: Tuple[pybel.Molecule, int, int]) -> int:
    """
    Calculate the structural recovery score of a molecule, the lower the score the better the structural recovery.
    The criteria are 2 points for each unbonded electron and 1 point for each absolute value of charge.

    This function can only be used for comparison between isomers recovered from the same set of atomic coordinates.

    Args:
        omol_tuple (Tuple[pybel.Molecule, int, int]): A tuple containing a pybel.Molecule object and its charge and radical number.

    Returns:
        The structural recovery score of a molecule.
    """
    score = 0
    score += 2 * abs(
        sum(abs(get_under_bonded_number(atom.OBAtom)) for atom in omol_tuple[0].atoms)
        - omol_tuple[2]
    )
    score += sum(abs(atom.OBAtom.GetFormalCharge()) for atom in omol_tuple[0].atoms)
    return score


def clean_resonances_1(omol: pybel.Molecule) -> pybel.Molecule:
    """
    `[#8]=[#6](-[!-])-[*]=[*]-[#7-,#6-]>>[#8-]-[#6](-[!-])=[*]-[*]=[#7,#6]`

    Parameters:
        omol (pybel.Molecule): The input molecule object

    Returns:
        The cleaned molecule.
    """
    smarts = pybel.Smarts("[#8]=[#6](-[!-])-[*]=[*]-[#7-,#6-]")
    res = smarts.findall(omol)
    while len(res):
        logger.debug(f"Cleaning resonance 1: {res}")
        idxs = res[0]
        omol.OBMol.GetBond(idxs[4], idxs[5]).SetBondOrder(2)
        omol.OBMol.GetBond(idxs[3], idxs[4]).SetBondOrder(1)
        omol.OBMol.GetBond(idxs[1], idxs[3]).SetBondOrder(2)
        omol.OBMol.GetBond(idxs[0], idxs[1]).SetBondOrder(1)
        omol.atoms[idxs[0] - 1].OBAtom.SetFormalCharge(-1)
        omol.atoms[idxs[-1] - 1].OBAtom.SetFormalCharge(0)
        res = smarts.findall(omol)
    return omol


def clean_resonances_2(omol: pybel.Molecule) -> pybel.Molecule:
    """
    `[#7+]=[*]-[*]=[*]-[#8-]>>[#7]-[*]=[*]-[*]=[#8]`

    Parameters:
        omol (pybel.Molecule): The input molecule object

    Returns:
        The cleaned molecule.
    """
    smarts = pybel.Smarts("[#7+]=[*]-[*]=[*]-[#8-]")
    res = smarts.findall(omol)
    while len(res):
        logger.debug(f"Cleaning resonance 2: {res}")
        idxs = res[0]
        omol.OBMol.GetBond(idxs[3], idxs[4]).SetBondOrder(2)
        omol.OBMol.GetBond(idxs[2], idxs[3]).SetBondOrder(1)
        omol.OBMol.GetBond(idxs[1], idxs[2]).SetBondOrder(2)
        omol.OBMol.GetBond(idxs[0], idxs[1]).SetBondOrder(1)
        omol.atoms[idxs[0] - 1].OBAtom.SetFormalCharge(0)
        omol.atoms[idxs[-1] - 1].OBAtom.SetFormalCharge(0)
        res = smarts.findall(omol)
    return omol


def clean_resonances_3(omol: pybel.Molecule) -> pybel.Molecule:
    """
    `[#7+]=[*]-[#6-,#8-]>>[#7]-[*]=[#6,#8]`

    Parameters:
        omol (pybel.Molecule): The input molecule object

    Returns:
        The cleaned molecule.
    """
    smarts = pybel.Smarts("[#7+,#8+]=[*]-[#6-,#7-,#8-]")
    res = smarts.findall(omol)
    while len(res):
        logger.debug(f"Cleaning resonance 3: {res}")
        idxs = res[0]
        omol.OBMol.GetBond(idxs[1], idxs[2]).SetBondOrder(2)
        omol.OBMol.GetBond(idxs[0], idxs[1]).SetBondOrder(1)
        omol.atoms[idxs[0] - 1].OBAtom.SetFormalCharge(0)
        omol.atoms[idxs[-1] - 1].OBAtom.SetFormalCharge(0)
        res = smarts.findall(omol)
    return omol


def clean_resonances_4(omol: pybel.Molecule) -> pybel.Molecule:
    """
    `[#8,#16]=[#6]-[#6-,#7-]>>[#8-,#16-]-[#6]=[#6,#7]`

    Parameters:
        omol (pybel.Molecule): The input molecule object

    Returns:
        The cleaned molecule.
    """
    smarts = pybel.Smarts("[#8+0,#16+0]=[#6+0]-[#6-,#7-]")
    res = smarts.findall(omol)
    while len(res):
        logger.debug(f"Cleaning resonance 4: {res}")
        idxs = res[0]
        omol.OBMol.GetBond(idxs[1], idxs[2]).SetBondOrder(2)
        omol.OBMol.GetBond(idxs[0], idxs[1]).SetBondOrder(1)
        omol.atoms[idxs[0] - 1].OBAtom.SetFormalCharge(-1)
        omol.atoms[idxs[-1] - 1].OBAtom.SetFormalCharge(0)
        res = smarts.findall(omol)
    return omol


def clean_resonances_5(omol: pybel.Molecule) -> pybel.Molecule:
    """
    `[#6]=[#6]=[#6-,#7-]>>[#6-]-[#6]#[#6,#7]`

    Parameters:
        omol (pybel.Molecule): The input molecule object

    Returns:
        The cleaned molecule.
    """
    smarts = pybel.Smarts("[#6]=[#6]=[#6-,#7-]")
    res = smarts.findall(omol)
    while len(res):
        logger.debug(f"Cleaning resonance 5: {res}")
        idxs = res[0]
        # print(idxs)
        omol.OBMol.GetBond(idxs[1], idxs[2]).SetBondOrder(3)
        omol.OBMol.GetBond(idxs[0], idxs[1]).SetBondOrder(1)
        omol.atoms[idxs[0] - 1].OBAtom.SetFormalCharge(-1)
        omol.atoms[idxs[-1] - 1].OBAtom.SetFormalCharge(0)
        res = smarts.findall(omol)
    return omol


def clean_resonances_6(omol: pybel.Molecule) -> pybel.Molecule:
    """
    `[*-]=[*+]=[*]>>[#6]#[*+]-[*-]`

    Parameters:
        omol (pybel.Molecule): The input molecule object

    Returns:
        The cleaned molecule.
    """
    smarts = pybel.Smarts("[*-]=[*+]=[*]")
    res = smarts.findall(omol)
    while len(res):
        logger.debug(f"Cleaning resonance 6: {res}")
        idxs = res[0]
        omol.OBMol.GetBond(idxs[1], idxs[2]).SetBondOrder(1)
        omol.OBMol.GetBond(idxs[0], idxs[1]).SetBondOrder(3)
        omol.atoms[idxs[0] - 1].OBAtom.SetFormalCharge(0)
        omol.atoms[idxs[-1] - 1].OBAtom.SetFormalCharge(-1)
        res = smarts.findall(omol)
    return omol


def clean_resonances_7(omol: pybel.Molecule) -> pybel.Molecule:
    """
    `[*-]1-[*](=[*])-[*]=[*]-[*]=[*]1>>[*]1=[*](-[*-])-[*]=[*]-[*]=[*]1`

    Parameters:
        omol (pybel.Molecule): The input molecule object

    Returns:
        The cleaned molecule.
    """
    smarts = pybel.Smarts("[*-]1-[*](=[*])-[*]=[*]-[*]=[*]1")
    res = smarts.findall(omol)
    while len(res):
        logger.debug(f"Cleaning resonance 7: {res}")
        idxs = res[0]
        omol.OBMol.GetBond(idxs[1], idxs[2]).SetBondOrder(1)
        omol.OBMol.GetBond(idxs[0], idxs[1]).SetBondOrder(2)
        omol.atoms[idxs[0] - 1].OBAtom.SetFormalCharge(0)
        omol.atoms[idxs[-1] - 1].OBAtom.SetFormalCharge(-1)
        res = smarts.findall(omol)
    return omol


def clean_resonances_8(omol: pybel.Molecule) -> pybel.Molecule:
    """
    `[*-]1-[*]=[*]-[*](=[*])-[*]=[*]1>>[*]1=[*]-[*]=[*](-[*-])-[*]=[*]1`

    Parameters:
        omol (pybel.Molecule): The input molecule object

    Returns:
        The cleaned molecule.
    """
    smarts = pybel.Smarts("[*-]1-[*]=[*]-[*](=[*])-[*]=[*]1")
    res = smarts.findall(omol)
    while len(res):
        logger.debug(f"Cleaning resonance 8: {res}")
        idxs = res[0]
        omol.OBMol.GetBond(idxs[3], idxs[4]).SetBondOrder(1)
        omol.OBMol.GetBond(idxs[2], idxs[3]).SetBondOrder(2)
        omol.OBMol.GetBond(idxs[1], idxs[2]).SetBondOrder(1)
        omol.OBMol.GetBond(idxs[0], idxs[1]).SetBondOrder(2)
        omol.atoms[idxs[0] - 1].OBAtom.SetFormalCharge(0)
        omol.atoms[idxs[-1] - 1].OBAtom.SetFormalCharge(-1)
        res = smarts.findall(omol)
    return omol


def clean_resonances(omol: pybel.Molecule) -> pybel.Molecule:
    """
    Cleaning up resonance structures in molecules.

    Parameters:
        omol (pybel.Molecule): The input molecule object

    Returns:
        The cleaned molecule.
    """
    processes = [
        clean_resonances_6,
        clean_resonances_1,
        clean_resonances_2,
        clean_resonances_3,
        clean_resonances_4,
        clean_resonances_5,
        clean_resonances_7,
        clean_resonances_8,
    ]
    for process in processes:
        omol = process(omol)
    return omol


def fix_N_BCP(omol: pybel.Molecule):
    """
    Fix N-BCP structure.

    N-BCP:
    ```
          N
       /  |  \\
      /   |   \\
     C    C    C
      \   |   /
       \  |  /
          C
    ```
    N-BCP in Openbabel, will make one more bond between top N and bottom C:
    ```
          N
       / /| \\
      / / |  \\
     C |  C   C
      \ \ |  /
       \ \| /
          C
    ```

    Parameters:
        omol (pybel.Molecule): The input molecule object
    """
    smarts = pybel.Smarts("[#6]1([#6]2)([#6]3)[#7]23[#6]1")
    n_bcp = smarts.findall(omol)
    for n_bcp_idxs in n_bcp:
        for idx in n_bcp_idxs:
            indexs = set(n_bcp_idxs) - set([idx])
            if all(
                omol.atoms[idx - 1].OBAtom.GetBond(omol.atoms[idx_2 - 1].OBAtom)
                for idx_2 in indexs
            ):
                if omol.atoms[idx - 1].atomicnum == 7:
                    bcp_n = idx
                if omol.atoms[idx - 1].atomicnum == 6:
                    bcp_c = idx
        omol.OBMol.DeleteBond(omol.OBMol.GetBond(bcp_n, bcp_c))
        logger.debug(f"{DEBUG_TAG} Fix N-BCP: {bcp_n} - {bcp_c}")


def fix_cyclobutylamine(omol: pybel.Molecule):
    """
    Fix cyclobutylamine structure.
    """
    smarts = pybel.Smarts("[#6]1([#6]2)[#7]2[#6]1")
    cyclobutylamine = smarts.findall(omol)
    for cyclobutylamine_idxs in cyclobutylamine:
        for idx in cyclobutylamine_idxs:
            indexs = set(cyclobutylamine_idxs) - set([idx])
            if all(
                omol.atoms[idx - 1].OBAtom.GetBond(omol.atoms[idx_2 - 1].OBAtom)
                for idx_2 in indexs
            ):
                if omol.atoms[idx - 1].atomicnum == 7:
                    amine_n = idx
                if omol.atoms[idx - 1].atomicnum == 6:
                    butyl_c = idx
        omol.OBMol.DeleteBond(omol.OBMol.GetBond(amine_n, butyl_c))
        logger.debug(f"{DEBUG_TAG} Fix N-BCP: {amine_n} - {butyl_c}")


def fix_over_bonded_C(given_charge: int, omol: pybel.Molecule):
    """
    Sometimes openbabel will make C atom with 5 bonds.
    This function fixes the over bonded C atom in the molecule.

    Conditions:
        - `[Cv5]=[N]=[*]>>[Cv4]-[N]=[*]`
        - `[Cv5]=[X]>>[Cv4]-[X-]`, make a new negative charge
        - `[Cv5]=[C]>>[Cv4]-[C]`, make a new radical

    *This function can not be used separately.*

    Parameters:
        given_charge (int): The given charge to be fixed.
        omol (OBMol): The molecule object.
    """
    for atom in omol.atoms:
        if (
            atom.atomicnum == 6
            and atom.OBAtom.GetExplicitValence() == 5
            and atom.OBAtom.GetFormalCharge() == 0
        ):
            for neighbour_atom in ob.OBAtomAtomIter(atom.OBAtom):
                if (
                    neighbour_atom.GetAtomicNum() in HETEROATOM
                    and atom.OBAtom.GetBond(neighbour_atom).GetBondOrder() > 1
                    and atom.OBAtom.GetExplicitValence() == 5
                    and neighbour_atom.GetExplicitValence()
                    - pt.GetDefaultValence(neighbour_atom.GetAtomicNum())
                    == 1
                ):
                    atom.OBAtom.GetBond(neighbour_atom).SetBondOrder(
                        atom.OBAtom.GetBond(neighbour_atom).GetBondOrder() - 1
                    )
                    logger.debug(
                        f"{DEBUG_TAG} Fix over bonded C: {atom.OBAtom.GetIdx()} - {neighbour_atom.GetIdx()}, charge to be allocated: {given_charge}"
                    )
                    break
            for neighbour_atom in ob.OBAtomAtomIter(atom.OBAtom):
                if (
                    neighbour_atom.GetAtomicNum() in HETEROATOM
                    and atom.OBAtom.GetBond(neighbour_atom).GetBondOrder() > 1
                    and atom.OBAtom.GetExplicitValence() == 5
                ):
                    atom.OBAtom.GetBond(neighbour_atom).SetBondOrder(
                        atom.OBAtom.GetBond(neighbour_atom).GetBondOrder() - 1
                    )
                    if given_charge < 0:
                        neighbour_atom.SetFormalCharge(-1)
                        given_charge += 1
                    logger.debug(
                        f"{DEBUG_TAG} Fix over bonded C: {atom.OBAtom.GetIdx()} - {neighbour_atom.GetIdx()}, charge to be allocated: {given_charge}"
                    )
                    break
            for neighbour_atom in ob.OBAtomAtomIter(atom.OBAtom):
                if (
                    neighbour_atom.GetAtomicNum() in (6,)
                    and atom.OBAtom.GetBond(neighbour_atom).GetBondOrder() > 1
                    and atom.OBAtom.GetExplicitValence() == 5
                ):
                    atom.OBAtom.GetBond(neighbour_atom).SetBondOrder(
                        atom.OBAtom.GetBond(neighbour_atom).GetBondOrder() - 1
                    )
                    logger.debug(
                        f"{DEBUG_TAG} Fix over bonded C: {atom.OBAtom.GetIdx()} - {neighbour_atom.GetIdx()}, charge to be allocated: {given_charge}"
                    )
                    break
    return given_charge


def fix_fake_dipole(given_charge: int, omol: pybel.Molecule):
    """
    *This function can not be used separately.*

    Parameters:
        given_charge (int): The given charge to be fixed.
        omol (OBMol): The molecule object.
    """
    if given_charge >= 0:
        for atom in omol.atoms:
            if (
                atom.atomicnum == 7
                and atom.OBAtom.GetExplicitValence() == 4
                and atom.OBAtom.GetFormalCharge() == 0
                and all(
                    atom.OBAtom.GetBond(neighbour_atom).GetBondOrder() == 2
                    for neighbour_atom in ob.OBAtomAtomIter(atom.OBAtom)
                )
            ):
                for neighbour_atom in ob.OBAtomAtomIter(atom.OBAtom):
                    if (
                        atom.atomicnum == 7
                        and atom.OBAtom.GetExplicitValence() == 4
                        and atom.OBAtom.GetFormalCharge() == 0
                        and all(
                            atom.OBAtom.GetBond(neighbour_atom).GetBondOrder() == 2
                            for neighbour_atom in ob.OBAtomAtomIter(atom.OBAtom)
                        )
                    ):
                        atom.OBAtom.GetBond(neighbour_atom).SetBondOrder(
                            atom.OBAtom.GetBond(neighbour_atom).GetBondOrder() - 1
                        )
                        logger.debug(
                            f"{DEBUG_TAG} Fix fake dipole: {atom.OBAtom.GetIdx()} - {neighbour_atom.GetIdx()}, charge to be allocated: {given_charge}"
                        )
    return given_charge


def fix_CN_in_doubt(given_charge: int, omol: pybel.Molecule):
    """
    *This function can not be used separately.*

    Parameters:
        given_charge (int): The given charge to be fixed.
        omol (OBMol): The molecule object.
    """
    CN_in_doubt = 0
    doubt_pair = []
    for atom in omol.atoms:
        if (
            atom.atomicnum == 6
            and atom.OBAtom.GetFormalCharge() == 0
            and atom.OBAtom.GetExplicitValence() == 4
        ):
            for neighbour_atom in ob.OBAtomAtomIter(atom.OBAtom):
                if (
                    atom.OBAtom.GetBond(neighbour_atom).GetBondOrder() > 1
                    and neighbour_atom.GetAtomicNum() == 7
                    and neighbour_atom.GetExplicitValence() == 4
                ):
                    CN_in_doubt += 1
                    doubt_pair.append((atom.OBAtom, neighbour_atom))
    if CN_in_doubt > 0:
        logger.debug(
            f"{DEBUG_TAG} Fixing CN in doubt, number of pairs: {CN_in_doubt}, charge to be allocated: {given_charge}"
        )

    if CN_in_doubt % 2 == 0 and CN_in_doubt > 0:
        for atom_1, atom_2 in doubt_pair[: CN_in_doubt // 2]:
            atom_1.SetFormalCharge(-1)
            atom_1.GetBond(atom_2).SetBondOrder(
                atom_1.GetBond(atom_2).GetBondOrder() - 1
            )
            given_charge += 1
            logger.debug(
                f"{DEBUG_TAG} Fix CN in doubt: {atom_1.GetIdx()} - {atom_2.GetIdx()}, charge to be allocated: {given_charge}"
            )
    if CN_in_doubt == 1 and given_charge == 0:
        for atom_1, atom_2 in doubt_pair:
            atom_1.GetBond(atom_2).SetBondOrder(
                atom_1.GetBond(atom_2).GetBondOrder() - 1
            )
            logger.debug(
                f"{DEBUG_TAG} Fix CN in doubt: {atom_1.GetIdx()} - {atom_2.GetIdx()}, charge to be allocated: {given_charge}"
            )
    return given_charge


def fix_over_bonded_N_in_possitive(given_charge: int, omol: pybel.Molecule):
    """
    This function fixes the over-bonded nitrogen atoms in a molecule with a positive charge.

    *This function can not be used separately.*

    Args:
        given_charge (int): The given charge of the molecule.
        omol (pybel.Molecule): The molecule object to be fixed.
    """
    if given_charge > 0:
        for atom in omol.atoms:
            if (
                atom.atomicnum == 7  # Check if the atom is nitrogen
                and atom.OBAtom.GetExplicitValence()
                == 4  # Check if the atom is over-bonded
                and atom.OBAtom.GetFormalCharge()
                == 0  # Check if the atom has no formal charge
            ):
                atom.OBAtom.SetFormalCharge(1)  # Set the formal charge of the atom to 1
                given_charge -= 1  # Decrease the given charge by 1
                logger.debug(
                    f"{DEBUG_TAG} Fix over-bonded N: {atom.OBAtom.GetIdx()}, charge to be allocated: {given_charge}"
                )
    return given_charge


def fix_convinced_possitive_N(given_charge: int, omol: pybel.Molecule):
    """
    This function fixes the convinced possitive N atom in the molecule.

    Conditions:
        - `[*]-[N](-[*])(-[*])-[*]>>[*]-[N+](-[*])(-[*])-[*]`, make a new possitive charge

    *This function can not be used separately.*

    Parameters:
        given_charge (int): The given charge to be fixed.
        omol (OBMol): The molecule object.
    """
    for atom in omol.atoms:
        if (
            atom.atomicnum in (7, 15)
            and atom.OBAtom.GetExplicitValence() == 4
            and atom.OBAtom.GetFormalCharge() == 0
            and all(
                neighbour_atom.GetBond(atom.OBAtom).GetBondOrder() == 1
                for neighbour_atom in ob.OBAtomAtomIter(atom.OBAtom)
            )
        ):
            for neighbour_atom in ob.OBAtomAtomIter(atom.OBAtom):
                if atom.OBAtom.GetDistance(neighbour_atom) >= 1.1 * (
                    pt.GetRcovalent(atom.OBAtom.GetAtomicNum())
                    + pt.GetRcovalent(neighbour_atom.GetAtomicNum())
                ):
                    omol.OBMol.DeleteBond(atom.OBAtom.GetBond(neighbour_atom))
                    logger.debug(
                        f"{DEBUG_TAG} Fix un-convinced possitive N: {atom.OBAtom.GetIdx()}, charge to be allocated: {given_charge}"
                    )
                    break
            else:
                atom.OBAtom.SetFormalCharge(1)
                given_charge -= 1
                logger.debug(
                    f"{DEBUG_TAG} Fix convinced possitive N: {atom.OBAtom.GetIdx()}, charge to be allocated: {given_charge}"
                )
    return given_charge


def fix_charge_on_metal(
    resonance: pybel.Molecule, charge: int, charge_to_be_allocated: int
):
    """
    If metal found, allocate all the positive charge to it. Suppose only one metal atom.
    This function fixes the charge on a metal atom in a resonance structure.

    *This function can not be used separately.*

    Parameters:
        resonance (Resonance): The resonance structure containing the metal atom.
        charge (int): The charge to be allocated.
        charge_to_be_allocated (int): The remaining charge to be allocated.
    """
    if charge > 0:
        logger.debug(f"{DEBUG_TAG} Fixing charge on metal")
        for atom in resonance.atoms:
            if charge_to_be_allocated <= 0:
                break
            if atom.OBAtom.IsMetal():
                atom.OBAtom.SetFormalCharge(
                    atom.OBAtom.GetFormalCharge() + charge_to_be_allocated
                )
                charge_to_be_allocated -= charge_to_be_allocated
                logger.debug(
                    f"{DEBUG_TAG} Fix charge on metal: {atom.OBAtom.GetIdx()}, charge to be allocated: {charge_to_be_allocated}"
                )
    return charge_to_be_allocated


def fix_over_bonded_heteroatom(
    resonance: pybel.Molecule, charge: int, charge_to_be_allocated: int
):
    """
    Step 2.2.2: If no metal found, try to find the heteroatom with positive charge. e.g. [N+]R4
    This function fixes the charge on a metal atom in a resonance structure.

    *This function can not be used separately.*

    Parameters:
        resonance (Resonance): The resonance structure containing the metal atom.
        charge (int): The charge to be allocated.
        charge_to_be_allocated (int): The remaining charge to be allocated.
    """
    if charge > 0:
        logger.debug(f"{DEBUG_TAG} Fixing over-bonded heteroatom")
        for atom in resonance.atoms:
            if charge_to_be_allocated <= 0:
                break
            if atom.atomicnum in HETEROATOM and atom.OBAtom.GetFormalCharge() == 0:
                over_valence = atom.OBAtom.GetExplicitValence() - pt.GetDefaultValence(
                    atom.atomicnum
                )
                if over_valence > 0:
                    atom.OBAtom.SetFormalCharge(over_valence)
                    charge_to_be_allocated -= over_valence
                    logger.debug(
                        f"{DEBUG_TAG} Fix over-bonded heteroatom: {atom.OBAtom.GetIdx()}, charge to be allocated: {charge_to_be_allocated}"
                    )
    return charge_to_be_allocated


def fix_unbonded_heteroatom(
    resonance: pybel.Molecule, charge: int, charge_to_be_allocated: int
):
    """
    If no metal found, try to find the heteroatom with positive charge.

    In some cases, OpenBabel will split the heteroatom and one substituent into a neutral heteroatom and a substituent with a free radical.
    Therefore, try to find the heteroatom with a free radical and make sure the distance is less (or close to equal) than covalent bond length.
    If found, allocate the positive charge to the heteroatom and build a single bond to connect the two atoms.

    *This function can not be used separately.*

    Parameters:
        resonance (Resonance): The resonance structure containing the metal atom.
        charge (int): The charge to be allocated.
        charge_to_be_allocated (int): The remaining charge to be allocated.
    """
    if charge > 0:
        logger.debug(f"{DEBUG_TAG} Fixing unbonded heteroatom")
        for atom in resonance.atoms:
            if charge_to_be_allocated <= 0:
                break
            if atom.atomicnum in HETEROATOM and atom.OBAtom.GetFormalCharge() == 0:
                radical_atoms = [
                    satom
                    for satom in resonance.atoms
                    if get_under_bonded_number(satom.OBAtom) >= 1
                    and satom.idx != atom.idx
                ]
                if len(radical_atoms):
                    distances = {
                        radical_atom: atom.OBAtom.GetDistance(radical_atom.OBAtom)
                        for radical_atom in radical_atoms
                    }
                    closet_radical_atom: pybel.Atom = min(distances, key=distances.get)
                    logger.debug(
                        f"{DEBUG_TAG} closet_radical_atoms: {atom.OBAtom.GetIdx()} with {closet_radical_atom.OBAtom.GetIdx()} in distance: {distances[closet_radical_atom]}"
                    )
                    if (
                        distances[closet_radical_atom]
                        <= 1.4
                        * (
                            pt.GetRcovalent(atom.atomicnum)
                            + pt.GetRcovalent(closet_radical_atom.atomicnum)
                        )
                        and atom.OBAtom.GetBond(closet_radical_atom.OBAtom) is None
                    ):
                        if resonance.OBMol.AddBond(
                            atom.idx, closet_radical_atom.idx, 1
                        ):
                            atom.OBAtom.SetFormalCharge(1)
                            charge_to_be_allocated -= 1
                            logger.debug(
                                f"{DEBUG_TAG} Fix unbonded heteroatom: {atom.OBAtom.GetIdx()} with {closet_radical_atom.OBAtom.GetIdx()}"
                            )
    return charge_to_be_allocated


def transform_alkene(resonance: pybel.Molecule):
    for atom in resonance.atoms:
        if atom.atomicnum == 6:
            if all(
                atom.OBAtom.GetBond(neighbour_atom).GetBondOrder() == 2
                for neighbour_atom in ob.OBAtomAtomIter(atom.OBAtom)
            ):
                for neighbour_atom in ob.OBAtomAtomIter(atom.OBAtom):
                    for neighbour_2_atom in ob.OBAtomAtomIter(neighbour_atom):
                        if get_under_bonded_number(neighbour_2_atom) == 1 and all(
                            atom.OBAtom.GetBond(neighbour_atom).GetBondOrder() == 2
                            for neighbour_atom in ob.OBAtomAtomIter(atom.OBAtom)
                        ):
                            atom.OBAtom.GetBond(neighbour_atom).SetBondOrder(
                                atom.OBAtom.GetBond(neighbour_atom).GetBondOrder() - 1
                            )
                            neighbour_atom.GetBond(neighbour_2_atom).SetBondOrder(
                                neighbour_2_atom.GetBond(neighbour_atom).GetBondOrder()
                                + 1
                            )


def attribute_radicals_to_d_orbitals(resonance: pybel.Molecule):
    for atom in resonance.atoms:
        if atom.OBAtom.IsMetal():
            for neighbour_atom in list(ob.OBAtomAtomIter(atom.OBAtom)):
                if neighbour_atom.GetAtomicNum() in HETEROATOM:
                    if any(
                        get_under_bonded_number(neighbour_2_atom) == 1
                        for neighbour_2_atom in ob.OBAtomAtomIter(neighbour_atom)
                        if neighbour_2_atom.GetIndex() != atom.OBAtom.GetIndex()
                    ) and atom.OBAtom.GetBond(neighbour_atom):
                        resonance.OBMol.DeleteBond(atom.OBAtom.GetBond(neighbour_atom))
                        logger.debug(
                            f"{DEBUG_TAG} Delete dative bond between {atom.OBAtom.GetIdx()} and {neighbour_atom.GetIdx()}"
                        )
        if atom.atomicnum == 16:
            for neighbour_atom in ob.OBAtomAtomIter(atom.OBAtom):
                if get_under_bonded_number(neighbour_atom) == 1:
                    atom.OBAtom.GetBond(neighbour_atom).SetBondOrder(
                        atom.OBAtom.GetBond(neighbour_atom).GetBondOrder() + 1
                    )


def fix_over_bonded_N_remove_bond(
    resonance: pybel.Molecule, charge_to_be_allocated: int
):
    for atom in resonance.atoms:
        if charge_to_be_allocated <= 0:
            break
        if (
            atom.atomicnum in (7,)
            and atom.OBAtom.GetFormalCharge() == 0
            and atom.OBAtom.GetExplicitValence() == 4
        ):
            for neighbour_atom in ob.OBAtomAtomIter(atom.OBAtom):
                if (
                    atom.OBAtom.GetBond(neighbour_atom).GetBondOrder() >= 2
                    and neighbour_atom.GetFormalCharge() == 0
                ):
                    atom.OBAtom.GetBond(neighbour_atom).SetBondOrder(
                        atom.OBAtom.GetBond(neighbour_atom).GetBondOrder() - 1
                    )
                    neighbour_atom.SetFormalCharge(-1)
                    charge_to_be_allocated -= 1
    return charge_to_be_allocated


def fix_under_bonded_with_negative(
    resonance: pybel.Molecule, charge_to_be_allocated: int
):
    # Step 2.3.1: Try to find the heteroatom with negative charge first.
    for atom in resonance.atoms:
        if charge_to_be_allocated <= 0:
            break
        if atom.atomicnum in HETEROATOM and atom.OBAtom.GetFormalCharge() == 0:
            under_valence = atom.OBAtom.GetExplicitValence() - pt.GetDefaultValence(
                atom.atomicnum
            )
            if under_valence < 0:
                atom.OBAtom.SetFormalCharge(under_valence)
                charge_to_be_allocated -= abs(under_valence)

            # Step 2.3.2: Try to find the heteroatom with negative charge first.
    for atom in resonance.atoms:
        if charge_to_be_allocated <= 0:
            break
        if (
            atom.atomicnum in (6,)
            and get_under_bonded_number(atom.OBAtom) == 1
            and atom.OBAtom.GetFormalCharge() == 0
        ):
            atom.OBAtom.SetFormalCharge(-1)
            charge_to_be_allocated -= 1
    for atom in resonance.atoms:
        if charge_to_be_allocated <= 0:
            break
        if (
            atom.atomicnum in (1,)
            and get_under_bonded_number(atom.OBAtom) == 1
            and atom.OBAtom.GetFormalCharge() == 0
        ):
            atom.OBAtom.SetFormalCharge(-1)
            charge_to_be_allocated -= 1
    for atom in resonance.atoms:
        if charge_to_be_allocated <= 0:
            break
        if (
            atom.atomicnum in (6,)
            and get_under_bonded_number(atom.OBAtom) > 1
            and atom.OBAtom.GetFormalCharge() == 0
        ):
            atom.OBAtom.SetFormalCharge(-1)
            charge_to_be_allocated -= 1
    return charge_to_be_allocated


def fix_under_bonded_with_positive(
    resonance: pybel.Molecule, charge_to_be_allocated: int
):
    # Step 2.2.4: If no heteroatom (with the neighboring free radical) found, allocate the positive charge to the carbon or hydrogen atom with free radical.
    # Explaination:
    # In this case, OpenBabel will consider the carbon atom lack of bonds, so it will set the radical to 1.
    # Thus, find the carbon atoms with free radical, set the charge to 1 and set the radical to 0.
    for atom in resonance.atoms:
        if charge_to_be_allocated <= 0:
            break
        if (
            atom.atomicnum in (1, 6)
            and get_under_bonded_number(atom.OBAtom) == 1
            and atom.OBAtom.GetFormalCharge() == 0
        ):
            atom.OBAtom.SetFormalCharge(1)
            charge_to_be_allocated -= 1
            # allocate postive charge to the other elements
    for atom in resonance.atoms:
        if charge_to_be_allocated <= 0:
            break
        if (
            get_under_bonded_number(atom.OBAtom) == 1
            and atom.OBAtom.GetFormalCharge() == 0
        ):
            atom.OBAtom.SetFormalCharge(1)
            charge_to_be_allocated -= 1
    return charge_to_be_allocated


def has_metal(omol: pybel.Molecule):
    for atom in omol.atoms:
        if atom.OBAtom.IsMetal():
            return True
    return False


def final_check(omol: pybel.Molecule, total_charge: int = 0, total_radical: int = 0):
    logger.debug(
        f"{DEBUG_TAG} Final check, target charge {total_charge}, target radical {total_radical}"
    )
    charge, radical = 0, 0
    for atom in omol.atoms:
        charge += atom.OBAtom.GetFormalCharge()
        atom_radical = get_under_bonded_number(atom.OBAtom)
        if atom_radical != 2:
            radical += abs(atom_radical)
        logger.debug(
            f"{DEBUG_TAG} {atom.OBAtom.GetIdx()} {atom.OBAtom.GetAtomicNum()} {atom.OBAtom.GetFormalCharge()} {atom_radical}"
        )
    logger.debug(f"{DEBUG_TAG} total charge {charge}, total radical {radical}")
    if has_metal(omol):
        if total_charge == charge:
            logger.debug(f"{DEBUG_TAG} final check passed")
            logger.debug(f"{DEBUG_TAG} {omol.write('mol2')}")
            return True
    else:
        if total_charge == charge and total_radical == radical:
            logger.debug(f"{DEBUG_TAG} final check passed")
            logger.debug(f"{DEBUG_TAG} {omol.write('mol2')}")
            return True
    logger.debug(f"{DEBUG_TAG} final check failed")
    return False


def split_charge(omol: pybel.Molecule, given_charge: int = 0, given_radical: int = 0):
    if (
        given_charge == 0
        and all(atom.OBAtom.GetFormalCharge() == 0 for atom in omol.atoms)
        and sum(get_under_bonded_number(atom.OBAtom) for atom in omol.atoms) > 2
    ):
        logger.debug(f"{DEBUG_TAG} Spliting charges")
        radical_atoms = [
            atom for atom in omol.atoms if get_under_bonded_number(atom.OBAtom)
        ]
        for i in range(len(radical_atoms) // 2 - abs(given_charge)):
            for atom in radical_atoms:
                if atom.atomicnum == 8 and atom.OBAtom.GetFormalCharge() == 0:
                    atom.OBAtom.SetFormalCharge(atom.OBAtom.GetFormalCharge() - 1)
                    given_charge += 1
                    break
        for i in range(len(radical_atoms) // 2 - abs(given_charge)):
            for atom in radical_atoms:
                if atom.atomicnum == 7 and atom.OBAtom.GetFormalCharge() == 0:
                    atom.OBAtom.SetFormalCharge(atom.OBAtom.GetFormalCharge() - 1)
                    given_charge += 1
                    break
        for i in range(len(radical_atoms) // 2 - abs(given_charge)):
            for atom in radical_atoms:
                if atom.atomicnum == 6 and atom.OBAtom.GetFormalCharge() == 0:
                    atom.OBAtom.SetFormalCharge(atom.OBAtom.GetFormalCharge() - 1)
                    given_charge += 1
                    break
    return given_charge


def fix_carbine_neighbor_heteroatom(omol: pybel.Molecule, given_charge: int = 0):
    for atom in omol.atoms:
        if get_under_bonded_number(atom.OBAtom) == 2:
            for neighbor in ob.OBAtomAtomIter(atom.OBAtom):
                if (
                    neighbor.GetAtomicNum() in HETEROATOM
                    and neighbor.GetFormalCharge() == 0
                    and get_under_bonded_number(neighbor) == 0
                ):
                    atom.OBAtom.GetBond(neighbor).SetBondOrder(
                        atom.OBAtom.GetBond(neighbor).GetBondOrder() + 1
                    )
                    atom.OBAtom.SetFormalCharge(atom.OBAtom.GetFormalCharge() - 1)
                    neighbor.SetFormalCharge(1)
                    logger.debug(
                        f"{DEBUG_TAG} fixing carbine neighbor: {atom.OBAtom.GetIdx()} and {neighbor.GetIdx()}, charge to be allocated: {given_charge}"
                    )
                    break
    return given_charge


def fix_carbine_neighbor_unsaturated(omol: pybel.Molecule):
    for atom in omol.atoms:
        if get_under_bonded_number(atom.OBAtom) == 2:
            for neighbor_1 in ob.OBAtomAtomIter(atom.OBAtom):
                for neighbor_2 in ob.OBAtomAtomIter(neighbor_1):
                    if neighbor_2.GetIdx() != atom.OBAtom.GetIdx():
                        if neighbor_1.GetBond(neighbor_2).GetBondOrder() > 1:
                            neighbor_1.GetBond(neighbor_2).SetBondOrder(
                                neighbor_1.GetBond(neighbor_2).GetBondOrder() - 1
                            )
                            neighbor_1.GetBond(atom.OBAtom).SetBondOrder(
                                neighbor_1.GetBond(atom.OBAtom).GetBondOrder() + 1
                            )
                            logger.debug(
                                f"{DEBUG_TAG} fixing carbine neighbor: {atom.OBAtom.GetIdx()} and {neighbor_1.GetIdx()} {neighbor_2.GetIdx()}"
                            )


def polarisation(omol: pybel.Molecule, given_charge: int = 0, given_radical: int = 0):
    radicals = [atom for atom in omol.atoms if get_under_bonded_number(atom.OBAtom)]
    single_radicals = [
        atom for atom in radicals if get_under_bonded_number(atom.OBAtom) == 1
    ]
    if given_charge == 0 and given_radical == 0:
        if len(single_radicals) == 2 and len(radicals) == 2:
            for atom_1, atom_2 in itertools.permutations(single_radicals):
                if atom_1.atomicnum in HETEROATOM and atom_2.atomicnum in HETEROATOM:
                    if HETEROATOM.index(atom_1.atomicnum) <= HETEROATOM.index(
                        atom_2.atomicnum
                    ):
                        atom_1.OBAtom.SetFormalCharge(-1)
                        atom_2.OBAtom.SetFormalCharge(1)
                        return
                elif (
                    atom_1.atomicnum in HETEROATOM
                    and atom_2.atomicnum not in HETEROATOM
                ):
                    atom_1.OBAtom.SetFormalCharge(-1)
                    atom_2.OBAtom.SetFormalCharge(1)
                    return


def xyz_block_to_omol(
    xyz_block: str,
    total_charge: int = 0,
    total_radical: int = 0,
    greed_search=True,
    return_all_resonances=False,
) -> pybel.Molecule:
    """
    Convert XYZ block to pybel molecule object.
    Overcome the lack of original Openbabel function that can not consider charges.

    Parameters:
        xyz_block (str): The XYZ block to be converted.
        total_charge (int): The given charge to be fixed.
        total_radical (int): The given radical to be fixed.
        greed_search (bool): Whether to use greedy search the resonances to find the best structure.
        return_all_resonances (bool): Whether to return all resonances.

    Returns:
        The pybel molecule object.
    """
    given_charge = total_charge
    if abs(given_charge) > 3:
        raise ValueError("Charge must be between -3 and 3")

    logger.debug(f"{DEBUG_TAG} charge: {given_charge}, radicals_num: {total_radical}")

    # Use openbabel to initialize molecule without charge and radical.
    # OpenBabel can use XYZ file to recover a molecule with proper bonds.
    # If the molecule is a neutral molecule, the recovery will almostly always be true (dipole except).
    # Howerver, the lack is that OpenBabel does not consider the charge and radical information.
    # Therfore, the steps following will try to fix the charge and radical information.
    omol = pybel.readstring("xyz", xyz_block)
    logger.debug(f"{DEBUG_TAG} omol smiles: {omol.write('smi')}")
    time_start = time.time()
    # N-BCP
    fix_N_BCP(omol)
    fix_cyclobutylamine(omol)
    given_charge = fix_over_bonded_C(given_charge, omol)
    given_charge = fix_fake_dipole(given_charge, omol)
    given_charge = fix_CN_in_doubt(given_charge, omol)
    given_charge = fix_over_bonded_N_in_possitive(given_charge, omol)
    given_charge = fix_convinced_possitive_N(given_charge, omol)
    fix_carbine_neighbor_unsaturated(omol)
    given_charge = fix_carbine_neighbor_heteroatom(omol, given_charge)
    break_one_bond(omol, given_charge, total_radical)
    attribute_radicals_to_d_orbitals(omol)
    clean_neighbor_radicals(omol, given_charge, total_radical)
    given_charge = fix_carbine_neighbor_heteroatom(omol, given_charge)
    given_charge = fix_over_bonded_N_in_possitive(given_charge, omol)
    given_charge = fix_convinced_possitive_N(given_charge, omol)

    if final_check(omol, total_charge, total_radical):
        return omol
    given_charge = split_charge(omol, given_charge, total_radical)
    omol.OBMol.MakeDativeBonds()

    if final_check(omol, total_charge, total_radical):
        return omol

    if greed_search:
        possible_resonances = get_radical_resonances(omol)
    else:
        possible_resonances = [omol]

    logger.debug(f"{DEBUG_TAG} possible_resonances: {len(possible_resonances)}")
    for resonance in possible_resonances:
        logger.debug(f"{DEBUG_TAG} resonance smiles: {resonance.write('smi')}")
    recovered_resonances = []
    for resonance in possible_resonances:
        time_point = time.time()
        if time_point - time_start > molopconfig.max_structure_recovery_time:
            break
        charge = given_charge
        clean_neighbor_radicals(resonance, given_charge, total_radical)
        logger.debug(
            f"{DEBUG_TAG} charge to be allocated: {charge}, cleaned resonance smiles: {resonance.write('smi')}"
        )

        if charge == 0:
            fix_dipole_type_a(resonance)
            fix_dipole_type_b(resonance)

        charge_to_be_allocated = abs(charge)
        clean_neighbor_radicals(resonance, given_charge, total_radical)
        charge_to_be_allocated = fix_charge_on_metal(
            resonance, charge, charge_to_be_allocated
        )
        charge_to_be_allocated = fix_over_bonded_heteroatom(
            resonance, charge, charge_to_be_allocated
        )
        charge_to_be_allocated = fix_unbonded_heteroatom(
            resonance, charge, charge_to_be_allocated
        )

        if charge > 0:
            charge_to_be_allocated = fix_under_bonded_with_positive(
                resonance, charge_to_be_allocated
            )

        fix_under_bonded_dipole(resonance)
        clean_neighbor_radicals(resonance, given_charge, total_radical)
        if charge < 0:
            charge_to_be_allocated = fix_under_bonded_with_negative(
                resonance, charge_to_be_allocated
            )
            charge_to_be_allocated = fix_over_bonded_N_remove_bond(
                resonance, charge_to_be_allocated
            )

        attribute_radicals_to_d_orbitals(resonance)
        fix_under_bonded_dipole(resonance)
        clean_neighbor_radicals(resonance, given_charge, total_radical)
        if charge >= 0:
            polarisation(resonance, charge_to_be_allocated, total_radical)
        else:
            polarisation(resonance, -charge_to_be_allocated, total_radical)
        transform_alkene(resonance)
        resonance.OBMol.MakeDativeBonds()
        if charge_to_be_allocated == 0 and final_check(
            resonance, total_charge, total_radical
        ):
            recovered_resonances.append(
                (clean_resonances(resonance), charge_to_be_allocated, total_radical)
            )

    if len(recovered_resonances) == 0:
        raise ValueError("No legal molecule resonance found")
    if return_all_resonances:
        return [item[0] for item in recovered_resonances]

    recovered_resonances.sort(key=omol_score)
    final_omol = recovered_resonances[0][0]
    return final_omol


def rdmol_check(
    rdmol: Chem.rdchem.Mol, total_atom_num, total_charge=0, total_radical=0
):
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
    if has_metal:
        return charge == total_charge
    else:
        return charge == total_charge and radical == total_radical


def omol_to_rdmol(
    omol: pybel.Molecule, total_charge=0, total_radical=0
) -> Chem.rdchem.Mol:
    total_atom_num = omol.OBMol.NumAtoms()
    rdmol = Chem.MolFromMol2Block(omol.write("mol2"), removeHs=False)
    if rdmol is not None:
        rdmol = Chem.AddHs(rdmol)
        if rdmol_check(rdmol, total_atom_num, total_charge, total_radical):
            return rdmol
    rdmols = Chem.MolsFromCDXML(omol.write("cdxml"), removeHs=False)
    if len(rdmols) > 0:
        rdmol = Chem.AddHs(rdmols[0])
        if rdmol_check(rdmol, total_atom_num, total_charge, total_radical):
            rdmol.ClearProp("CDXML_FRAG_ID")
            return rdmol
    rdmol = Chem.MolFromMolBlock(omol.write("sdf"), removeHs=False)
    if rdmol is not None:
        rdmol = Chem.AddHs(rdmol)
        if rdmol_check(rdmol, total_atom_num, total_charge, total_radical):
            return rdmol
    raise ValueError("Failed to convert omol to rdmol")


def xyz_block_to_rdmol(
    xyz_block: str, total_charge: int = 0, total_radical: int = 0, greed_search=True
):
    return omol_to_rdmol(
        xyz_block_to_omol(xyz_block, total_charge, total_radical, greed_search),
        total_charge,
        total_radical,
    )
