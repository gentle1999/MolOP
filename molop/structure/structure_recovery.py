"""
This module contains functions for recovering the structure of the molecule.

Offer a molecular graph recovery algorithm from the simple coodinates of atoms based on the initial work by [OpenBabel](https://openbabel.org/index.html), which can be easily used in the file reading process. This algorithm is different from the rdDetermineBonds (Original code implemented by [Jensen group](https://github.com/jensengroup/xyz2mol) and integrated in RDKit from the 2022.09 release, which is not suitable for the free radicals and complex containing metal.
  
Although our algorithm overcome the free radicals and metal problem and tested, it is still not perfect. There is no denying that, rdDetermineBonds works well for normal organic molecules. Thus, we would give molecule structure recovered by rdDetermineBonds first, if error happens, we will use our algorithm to recover the molecule structure instead. We hope that this strategy can take advantage of both approaches.
"""

import itertools
from typing import List, Tuple
from molop.logger.logger import logger

from openbabel import openbabel as ob
from openbabel import pybel
from rdkit import Chem

pt = Chem.GetPeriodicTable()
HETEROATOM = (7, 8, 9, 15, 16, 17, 35, 53)


def get_under_bonded_number(atom: ob.OBAtom) -> int:
    """
    OpenBabel determines that there are a large number of UNDER BONDED atoms in the molecule, which may be actual free radicals, positively or negatively charged.

    This function will finds out the UNDER BONDED number of the atom.

    Parameters:
        atom (ob.OBAtom): An object containing atomic information.

    Returns:
        The UNDER BONDED number of the atom.
    """
    if atom.GetAtomicNum() in (7, 15):
        if atom.GetExplicitValence() == 5:
            # For nitrogen and phosphorus atoms, the spin number is 0 when the explicit valence is 5
            return 0
    elif atom.GetAtomicNum() == 16:
        if atom.GetExplicitValence() in (4, 6):
            # For sulfur atoms, the spin number is 0 when the explicit valence is 4 or 6
            return 0
        if atom.GetExplicitValence() == 5:
            # For sulfur atoms, the spin number is 0 when the explicit valence is 4 or 6
            return 1
    # For other atoms, the spin number is the default valence electrons plus formal charge minus the explicit valence
    return (
        pt.GetDefaultValence(atom.GetAtomicNum())
        + atom.GetFormalCharge()
        - atom.GetExplicitValence()
    )


def clean_neighbor_radicals(omol: pybel.Molecule):
    """
    This function cleans up the neighbor radicals of a given molecule.

    Parameters:
        omol (pybel.Molecule): The molecule to clean.
    """
    radical_atoms = [
        ratom for ratom in omol.atoms if get_under_bonded_number(ratom.OBAtom) >= 1
    ]
    for radical_atom_1, radical_atom_2 in itertools.combinations(radical_atoms, 2):
        if (
            radical_atom_1.OBAtom.GetBond(
                radical_atom_2.OBAtom
            )  # the two atoms should be bonded
            # and pt.GetDefaultValence(radical_atom_1.OBAtom.GetAtomicNum())
            # - radical_atom_1.OBAtom.GetExplicitValence()
            # > 0
            # and pt.GetDefaultValence(radical_atom_2.OBAtom.GetAtomicNum())
            # - radical_atom_2.OBAtom.GetExplicitValence()
            # > 0
        ):
            """
            This code block checks if two spin atoms are bonded and if they have a valence greater than their explicit valence.
            If they do, it increases the bond order between them by the minimum spin value of the two atoms.
            """
            radical_atom_1.OBAtom.GetBond(radical_atom_2.OBAtom).SetBondOrder(
                radical_atom_1.OBAtom.GetBond(radical_atom_2.OBAtom).GetBondOrder()
                + min(
                    get_under_bonded_number(radical_atom_1.OBAtom),
                    get_under_bonded_number(radical_atom_2.OBAtom),
                ),
            )
        if (
            radical_atom_1.OBAtom.GetDistance(radical_atom_2.OBAtom)
            <= 1.3
            * (
                pt.GetRcovalent(radical_atom_1.atomicnum)
                + pt.GetRcovalent(radical_atom_2.atomicnum)
            )
            and radical_atom_1.OBAtom.GetBond(radical_atom_2.OBAtom) is None
        ):
            """
            This code block checks if two spin atoms are within a certain distance and if they are not already bonded.
            If they are, it adds a bond between them with order 1.
            """
            omol.OBMol.AddBond(
                radical_atom_1.OBAtom.GetIdx(), radical_atom_2.OBAtom.GetIdx(), 1
            )


def fix_under_bonded_dipole(omol: pybel.Molecule):
    """
    Sometimes, the neighbored atoms with posstive and negative charges are under bonded.

    e.g. `[H][C+]=[C-][H]>>[H][C]#[C][H]`

    This function fixes under-bonded dipole moments in a molecule.

    Parameters:
        omol (pybel.Molecule): The input molecule.
    """

    for atom_1 in omol.atoms:
        if atom_1.OBAtom.GetFormalCharge() == 1:
            for atom_2 in ob.OBAtomAtomIter(atom_1.OBAtom):
                if (
                    atom_2.GetFormalCharge()
                    == -1
                    # and pt.GetDefaultValence(atom_1.OBAtom.GetAtomicNum())
                    # - atom_1.OBAtom.GetExplicitValence()
                    # > 0
                    # and pt.GetDefaultValence(atom_2.GetAtomicNum())
                    # - atom_2.GetExplicitValence()
                    # > 0
                ):
                    atom_1.OBAtom.GetBond(atom_2).SetBondOrder(
                        atom_1.OBAtom.GetBond(atom_2).GetBondOrder() + 1
                    )
                    atom_1.OBAtom.SetFormalCharge(0)
                    atom_2.SetFormalCharge(0)


def fix_dipole_type_a(mol: pybel.Molecule) -> pybel.Molecule:
    """
    (Type A): Dipole like `[CH2]=[O+]-[NH-]`.

    OpenBabel will set the central atom to be neutral and the other two atoms neutral with spin 1 like `[CH2]-[O]-[NH]`.
    Thus, find the combination of the two atoms with spin 1 and one heteroatom between them with spin 0.
    Negative charges are mutually resonant at any position of the atoms on either side.
    So it is straightforward to set one of the atoms as neutral and the other as charge -1.

    Known issues: Molecule like `C=[O+][N-]N[CH-]C(=O)C` can not distinguish between the part CON and NNC.

    Parameters:
        mol (pybel.Molecule): The molecule to be fixed.

    Returns:
        pybel.Molecule: The fixed molecule.
    """

    # Get all atoms with spin 1
    spin_atoms = [
        atom for atom in mol.atoms if get_under_bonded_number(atom.OBAtom) == 1
    ]

    # Iterate through all combinations of spin atoms
    for atom_1, atom_2 in itertools.combinations(spin_atoms, 2):
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

            # Check the bond orders between the center atom and the spin atoms
            if (
                mol.OBMol.GetBond(atom_1.idx, center_idx + 1).GetBondOrder() == 1
                and mol.OBMol.GetBond(atom_2.idx, center_idx + 1).GetBondOrder() == 1
            ):
                mol.OBMol.GetBond(atom_1.idx, center_idx + 1).SetBondOrder(2)
                mol.atoms[center_idx].OBAtom.SetFormalCharge(1)
                atom_2.OBAtom.SetFormalCharge(-1)
            elif (
                mol.OBMol.GetBond(atom_1.idx, center_idx + 1).GetBondOrder() == 2
                and mol.OBMol.GetBond(atom_2.idx, center_idx + 1).GetBondOrder() == 1
            ):
                mol.OBMol.GetBond(atom_1.idx, center_idx + 1).SetBondOrder(3)
                mol.atoms[center_idx].OBAtom.SetFormalCharge(1)
                atom_2.OBAtom.SetFormalCharge(-1)
            elif (
                mol.OBMol.GetBond(atom_1.idx, center_idx + 1).GetBondOrder() == 1
                and mol.OBMol.GetBond(atom_2.idx, center_idx + 1).GetBondOrder() == 2
            ):
                mol.OBMol.GetBond(atom_2.idx, center_idx + 1).SetBondOrder(3)
                mol.atoms[center_idx].OBAtom.SetFormalCharge(1)
                atom_1.OBAtom.SetFormalCharge(-1)

    return mol


def fix_dipole_type_b(mol: pybel.Molecule) -> pybel.Molecule:
    """
    (Type B): Dipole like `[CH]#[N+]-[O-]`. This type only allow N (maybe P) to be the center.

    OpenBabel will set the central atom N (maybe P) with spin 1, the neutral atom with spin 2, and the negative atom with spin 1 like `[CH]-[N]-[O]`
    Thus, find the combination of the three atoms follow rule above.
    Negative charges are mutually resonant at any position of the atoms on either side.
    Nevertheless, I pact that atoms with spin 1 carry a formal charge -1.

    Parameters:
        mol (pybel.Molecule): The molecule to be fixed.

    Returns:
        The fixed molecule.
    """
    # Get spin 1 and spin 2 atoms
    spin_1_atoms = [
        atom for atom in mol.atoms if get_under_bonded_number(atom.OBAtom) == 1
    ]
    spin_2_atoms = [
        atom for atom in mol.atoms if get_under_bonded_number(atom.OBAtom) == 2
    ]

    # Iterate over all possible combinations of spin 1 atoms and spin 2 atoms
    for (atom_1, atom_2), atom_3 in itertools.product(
        itertools.permutations(spin_1_atoms, 2), spin_2_atoms
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

    return mol


def get_radical_resonances(omol: pybel.Molecule) -> List[pybel.Molecule]:
    """
    Retrieves a list of molecular structures that exhibit radical resonance.

    Parameters:
        omol (pybel.Molecule): The input molecule object

    Returns:
        A list of molecules with radical resonance structures
    """

    # Get a list of atoms with unpaired electrons (radicals)
    spin_atoms = [
        satom for satom in omol.atoms if get_under_bonded_number(satom.OBAtom) >= 1
    ]

    # Initialize the resonance molecules list with the input molecule
    resonances = [omol]

    # Iterate through the list of radical atoms
    for spin_atom in spin_atoms:
        # Iterate over neighboring atom 1 of the radical atom
        for neighbour_1_atom in ob.OBAtomAtomIter(spin_atom.OBAtom):
            # Iterate over neighboring atom 2 of neighboring atom 1
            for neighbour_2_atom in ob.OBAtomAtomIter(neighbour_1_atom):
                # If the bond order between neighbor 2 and neighbor 1 is greater than or equal to 2
                if neighbour_2_atom.GetBond(neighbour_1_atom).GetBondOrder() >= 2:
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
                        spin_atom.OBAtom.GetIdx(), neighbour_1_atom.GetIdx()
                    ).SetBondOrder(
                        new_omol.OBMol.GetBond(
                            spin_atom.OBAtom.GetIdx(), neighbour_1_atom.GetIdx()
                        ).GetBondOrder()
                        + 1
                    )
                    resonances.append(new_omol)

    # Return the list of resonance molecules
    return resonances


def omol_score(omol_tuple: Tuple[pybel.Molecule, int]) -> int:
    """
    Calculate the structural recovery score of a molecule, the lower the score the better the structural recovery.
    The criteria are 2 points for each unbonded electron and 1 point for each absolute value of charge.

    This function can only be used for comparison between isomers recovered from the same set of atomic coordinates.

    Args:
        omol_tuple (Tuple[pybel.Molecule, int]): A tuple containing a pybel.Molecule object and its charge.

    Returns:
        The structural recovery score of a molecule.
    """
    score = 0
    score += 2 * sum(
        get_under_bonded_number(atom.OBAtom) for atom in omol_tuple[0].atoms
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
    smarts = pybel.Smarts("[#7+]=[*]-[#6-,#8-]")
    res = smarts.findall(omol)
    while len(res):
        idxs = res[0]
        omol.OBMol.GetBond(idxs[1], idxs[2]).SetBondOrder(2)
        omol.OBMol.GetBond(idxs[0], idxs[1]).SetBondOrder(1)
        omol.atoms[idxs[0] - 1].OBAtom.SetFormalCharge(0)
        omol.atoms[idxs[-1] - 1].OBAtom.SetFormalCharge(0)
        res = smarts.findall(omol)
    return omol


def clean_resonances_4(omol: pybel.Molecule) -> pybel.Molecule:
    """
    `[#8]=[#6]-[#6-,#7-]>>[#8-]-[#6]=[#6,#7]`

    Parameters:
        omol (pybel.Molecule): The input molecule object

    Returns:
        The cleaned molecule.
    """
    smarts = pybel.Smarts("[#8]=[#6]-[#6-,#7-]")
    res = smarts.findall(omol)
    while len(res):
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
        idxs = res[0]
        omol.OBMol.GetBond(idxs[1], idxs[2]).SetBondOrder(3)
        omol.OBMol.GetBond(idxs[0], idxs[1]).SetBondOrder(1)
        omol.atoms[idxs[0] - 1].OBAtom.SetFormalCharge(-1)
        omol.atoms[idxs[-1] - 1].OBAtom.SetFormalCharge(0)
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
        clean_resonances_1,
        clean_resonances_2,
        clean_resonances_3,
        clean_resonances_4,
        clean_resonances_5,
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
                    break


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
            atom.atomicnum == 7
            and atom.OBAtom.GetExplicitValence() == 4
            and atom.OBAtom.GetFormalCharge() == 0
            and all(
                neighbour_atom.GetBond(atom.OBAtom).GetBondOrder() == 1
                for neighbour_atom in ob.OBAtomAtomIter(atom.OBAtom)
            )
        ):
            atom.OBAtom.SetFormalCharge(1)
            given_charge -= 1


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


def fix_CN_in_doubt(given_charge: int, omol: pybel.Molecule):
    """
    *This function can not be used separately.*

    Parameters:
        given_charge (int): The given charge to be fixed.
        omol (OBMol): The molecule object.
    """
    if given_charge >= 0:
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
        if CN_in_doubt % 2 == 0 and CN_in_doubt > 0:
            for atom_1, atom_2 in doubt_pair[: CN_in_doubt // 2]:
                atom_1.SetFormalCharge(-1)
                atom_1.GetBond(atom_2).SetBondOrder(
                    atom_1.GetBond(atom_2).GetBondOrder() - 1
                )
                given_charge += 1
        if CN_in_doubt == 1 and given_charge == 0:
            for atom_1, atom_2 in doubt_pair:
                atom_1.GetBond(atom_2).SetBondOrder(
                    atom_1.GetBond(atom_2).GetBondOrder() - 1
                )


def fix_over_bonded_N_in_possitive(given_charge: int, omol: pybel.Molecule):
    """
    This function fixes the over-bonded nitrogen atoms in a molecule with a positive charge.

    *This function can not be used separately.*

    Args:
        given_charge (int): The given charge of the molecule.
        omol (pybel.Molecule): The molecule object to be fixed.
    """
    if given_charge >= 0:
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
        for atom in resonance.atoms:
            if charge_to_be_allocated <= 0:
                break
            if atom.OBAtom.IsMetal():
                atom.OBAtom.SetFormalCharge(
                    atom.OBAtom.GetFormalCharge() + charge_to_be_allocated
                )
                charge_to_be_allocated -= charge_to_be_allocated


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


def fix_unbonded_heteroatom(resonance: pybel.Molecule, charge: int, charge_to_be_allocated: int):
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
        for atom in resonance.atoms:
            if charge_to_be_allocated <= 0:
                break
            if atom.atomicnum in HETEROATOM and atom.OBAtom.GetFormalCharge() == 0:
                spin_atoms = [
                    satom
                    for satom in resonance.atoms
                    if get_under_bonded_number(satom.OBAtom) == 1
                ]
                if len(spin_atoms):
                    distances = {
                        spin_atom: atom.OBAtom.GetDistance(spin_atom.OBAtom)
                        for spin_atom in spin_atoms
                    }
                    closet_spin_atom: pybel.Atom = min(distances, key=distances.get)
                    if (
                        distances[closet_spin_atom]
                        <= 1.3
                        * (
                            pt.GetRcovalent(atom.atomicnum)
                            + pt.GetRcovalent(closet_spin_atom.atomicnum)
                        )
                        and atom.OBAtom.GetBond(closet_spin_atom.OBAtom) is None
                    ):
                        resonance.OBMol.AddBond(atom.idx, closet_spin_atom.idx, 1)
                        atom.OBAtom.SetFormalCharge(1)
                        charge_to_be_allocated -= 1


def xyz_block_to_omol(
    xyz_block: str,
    given_charge: int = 0,
    greed_search=True,
) -> pybel.Molecule:
    """
    Convert XYZ block to pybel molecule object.
    Overcome the lack of original Openbabel function that can not consider charges.

    Parameters:
        xyz_block (str): The XYZ block to be converted.
        given_charge (int): The given charge to be fixed.
        greed_search (bool): Whether to use greedy search to find the best structure.

    Returns:
        The pybel molecule object.
    """
    if abs(given_charge) > 3:
        raise ValueError("Charge must be between -3 and 3")

    logger.debug(f"charge: {given_charge}")

    # Use openbabel to initialize molecule without charge and spin.
    # OpenBabel can use XYZ file to recover a molecule with proper bonds.
    # If the molecule is a neutral molecule, the recovery will almostly always be true (dipole except).
    # Howerver, the lack is that OpenBabel does not consider the charge and spin information.
    # Therfore, the steps following will try to fix the charge and spin information.
    omol = pybel.readstring("xyz", xyz_block)
    logger.debug(f"omol smiles: {omol.write('smi')}")

    # N-BCP
    fix_N_BCP(omol)
    fix_over_bonded_C(given_charge, omol)
    fix_fake_dipole(given_charge, omol)
    fix_CN_in_doubt(given_charge, omol)
    fix_over_bonded_N_in_possitive(given_charge, omol)
    fix_convinced_possitive_N(given_charge, omol)
    omol.OBMol.MakeDativeBonds()
    clean_neighbor_radicals(omol)

    if greed_search:
        possible_resonances = get_radical_resonances(omol)
    else:
        possible_resonances = [omol]

    recovered_resonances = []
    for resonance in possible_resonances:
        charge = given_charge

        if charge >= 0:
            resonance = fix_dipole_type_a(resonance)
            resonance = fix_dipole_type_b(resonance)

        charge_to_be_allocated = abs(charge)
        clean_neighbor_radicals(resonance)
        fix_charge_on_metal(resonance, charge, charge_to_be_allocated)
        fix_over_bonded_heteroatom(resonance, charge, charge_to_be_allocated)
        fix_unbonded_heteroatom(resonance, charge, charge_to_be_allocated)

        if charge > 0:
            # Step 2.2.4: If no heteroatom (with the neighboring free radical) found, allocate the positive charge to the carbon or hydrogen atom with free radical.
            # Explaination:
            # In this case, OpenBabel will consider the carbon atom lack of bonds, so it will set the spin to 1.
            # Thus, find the carbon atoms with free radical, set the charge to 1 and set the spin to 0.
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

        fix_under_bonded_dipole(resonance)
        clean_neighbor_radicals(resonance)
        if charge < 0:
            # Step 2.3.1: Try to find the heteroatom with negative charge first.
            for atom in resonance.atoms:
                if charge_to_be_allocated <= 0:
                    break
                if atom.atomicnum in HETEROATOM and atom.OBAtom.GetFormalCharge() == 0:
                    under_valence = (
                        atom.OBAtom.GetExplicitValence()
                        - pt.GetDefaultValence(atom.atomicnum)
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

        for atom in resonance.atoms:
            if atom.OBAtom.IsMetal():
                for neighbour_atom in list(ob.OBAtomAtomIter(atom.OBAtom)):
                    if neighbour_atom.GetAtomicNum() in HETEROATOM:
                        if any(
                            get_under_bonded_number(neighbour_2_atom)
                            for neighbour_2_atom in ob.OBAtomAtomIter(neighbour_atom)
                            if neighbour_2_atom.GetIndex() != atom.OBAtom.GetIndex()
                        ) and atom.OBAtom.GetBond(neighbour_atom):
                            resonance.OBMol.DeleteBond(
                                atom.OBAtom.GetBond(neighbour_atom)
                            )
            if atom.atomicnum == 16:
                for neighbour_atom in ob.OBAtomAtomIter(atom.OBAtom):
                    if get_under_bonded_number(neighbour_atom) == 1:
                        atom.OBAtom.GetBond(neighbour_atom).SetBondOrder(
                            atom.OBAtom.GetBond(neighbour_atom).GetBondOrder() + 1
                        )
        fix_under_bonded_dipole(resonance)
        clean_neighbor_radicals(resonance)

        resonance.OBMol.MakeDativeBonds()
        if charge_to_be_allocated == 0:
            recovered_resonances.append((resonance, charge_to_be_allocated))

    recovered_resonances = [item for item in recovered_resonances if item[1] == 0]
    if len(recovered_resonances) == 0:
        raise ValueError("No legal molecule resonance found")

    recovered_resonances.sort(key=omol_score)
    final_omol = recovered_resonances[0][0]
    return final_omol
