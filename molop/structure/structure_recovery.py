import itertools
from typing import List, Tuple

from openbabel import openbabel as ob
from openbabel import pybel
from rdkit import Chem

pt = Chem.GetPeriodicTable()
HETEROATOM = (7, 8, 9, 15, 16, 17, 35, 53)


def get_spin(atom: ob.OBAtom) -> int:
    if atom.GetAtomicNum() == 6:
        if atom.GetFormalCharge() > 0:
            return pt.GetDefaultValence(atom.GetAtomicNum()) - atom.GetExplicitValence()
        if atom.GetFormalCharge() <= 0:
            return (
                pt.GetDefaultValence(atom.GetAtomicNum())
                + atom.GetFormalCharge()
                - atom.GetExplicitValence()
            )
    elif atom.GetAtomicNum() in (7, 15):
        if atom.GetExplicitValence() == 5:
            return 0
    elif atom.GetAtomicNum() == 16:
        if atom.GetExplicitValence() in (4, 6):
            return 0
    return (
        pt.GetDefaultValence(atom.GetAtomicNum())
        + atom.GetFormalCharge()
        - atom.GetExplicitValence()
    )


def fix_dipole_type_a(mol: pybel.Molecule):
    spin_atoms = [atom for atom in mol.atoms if get_spin(atom.OBAtom) == 1]
    for atom_1, atom_2 in itertools.combinations(spin_atoms, 2):
        intersections = set(
            (atom.idx for atom in mol.atoms if atom_1.OBAtom.GetBond(atom.OBAtom))
        ) & set((atom.idx for atom in mol.atoms if atom_2.OBAtom.GetBond(atom.OBAtom)))
        if len(intersections) == 1:
            # atom idx in OpenBabel is 1-based
            center_idx = intersections.pop() - 1
            if (
                mol.atoms[center_idx].atomicnum not in (7, 8)
                or get_spin(atom_1.OBAtom) != 1
                or get_spin(atom_2.OBAtom) != 1
            ):
                continue
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


def fix_dipole_type_b(mol: pybel.Molecule):
    spin_1_atoms = [atom for atom in mol.atoms if get_spin(atom.OBAtom) == 1]
    spin_2_atoms = [atom for atom in mol.atoms if get_spin(atom.OBAtom) == 2]
    for (atom_1, atom_2), atom_3 in itertools.product(
        itertools.permutations(spin_1_atoms, 2), spin_2_atoms
    ):
        if atom_2.atomicnum not in (7, 15):
            continue
        if mol.OBMol.GetBond(atom_1.idx, atom_2.idx) and mol.OBMol.GetBond(
            atom_2.idx, atom_3.idx
        ):
            atom_1.OBAtom.SetFormalCharge(-1)
            atom_2.OBAtom.SetFormalCharge(1)
            mol.OBMol.GetBond(atom_2.idx, atom_3.idx).SetBondOrder(3)
    return mol


def xyz_block_to_omol(xyz_block: str, charge: int = 0, spin: int = 0, check_spin=True):
    if abs(charge) > 2:
        raise ValueError("Charge must be between -2 and 2")
    if spin > 2 or spin < 0:
        raise ValueError("Spin must be between 0 and 2")

    # Step 1: Use openbabel to initialize molecule without charge and spin.
    # OpenBabel can use XYZ file to recover a molecule with proper bonds.
    # If the molecule is a neutral molecule, the recovery will almostly always be true (dipole except).
    # Howerver, the lack is that OpenBabel does not consider the charge and spin information.
    # Therfore, the steps following will try to fix the charge and spin information.
    omol = pybel.readstring("xyz", xyz_block)

    # Step 2: Process the molecule to allocate the charge.
    # The allocation will consume the number `charge_to_be_allocated`, and may reduce the spin in molecule.
    # Inner salts other than dipoles are not considered.
    # That means, you will not get molecules that have both positive and negative charges, except for the dipole.
    # Anyway, if charges have been all allocated, the spin remained will be used in the next step.

    # Step 2.1: Whatever the charge is, there is possibility that the dipole exists. Find them first.
    # Hope not to have dipoles along with other charges.

    # Step 2.1.1(Type A): Dipole like [CH2]=[O+]-[NH-].
    # OpenBabel will set the central atom to be neutral and the other two atoms neutral with spin 1.
    # Thus, find the combination of the two atoms with spin 1 and one heteroatom between them with spin 0.
    # Negative charges are mutually resonant at any position of the atoms on either side.
    # So it is straightforward to set one of the atoms as neutral and the other as charge -1.
    # Known issues: Molecule like `C=[O+][N-]N[CH-]C(=O)C` can not distinguish between the part CON and NNC.
    # To maintain the robustness of the script, the script will only process the first ternary it recognizes as a dipole.
    # This type of case is rare in real-world. Hope not to see it. :(
    omol = fix_dipole_type_a(omol)

    # Step 2.1.2(Type B): Dipole like [CH]#[N+]-[O-]. This type only allow N (maybe P) to be the center.
    # OpenBabel will set the central atom N (maybe P) with spin 1, the neutral atom with spin 2, and the negative atom with spin 1.
    # Thus, find the combination of the three atoms follow rule above.
    # Negative charges are mutually resonant at any position of the atoms on either side.
    # Nevertheless, I pact that atoms with spin 1 carry a formal charge -1.
    omol = fix_dipole_type_b(omol)

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
            charge += 1

    for atom in omol.atoms:
        if (
            atom.atomicnum == 7
            and atom.OBAtom.GetExplicitValence() == 4
            and atom.OBAtom.GetFormalCharge() == 0
            and all(
                atom.OBAtom.GetBond(neighbour_atom).GetBondOrder() == 1
                for neighbour_atom in ob.OBAtomAtomIter(atom.OBAtom)
            )
        ):
            atom.OBAtom.SetFormalCharge(1)
            charge -= 1

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
                ):
                    atom.OBAtom.GetBond(neighbour_atom).SetBondOrder(
                        atom.OBAtom.GetBond(neighbour_atom).GetBondOrder() - 1
                    )
                    neighbour_atom.SetFormalCharge(-1)
                    charge += 1
                    break

    charge_to_be_allocated = abs(charge)

    if charge > 0:
        # Step 2.2.1: If metal found, allocate all the positive charge to it.
        # Suppose only one metal atom.
        for atom in omol.atoms:
            if charge_to_be_allocated <= 0:
                break
            if atom.OBAtom.IsMetal():
                atom.OBAtom.SetFormalCharge(
                    atom.OBAtom.GetFormalCharge() + charge_to_be_allocated
                )
                charge_to_be_allocated -= charge_to_be_allocated

        # Step 2.2.2: If no metal found, try to find the heteroatom with positive charge. e.g. [N+]R4
        for atom in omol.atoms:
            if charge_to_be_allocated <= 0:
                break
            if atom.atomicnum in HETEROATOM and atom.OBAtom.GetFormalCharge() == 0:
                over_valence = atom.OBAtom.GetExplicitValence() - pt.GetDefaultValence(
                    atom.atomicnum
                )
                if over_valence > 0:
                    atom.OBAtom.SetFormalCharge(over_valence)
                    charge_to_be_allocated -= over_valence

        # Step 2.2.3: If no metal found, try to find the heteroatom with positive charge.
        # Explaination:
        # In some cases, OpenBabel will split the heteroatom and one substituent into a neutral heteroatom and a substituent with a free radical.
        # Therefore, try to find the heteroatom with a free radical and make sure the distance is less (or close to equal) than covalent bond length.
        # If found, allocate the positive charge to the heteroatom and build a single bond to connect the two atoms.
        for atom in omol.atoms:
            if charge_to_be_allocated <= 0:
                break
            if atom.atomicnum in HETEROATOM and atom.OBAtom.GetFormalCharge() == 0:
                spin_atoms = [
                    satom for satom in omol.atoms if get_spin(satom.OBAtom) == 1
                ]
                if len(spin_atoms):
                    distances = {
                        spin_atom: atom.OBAtom.GetDistance(spin_atom.OBAtom)
                        for spin_atom in spin_atoms
                    }
                    closet_spin_atom: pybel.Atom = min(distances, key=distances.get)
                    if distances[closet_spin_atom] <= 1.05 * (
                        pt.GetRcovalent(atom.atomicnum)
                        + pt.GetRcovalent(closet_spin_atom.atomicnum)
                    ):
                        omol.OBMol.AddBond(atom.idx, closet_spin_atom.idx, 1)
                        atom.OBAtom.SetFormalCharge(1)
                        charge_to_be_allocated -= 1

        # Step 2.2.4: If no heteroatom (with the neighboring free radical) found, allocate the positive charge to the carbon or hydrogen atom with free radical.
        # Explaination:
        # In this case, OpenBabel will consider the carbon atom lack of bonds, so it will set the spin to 1.
        # Thus, find the carbon atoms with free radical, set the charge to 1 and set the spin to 0.
        for atom in omol.atoms:
            if charge_to_be_allocated <= 0:
                break
            if (
                atom.atomicnum in (1, 6)
                and get_spin(atom.OBAtom) == 1
                and atom.OBAtom.GetFormalCharge() == 0
            ):
                atom.OBAtom.SetFormalCharge(1)
                charge_to_be_allocated -= 1

    if charge < 0:
        # Step 2.3.1: Try to find the heteroatom with negative charge first.
        for atom in omol.atoms:
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
        for atom in omol.atoms:
            if charge_to_be_allocated <= 0:
                break
            if (
                atom.atomicnum in (1, 6)
                and get_spin(atom.OBAtom) == 1
                and atom.OBAtom.GetFormalCharge() == 0
            ):
                atom.OBAtom.SetFormalCharge(-1)
                charge_to_be_allocated -= 1

    for atom in omol.atoms:
        if atom.OBAtom.IsMetal():
            for neighbour_atom in ob.OBAtomAtomIter(atom.OBAtom):
                if neighbour_atom.GetAtomicNum() in HETEROATOM:
                    if any(
                        get_spin(neighbour_2_atom)
                        for neighbour_2_atom in ob.OBAtomAtomIter(neighbour_atom)
                        if neighbour_2_atom.GetIndex() != atom.OBAtom.GetIndex()
                    ):
                        omol.OBMol.DeleteBond(atom.OBAtom.GetBond(neighbour_atom))
        if atom.atomicnum == 16:
            for neighbour_atom in ob.OBAtomAtomIter(atom.OBAtom):
                if get_spin(neighbour_atom) == 1:
                    atom.OBAtom.GetBond(neighbour_atom).SetBondOrder(
                        atom.OBAtom.GetBond(neighbour_atom).GetBondOrder() + 1
                    )

    spin_atoms = [satom for satom in omol.atoms if get_spin(satom.OBAtom) == 1]
    for spin_atom_1, spin_atom_2 in itertools.combinations(spin_atoms, 2):
        if spin_atom_1.OBAtom.GetBond(spin_atom_2.OBAtom):
            if (
                spin_atom_1.atomicnum == 6
                and spin_atom_1.OBAtom.GetExplicitValence() >= 4
            ):
                continue
            if (
                spin_atom_2.atomicnum == 6
                and spin_atom_2.OBAtom.GetExplicitValence() >= 4
            ):
                continue
            spin_atom_1.OBAtom.GetBond(spin_atom_2.OBAtom).SetBondOrder(
                spin_atom_1.OBAtom.GetBond(spin_atom_2.OBAtom).GetBondOrder() + 1,
            )
    omol.OBMol.MakeDativeBonds()
    totol_spin = sum(
        get_spin(atom.OBAtom) for atom in omol.atoms if not atom.OBAtom.IsMetal()
    )
    if check_spin:
        if charge_to_be_allocated > 0 or totol_spin != spin:
            raise ValueError(
                f"Charge {charge} and spin {spin} cannot be allocated to the molecule. Charge {abs(charge) - charge_to_be_allocated} and spin {totol_spin} found."
            )
    else:
        if charge_to_be_allocated > 0:
            raise ValueError(
                f"Charge {charge} cannot be allocated to the molecule. Charge {abs(charge) - charge_to_be_allocated} found."
            )
    return omol
