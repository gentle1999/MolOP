import itertools
from typing import List, Tuple
from molop.logger.logger import logger

from openbabel import openbabel as ob
from openbabel import pybel
from rdkit import Chem

pt = Chem.GetPeriodicTable()
HETEROATOM = (7, 8, 9, 15, 16, 17, 35, 53)


def get_spin(atom: ob.OBAtom) -> int:
    if atom.GetAtomicNum() == 6:
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


def clean_neighbor_spins(omol: pybel.Molecule):
    spin_atoms = [satom for satom in omol.atoms if get_spin(satom.OBAtom) >= 1]
    for spin_atom_1, spin_atom_2 in itertools.combinations(spin_atoms, 2):
        if (
            spin_atom_1.OBAtom.GetBond(spin_atom_2.OBAtom)
            and pt.GetDefaultValence(spin_atom_1.OBAtom.GetAtomicNum())
            - spin_atom_1.OBAtom.GetExplicitValence()
            > 0
            and pt.GetDefaultValence(spin_atom_2.OBAtom.GetAtomicNum())
            - spin_atom_2.OBAtom.GetExplicitValence()
            > 0
        ):
            spin_atom_1.OBAtom.GetBond(spin_atom_2.OBAtom).SetBondOrder(
                spin_atom_1.OBAtom.GetBond(spin_atom_2.OBAtom).GetBondOrder()
                + min(get_spin(spin_atom_1.OBAtom), get_spin(spin_atom_2.OBAtom)),
            )
        if (
            spin_atom_1.OBAtom.GetDistance(spin_atom_2.OBAtom)
            <= 1.3
            * (
                pt.GetRcovalent(spin_atom_1.atomicnum)
                + pt.GetRcovalent(spin_atom_2.atomicnum)
            )
            and spin_atom_1.OBAtom.GetBond(spin_atom_2.OBAtom) is None
        ):
            omol.OBMol.AddBond(
                spin_atom_1.OBAtom.GetIdx(), spin_atom_2.OBAtom.GetIdx(), 1
            )


def fix_under_bonded_dipole(omol: pybel.Molecule):
    for atom_1 in omol.atoms:
        if atom_1.OBAtom.GetFormalCharge() == 1:
            for atom_2 in ob.OBAtomAtomIter(atom_1.OBAtom):
                if (
                    atom_2.GetFormalCharge() == -1
                    and pt.GetDefaultValence(atom_1.OBAtom.GetAtomicNum())
                    - atom_1.OBAtom.GetExplicitValence()
                    > 0
                    and pt.GetDefaultValence(atom_2.GetAtomicNum())
                    - atom_2.GetExplicitValence()
                    > 0
                ):
                    atom_1.OBAtom.GetBond(atom_2).SetBondOrder(
                        atom_1.OBAtom.GetBond(atom_2).GetBondOrder() + 1
                    )
                    atom_1.OBAtom.SetFormalCharge(0)
                    atom_2.SetFormalCharge(0)


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
                or get_spin(mol.atoms[center_idx].OBAtom) != 0
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


def get_radical_resonances(omol: pybel.Molecule):
    spin_atoms = [satom for satom in omol.atoms if get_spin(satom.OBAtom) >= 1]
    resonances = [omol]
    for spin_atom in spin_atoms:
        for neighbour_1_atom in ob.OBAtomAtomIter(spin_atom.OBAtom):
            for neighbour_2_atom in ob.OBAtomAtomIter(neighbour_1_atom):
                if neighbour_2_atom.GetBond(neighbour_1_atom).GetBondOrder() >= 2:
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
    return resonances


def omol_score(omol_tuple: Tuple[pybel.Molecule, int]):
    score = 0
    score += 2 * sum(get_spin(atom.OBAtom) for atom in omol_tuple[0].atoms)
    score += sum(abs(atom.OBAtom.GetFormalCharge()) for atom in omol_tuple[0].atoms)
    return score


def clean_resonances_1(omol: pybel.Molecule):
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


def clean_resonances_2(omol: pybel.Molecule):
    smarts = pybel.Smarts("[#7+,#6+]=[*]-[*]=[*]-[#8-]")
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


def clean_resonances_3(omol: pybel.Molecule):
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


def clean_resonances_4(omol: pybel.Molecule):
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


def clean_resonances_5(omol: pybel.Molecule):
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


def clean_resonances(omol: pybel.Molecule):
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


def xyz_block_to_omol(
    xyz_block: str,
    given_charge: int = 0,
    given_spin: int = 0,
    greed_search=True,
):
    if abs(given_charge) > 3:
        raise ValueError("Charge must be between -3 and 3")
    if given_spin > 2 or given_spin < 0:
        raise ValueError("Spin must be between 0 and 2")
    
    logger.debug(f"charge: {given_charge}, spin: {given_spin}")

    # Step 1: Use openbabel to initialize molecule without charge and spin.
    # OpenBabel can use XYZ file to recover a molecule with proper bonds.
    # If the molecule is a neutral molecule, the recovery will almostly always be true (dipole except).
    # Howerver, the lack is that OpenBabel does not consider the charge and spin information.
    # Therfore, the steps following will try to fix the charge and spin information.
    omol = pybel.readstring("xyz", xyz_block)
    logger.debug(f"omol smiles: {omol.write('smi')}")

    # N-BCP
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

        for atom in omol.atoms:
            if (
                atom.atomicnum == 7
                and atom.OBAtom.GetExplicitValence() == 4
                and atom.OBAtom.GetFormalCharge() == 0
            ):
                atom.OBAtom.SetFormalCharge(1)
                given_charge -= 1

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
    omol.OBMol.MakeDativeBonds()
    clean_neighbor_spins(omol)

    if greed_search:
        possible_resonances = get_radical_resonances(omol)
    else:
        possible_resonances = [omol]

    recovered_resonances = []
    for resonance in possible_resonances:
        charge = given_charge
        spin = given_spin

        # Step 2: Process the molecule to allocate the charge.
        # The allocation will consume the number `charge_to_be_allocated`, and may reduce the spin in molecule.
        # Inner salts other than dipoles are not considered.
        # That means, you will not get molecules that have both positive and negative charges, except for the dipole.
        # Anyway, if charges have been all allocated, the spin remained will be used in the next step.

        # Step 2.1: Whatever the charge is, there is possibility that the dipole exists. Find them first.
        # Hope not to have dipoles along with other charges.

        if charge >= 0:
            # Step 2.1.1(Type A): Dipole like [CH2]=[O+]-[NH-].
            # OpenBabel will set the central atom to be neutral and the other two atoms neutral with spin 1.
            # Thus, find the combination of the two atoms with spin 1 and one heteroatom between them with spin 0.
            # Negative charges are mutually resonant at any position of the atoms on either side.
            # So it is straightforward to set one of the atoms as neutral and the other as charge -1.
            # Known issues: Molecule like `C=[O+][N-]N[CH-]C(=O)C` can not distinguish between the part CON and NNC.
            # To maintain the robustness of the script, the script will only process the first ternary it recognizes as a dipole.
            # This type of case is rare in real-world. Hope not to see it. :(
            resonance = fix_dipole_type_a(resonance)

            # Step 2.1.2(Type B): Dipole like [CH]#[N+]-[O-]. This type only allow N (maybe P) to be the center.
            # OpenBabel will set the central atom N (maybe P) with spin 1, the neutral atom with spin 2, and the negative atom with spin 1.
            # Thus, find the combination of the three atoms follow rule above.
            # Negative charges are mutually resonant at any position of the atoms on either side.
            # Nevertheless, I pact that atoms with spin 1 carry a formal charge -1.
            resonance = fix_dipole_type_b(resonance)

        charge_to_be_allocated = abs(charge)
        clean_neighbor_spins(resonance)

        if charge > 0:
            # Step 2.2.1: If metal found, allocate all the positive charge to it.
            # Suppose only one metal atom.
            for atom in resonance.atoms:
                if charge_to_be_allocated <= 0:
                    break
                if atom.OBAtom.IsMetal():
                    atom.OBAtom.SetFormalCharge(
                        atom.OBAtom.GetFormalCharge() + charge_to_be_allocated
                    )
                    charge_to_be_allocated -= charge_to_be_allocated

            # Step 2.2.2: If no metal found, try to find the heteroatom with positive charge. e.g. [N+]R4
            for atom in resonance.atoms:
                if charge_to_be_allocated <= 0:
                    break
                if atom.atomicnum in HETEROATOM and atom.OBAtom.GetFormalCharge() == 0:
                    over_valence = (
                        atom.OBAtom.GetExplicitValence()
                        - pt.GetDefaultValence(atom.atomicnum)
                    )
                    if over_valence > 0:
                        atom.OBAtom.SetFormalCharge(over_valence)
                        charge_to_be_allocated -= over_valence

            # Step 2.2.3: If no metal found, try to find the heteroatom with positive charge.
            # Explaination:
            # In some cases, OpenBabel will split the heteroatom and one substituent into a neutral heteroatom and a substituent with a free radical.
            # Therefore, try to find the heteroatom with a free radical and make sure the distance is less (or close to equal) than covalent bond length.
            # If found, allocate the positive charge to the heteroatom and build a single bond to connect the two atoms.
            for atom in resonance.atoms:
                if charge_to_be_allocated <= 0:
                    break
                if atom.atomicnum in HETEROATOM and atom.OBAtom.GetFormalCharge() == 0:
                    spin_atoms = [
                        satom
                        for satom in resonance.atoms
                        if get_spin(satom.OBAtom) == 1
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

            # Step 2.2.4: If no heteroatom (with the neighboring free radical) found, allocate the positive charge to the carbon or hydrogen atom with free radical.
            # Explaination:
            # In this case, OpenBabel will consider the carbon atom lack of bonds, so it will set the spin to 1.
            # Thus, find the carbon atoms with free radical, set the charge to 1 and set the spin to 0.
            for atom in resonance.atoms:
                if charge_to_be_allocated <= 0:
                    break
                if (
                    atom.atomicnum in (1, 6)
                    and get_spin(atom.OBAtom) == 1
                    and atom.OBAtom.GetFormalCharge() == 0
                ):
                    atom.OBAtom.SetFormalCharge(1)
                    charge_to_be_allocated -= 1

        fix_under_bonded_dipole(resonance)
        clean_neighbor_spins(resonance)
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
                    and get_spin(atom.OBAtom) == 1
                    and atom.OBAtom.GetFormalCharge() == 0
                ):
                    atom.OBAtom.SetFormalCharge(-1)
                    charge_to_be_allocated -= 1
            for atom in resonance.atoms:
                if charge_to_be_allocated <= 0:
                    break
                if (
                    atom.atomicnum in (1,)
                    and get_spin(atom.OBAtom) == 1
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
                            get_spin(neighbour_2_atom)
                            for neighbour_2_atom in ob.OBAtomAtomIter(neighbour_atom)
                            if neighbour_2_atom.GetIndex() != atom.OBAtom.GetIndex()
                        ) and atom.OBAtom.GetBond(neighbour_atom):
                            resonance.OBMol.DeleteBond(
                                atom.OBAtom.GetBond(neighbour_atom)
                            )
            if atom.atomicnum == 16:
                for neighbour_atom in ob.OBAtomAtomIter(atom.OBAtom):
                    if get_spin(neighbour_atom) == 1:
                        atom.OBAtom.GetBond(neighbour_atom).SetBondOrder(
                            atom.OBAtom.GetBond(neighbour_atom).GetBondOrder() + 1
                        )
        fix_under_bonded_dipole(resonance)
        clean_neighbor_spins(resonance)

        resonance.OBMol.MakeDativeBonds()
        if charge_to_be_allocated == 0:
            recovered_resonances.append((resonance, charge_to_be_allocated))

    if len(recovered_resonances) == 0:
        raise ValueError("No legal molecule resonance found")
    recovered_resonances.sort(key=omol_score)
    final_omol = recovered_resonances[0][0]
    return final_omol  # clean_resonances(final_omol)
