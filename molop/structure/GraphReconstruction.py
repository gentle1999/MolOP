"""
Author: TMJ
Date: 2025-01-15 23:01:22
LastEditors: TMJ
LastEditTime: 2025-11-03 00:08:47
Description: 请填写简介
"""

import itertools
from functools import lru_cache
from typing import List, Optional, Sequence, Tuple, Union

import numpy as np
import timeout_decorator
from openbabel import openbabel as ob
from openbabel import pybel
from rdkit import Chem
from rdkit.Geometry import Point3D

from molop.config import molopconfig, moloplogger
from molop.structure.FormatConverter import (
    omol_to_rdmol,
    omol_to_rdmol_by_graph,
    rdmol_to_omol,
    validate_omol,
)
from molop.structure.StructureTransformation import (
    make_dative_bonds,
    reset_atom_index,
)
from molop.structure.utils import HETEROATOM, estimate_bond_length, pt
from molop.utils.consts import (
    get_possible_metal_radicals,
    metal_valence_avialable_minor,
    metal_valence_avialable_prior,
)
from molop.utils.functions import is_metal

DEBUG_TAG = "[STRUCTURE RECOVERY]"


def get_under_bonded_number(atom: ob.OBAtom) -> int:
    """
    Get the number of atoms under the given atom.
    Suppose the atom is not a metal and follows the Octet rate (Suitable for most small organic molecules)

    Parameters:
        atom (ob.OBAtom): The atom to be checked.
    Returns:
        int: The number of atoms under the given atom.
    """
    atomic_number = atom.GetAtomicNum()
    if atomic_number <= 2:
        return (
            pt.GetDefaultValence(atomic_number)
            - atom.GetTotalValence()
            + atom.GetFormalCharge()
            - atom.GetSpinMultiplicity()
        )
    if pt.GetNOuterElecs(atomic_number) == 3 or (
        atom.GetFormalCharge() > 0
        and atom.GetTotalValence() < pt.GetDefaultValence(atomic_number)
    ):  # electron deficient
        return (
            pt.GetDefaultValence(atomic_number)
            - atom.GetTotalValence()
            - atom.GetFormalCharge()
            - atom.GetSpinMultiplicity()
        )
    return (
        pt.GetDefaultValence(atomic_number)
        - atom.GetTotalValence()
        + atom.GetFormalCharge()
        - atom.GetSpinMultiplicity()
    )  # electron abundant


def fresh_omol_charge_radical(omol: pybel.Molecule) -> pybel.Molecule:
    """
    Fresh the omol, set the spin multiplicity of radical atoms to 1.
    """
    for atom in omol.atoms:
        atom.OBAtom.SetSpinMultiplicity(0)
        if get_under_bonded_number(atom.OBAtom) < 0:
            # e.g. [B]R4 to [B-1]R4
            if pt.GetNOuterElecs(atom.OBAtom.GetAtomicNum()) == 3:
                atom.OBAtom.SetFormalCharge(get_under_bonded_number(atom.OBAtom))
            # e.g. [N]R4 to [N+1]R4
            else:
                atom.OBAtom.SetFormalCharge(-get_under_bonded_number(atom.OBAtom))
        elif get_under_bonded_number(atom.OBAtom):
            atom.OBAtom.SetSpinMultiplicity(get_under_bonded_number(atom.OBAtom))
    return omol


@lru_cache(maxsize=1000)
def xyz_to_separated_rwmol_no_metal(
    xyz_block: str, total_charge: int = 0, total_radical_electrons: int = 0
) -> Optional[Chem.rdchem.RWMol]:
    """
    Recover the structure of a molecule from its XYZ block.

    Suppose No metal in molecule

    Parameters:
        xyz_block (str): the XYZ block of the molecule
        total_charge (int): the total charge of the molecule
        total_radical_electrons (int): the total number of radical electrons of the molecule

    Returns:
        Optional[Chem.rdchem.RWMol]: the recovered RDKit molecule, or None if the recovery fails
    """
    if omol := xyz2omol(xyz_block, total_charge, total_radical_electrons):
        if rdmol := omol_to_rdmol(
            omol, total_charge=total_charge, total_radical=total_radical_electrons
        ):
            return Chem.RWMol(rdmol)
    return None


def combine_metal_with_rwmol(
    rwmol: Chem.rdchem.RWMol, metal_list: Sequence[Tuple[int, str, int, int, Point3D]]
) -> Chem.rdchem.RWMol:
    """
    Combine the metal atom with the rwmol.
    """
    sub_rwmol = Chem.RWMol(rwmol)
    for idx, symbol, valence, radical_num, position in metal_list:
        atom = Chem.Atom(symbol)
        atom.SetFormalCharge(valence)
        atom.SetNumRadicalElectrons(radical_num)
        atom.SetNoImplicit(True)
        metal_idx = sub_rwmol.AddAtom(atom)
        sub_rwmol.GetConformer().SetAtomPosition(metal_idx, position)
        reset_atom_index(
            sub_rwmol,
            list(range(0, idx)) + [metal_idx] + list(range(idx, metal_idx)),
        )
    return sub_rwmol


def xyz_to_separated_rwmol(
    xyz_block: str, total_charge: int = 0, total_radical_electrons: int = 0
) -> Optional[Chem.rdchem.RWMol]:
    """
    Convert the XYZ block of a molecule to a RDKit molecule, with the given total charge and radical electrons.

    Separated rwmol means that the metal atoms are set as cations or anions,
    and the remaining parts are reconstructed with the new total charge and radical electrons.

    Parameters:
        xyz_block (str): the XYZ block of the molecule
        total_charge (int): the total charge of the molecule
        total_radical_electrons (int): the total number of radical electrons of the molecule
    Returns:
        Optional[Chem.rdchem.RWMol]: the recovered RDKit molecule, or None if the recovery fails
    """
    rwmol = Chem.RWMol(Chem.MolFromXYZBlock(xyz_block))

    removable_metal_atom_idxs = [
        idx
        for idx in range(rwmol.GetNumAtoms())
        if is_metal(rwmol.GetAtomWithIdx(idx).GetAtomicNum())
        and pt.GetElementSymbol(rwmol.GetAtomWithIdx(idx).GetAtomicNum())
        in metal_valence_avialable_prior
    ]
    no_metal_rwmol = Chem.RWMol(rwmol)
    while idxs := [
        idx
        for idx in range(no_metal_rwmol.GetNumAtoms())
        if is_metal(no_metal_rwmol.GetAtomWithIdx(idx).GetAtomicNum())
        and pt.GetElementSymbol(no_metal_rwmol.GetAtomWithIdx(idx).GetAtomicNum())
        in metal_valence_avialable_prior
    ]:
        idx = idxs.pop(0)
        no_metal_rwmol.RemoveAtom(idx)
    no_metal_xyz_block = Chem.MolToXYZBlock(no_metal_rwmol)

    available_metal_valence_radical = [
        [
            (
                idx,
                rwmol.GetAtomWithIdx(idx).GetSymbol(),
                valence,
                radical_num,
                rwmol.GetConformer().GetAtomPosition(idx),
            )
            for valence in metal_valence_avialable_prior[
                rwmol.GetAtomWithIdx(idx).GetSymbol()
            ]
            + metal_valence_avialable_minor[rwmol.GetAtomWithIdx(idx).GetSymbol()]
            for radical_num in get_possible_metal_radicals(
                rwmol.GetAtomWithIdx(idx).GetSymbol(), valence
            )
        ]
        for idx in removable_metal_atom_idxs
    ]
    possible_metal_valence_radical_product = [
        product
        for product in itertools.product(*available_metal_valence_radical)
        if sum([radical_num for _, _, _, radical_num, _ in product])
        <= total_radical_electrons
    ]
    possible_rwmols = []
    for product in possible_metal_valence_radical_product:
        sub_rwmol = xyz_to_separated_rwmol_no_metal(
            no_metal_xyz_block,
            total_charge - sum([valence for _, _, valence, _, _ in product]),
            total_radical_electrons
            - sum([radical_num for _, _, _, radical_num, _ in product]),
        )
        if sub_rwmol is None:
            continue
        product_show = " | ".join(
            [
                f"{idx} {symbol} {valence} {radical_num}"
                for idx, symbol, valence, radical_num, position in product
            ]
        )
        moloplogger.debug(
            f"{DEBUG_TAG} | Trying Metal valence radical num with: {product_show}"
        )
        sub_rwmol = combine_metal_with_rwmol(sub_rwmol, product)
        possible_rwmols.append(sub_rwmol)
    if possible_rwmols:
        scores_pair = [(rwmol, structure_score(rwmol)) for rwmol in possible_rwmols]
        possible_rwmols = sorted(
            scores_pair,
            key=lambda x: x[1],
            reverse=True,
        )
        smiles_list = "\n".join(
            f"score = {score:.4f} | {Chem.CanonSmiles(Chem.MolToSmiles(rwmol))}"
            for rwmol, score in possible_rwmols
        )
        moloplogger.debug(
            f"{DEBUG_TAG} | Possible resonance structures: \n{smiles_list}"
        )
        return Chem.RWMol(possible_rwmols[0][0])
    return None


# Entrypoint for graph reconstruction
@timeout_decorator.timeout(molopconfig.max_structure_recovery_time)
def xyz_to_rdmol(
    xyz_block: str,
    total_charge: int = 0,
    total_radical_electrons: int = 0,
    make_dative: bool = True,
) -> Optional[Chem.rdchem.Mol]:
    """
    Convert the XYZ block of a molecule to a RDKit molecule, with the given total charge and radical electrons.

    Parameters:
        xyz_block (str):
            the XYZ block of the molecule
        total_charge (int):
            the total charge of the molecule
        total_radical_electrons (int):
            the total number of radical electrons of the molecule
        make_dative (bool):
            whether to make dative bonds or not
    Returns:
        Optional[Chem.rdchem.Mol]: the recovered RDKit molecule, or None if the recovery fails
    """
    if rwmol := xyz_to_separated_rwmol(
        xyz_block, total_charge, total_radical_electrons
    ):
        if make_dative:
            rwmol = make_dative_bonds(rwmol)
        moloplogger.debug(
            f"{DEBUG_TAG} | Final charge: {sum(atom.GetFormalCharge() for atom in rwmol.GetAtoms())} "
            f"total charge: {total_charge}; Final radical: {sum(atom.GetNumRadicalElectrons() for atom in rwmol.GetAtoms())}"
            f" total radical: {total_radical_electrons}; Final smiles: {Chem.MolToSmiles(rwmol)}"
        )
        Chem.SanitizeMol(rwmol, catchErrors=True)
        Chem.SetAromaticity(rwmol)
        Chem.DetectBondStereochemistry(rwmol)
        Chem.SetBondStereoFromDirections(rwmol)
        Chem.AssignStereochemistryFrom3D(rwmol)
        Chem.AssignCIPLabels(rwmol)
        return rwmol.GetMol()
    return None


@lru_cache(maxsize=1000)
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
    for atom in omol.atoms:
        atom.OBAtom.SetFormalCharge(0)  # Initialize the formal charge to 0
    omol = make_connections(omol)  # Connect the additional atoms if possible
    omol = pre_clean(omol)  # Clean the molecule before recovery
    moloplogger.debug(f"{DEBUG_TAG} | Input smiles: \n{omol.write('smi').strip()}")  # type: ignore

    # Allocate the posotive charge to the obviously charged atoms
    omol = fresh_omol_charge_radical(omol)
    given_charge = total_charge - sum(
        atom.OBAtom.GetFormalCharge() for atom in omol.atoms
    )

    omol, given_charge = eliminate_NNN(omol, given_charge)
    omol, given_charge = eliminate_high_positive_charge_atoms(omol, given_charge)
    omol, given_charge = eliminate_CN_in_doubt(omol, given_charge)
    omol, given_charge = eliminate_carboxyl(omol, given_charge)
    omol = clean_carbene_neighbor_unsaturated(omol)
    omol, given_charge = eliminate_carbene_neighbor_heteroatom(omol, given_charge)
    omol = clean_neighbor_radicals(omol)
    omol = clean_carbene_neighbor_unsaturated(omol)
    omol, given_charge = eliminate_charge_spliting(omol, given_charge)
    omol = break_deformed_ene(omol, given_charge, total_radical_electrons)
    omol, given_charge = break_one_bond(omol, given_charge, total_radical_electrons)

    moloplogger.debug(
        f"{DEBUG_TAG} | Given charge: {given_charge}, total charge: {total_charge}, actual charge: "
        f"{sum(atom.OBAtom.GetFormalCharge() for atom in omol.atoms)}, smiles: \n{omol.write('smi').strip()}"  # type: ignore
    )
    if validate_omol(omol, total_charge, total_radical_electrons):
        moloplogger.debug(
            f"{DEBUG_TAG} | Final charge: {sum(atom.OBAtom.GetFormalCharge() for atom in omol.atoms)} "
            f"total charge: {total_charge}; Final radical: {omol.OBMol.GetAtom(1).GetSpinMultiplicity()}"
            f" total radical: {total_radical_electrons}; Final smiles: \n{omol.write('smi').strip()}"  # type: ignore
        )
        return fresh_omol_charge_radical(omol)
    possible_resonances = get_radical_resonances(omol)
    moloplogger.debug(
        f"{DEBUG_TAG} | Possible resonance structures number: {len(possible_resonances)}"
    )
    recovered_resonances: List[Tuple[pybel.Molecule, int, int]] = []
    for idx, resonance in enumerate(possible_resonances):
        charge = given_charge
        moloplogger.debug(
            f"{DEBUG_TAG} | Resonance index: {idx}, recharge to be allocated: {charge}, smiles: \n{resonance.write('smi').strip()}"  # type: ignore
        )
        resonance, charge = eliminate_1_3_dipole(resonance, charge)
        resonance, charge = eliminate_positive_charges(resonance, charge)
        resonance, charge = eliminate_negative_charges(resonance, charge)
        resonance = clean_neighbor_radicals(resonance)
        moloplogger.debug(
            f"{DEBUG_TAG} | charge to be allocated: {charge}, pre-cleaned resonance smiles: \n{resonance.write('smi').strip()}"  # type: ignore
        )
        resonance = clean_resonances(clean_resonances(resonance))
        moloplogger.debug(
            f"{DEBUG_TAG} | charge to be allocated: {charge}, cleaned resonance smiles: \n{resonance.write('smi').strip()}"  # type: ignore
        )
        if validate_omol(resonance, total_charge, total_radical_electrons):
            recovered_resonances.append((resonance, charge, total_radical_electrons))
    if len(recovered_resonances) == 0:
        moloplogger.debug(
            f"{DEBUG_TAG} | Cannot recover the structure with given charge and radical"
        )
        return None

    recovered_resonances.sort(key=omol_score)
    scores_records = "\n".join(
        f"{omol_score(resonance):.4f} {resonance[0].write('smi').strip()}"  # type: ignore
        for resonance in recovered_resonances
    )
    moloplogger.debug(f"{DEBUG_TAG} | Scores with smiles: \n{scores_records}")
    final_omol = recovered_resonances[0][0]

    moloplogger.debug(
        f"{DEBUG_TAG} | Final charge: {sum(atom.OBAtom.GetFormalCharge() for atom in final_omol.atoms)} "
        f"total charge: {total_charge}; Final radical: {final_omol.OBMol.GetAtom(1).GetSpinMultiplicity()}"
        f" total radical: {total_radical_electrons}; Final smiles: \n{final_omol.write('smi').strip()}"  # type: ignore
    )
    return fresh_omol_charge_radical(final_omol)


def has_bridge(omol: pybel.Molecule, atom_idx_1: int, atom_idx_2: int) -> bool:
    atom_1 = omol.OBMol.GetAtom(atom_idx_1)
    atom_2 = omol.OBMol.GetAtom(atom_idx_2)
    neighbors_1 = [atom.GetIdx() for atom in ob.OBAtomAtomIter(atom_1)]
    neighbors_2 = [atom.GetIdx() for atom in ob.OBAtomAtomIter(atom_2)]
    for neighbor_1 in neighbors_1:
        if neighbor_1 in neighbors_2:
            return True
    return False


def make_connections(omol: pybel.Molecule, factor: float = 1.4) -> pybel.Molecule:
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
    donate_atoms: List[int] = list(itertools.chain(*donate_smarts.findall(omol)))
    accept_atoms: List[int] = list(itertools.chain(*accept_smarts.findall(omol)))

    for atom in donate_atoms:
        pairs = [(atom, accept_atom) for accept_atom in accept_atoms]
        pairs = sorted(
            pairs,
            key=lambda x: omol.OBMol.GetAtom(x[0]).GetDistance(
                omol.OBMol.GetAtom(x[1])
            ),
        )
        for pair_1, pair_2 in pairs:
            donate_atom = omol.OBMol.GetAtom(pair_1)
            accept_atom = omol.OBMol.GetAtom(pair_2)
            distance = donate_atom.GetDistance(accept_atom)
            if (
                distance
                < (
                    pt.GetRcovalent(donate_atom.GetAtomicNum())
                    + pt.GetRcovalent(accept_atom.GetAtomicNum())
                )
                * factor
                and pair_1 in donate_atoms
                and pair_2 in accept_atoms
                and pair_1 != pair_2
            ):
                if omol.OBMol.GetBond(pair_1, pair_2) is None:
                    omol.OBMol.AddBond(pair_1, pair_2, 1)
                    moloplogger.debug(
                        f"{DEBUG_TAG} | Add bond {pair_1} - {pair_2}, distance {distance}"
                    )
                    continue
                if omol.OBMol.GetBond(pair_1, pair_2).GetBondOrder() == 0:
                    omol.OBMol.GetBond(pair_1, pair_2).SetBondOrder(1)
                    moloplogger.debug(
                        f"{DEBUG_TAG} | Set bond order {pair_1} - {pair_2} to 1"
                    )
            donate_atoms: List[int] = list(
                itertools.chain(*donate_smarts.findall(omol))
            )
            accept_atoms: List[int] = list(
                itertools.chain(*accept_smarts.findall(omol))
            )
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

    smarts = pybel.Smarts("[S,P]=[*]")
    while res := smarts.findall(omol):
        idxs = res.pop(0)
        omol.OBMol.GetBond(idxs[0], idxs[1]).SetBondOrder(
            omol.OBMol.GetBond(idxs[0], idxs[1]).GetBondOrder() - 1
        )

    smarts = pybel.Smarts("[#6]1([#6]2)([#6]3)[#7]23[#6]1")
    while res := smarts.findall(omol):
        idxs = res.pop(0)
        for idx in idxs:
            indexs = set(idxs) - {idx}
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
            indexs = set(idxs) - {idx}
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

    return fresh_omol_charge_radical(omol)


def calculate_tetrahedron_volume(
    p1: Sequence[float], p2: Sequence[float], p3: Sequence[float], p4: Sequence[float]
) -> float:
    matrix = np.array([list(p1) + [1], list(p2) + [1], list(p3) + [1], list(p4) + [1]])
    determinant = np.linalg.det(matrix)
    return abs(determinant) / 6.0

def calculate_shape_quality(
    p1: Sequence[float], p2: Sequence[float], p3: Sequence[float], p4: Sequence[float]
) -> float:
    volume = calculate_tetrahedron_volume(p1, p2, p3, p4)
    _p1, _p2, _p3, _p4 = np.array(p1), np.array(p2), np.array(p3), np.array(p4)
    if np.isclose(volume, 0):
        return 0.0
    edges_sq = [
        np.linalg.norm(_p1 - _p2) ** 2,
        np.linalg.norm(_p1 - _p3) ** 2,
        np.linalg.norm(_p1 - _p4) ** 2,
        np.linalg.norm(_p2 - _p3) ** 2,
        np.linalg.norm(_p2 - _p4) ** 2,
        np.linalg.norm(_p3 - _p4) ** 2,
    ]
    sum_edges_sq = sum(edges_sq)
    l_rms_cubed = (sum_edges_sq / 6.0) ** 1.5
    NORMALIZATION_CONST = 6.0 * np.sqrt(2.0)
    quality = NORMALIZATION_CONST * (volume / l_rms_cubed)
    return max(0.0, min(1.0, quality))

def radical_deviation_scale(omol: pybel.Molecule, atom_idx: int) -> float:
    atom = omol.OBMol.GetAtom(atom_idx)
    neighbor_atoms = list(ob.OBAtomAtomIter(atom))
    if len(neighbor_atoms) == 2:
        return -abs(omol.OBMol.GetAngle(neighbor_atoms[0], atom, neighbor_atoms[1])-108) / 180
    if len(neighbor_atoms) == 3:
        return -10* calculate_shape_quality(
            [
                neighbor_atoms[0].GetVector().GetX(),
                neighbor_atoms[0].GetVector().GetY(),
                neighbor_atoms[0].GetVector().GetZ(),
            ],
            [
                neighbor_atoms[1].GetVector().GetX(),
                neighbor_atoms[1].GetVector().GetY(),
                neighbor_atoms[1].GetVector().GetZ(),
            ],
            [
                neighbor_atoms[2].GetVector().GetX(),
                neighbor_atoms[2].GetVector().GetY(),
                neighbor_atoms[2].GetVector().GetZ(),
            ],
            [atom.GetVector().GetX(), atom.GetVector().GetY(), atom.GetVector().GetZ()],
        )
    return 0.0


def negative_deviation_scale(omol: pybel.Molecule, atom_idx: int) -> float:
    return radical_deviation_scale(omol, atom_idx)


def positive_deviation_scale(omol: pybel.Molecule, atom_idx: int) -> float:
    return -radical_deviation_scale(omol, atom_idx)


def omol_score(omol_tuple: Tuple[pybel.Molecule, int, int]) -> float:
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
    omol, charge, radical_num = omol_tuple
    score = 0
    rdmol = Chem.MolFromSmiles(omol.write("smiles"))
    if rdmol is None:
        return 1000
    score -= len(rdmol.GetSubstructMatches(rdmol, uniquify=False))
    score += sum(
        radical_deviation_scale(omol, atom_idx)
        for atom_idx in range(1, omol.OBMol.NumAtoms() + 1)
        if omol.OBMol.GetAtom(atom_idx).GetSpinMultiplicity() > 0
    )
    score += sum(
        negative_deviation_scale(omol, atom_idx)
        for atom_idx in range(1, omol.OBMol.NumAtoms() + 1)
        if omol.OBMol.GetAtom(atom_idx).GetFormalCharge() < 0
    )
    score += sum(
        positive_deviation_scale(omol, atom_idx)
        for atom_idx in range(1, omol.OBMol.NumAtoms() + 1)
        if omol.OBMol.GetAtom(atom_idx).GetFormalCharge() > 0
    )
    score += 2 * (
        sum(atom.OBAtom.GetSpinMultiplicity() for atom in omol.atoms) - radical_num
    )
    score += sum(
        abs(atom.OBAtom.GetFormalCharge())
        for atom in omol.atoms
        if not is_metal(atom.OBAtom.GetAtomicNum())
    )
    return score


def structure_score(rwmol: Chem.rdchem.RWMol, ratio: float = 1.3) -> float:
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
    total_formal_charge = sum(
        abs(atom.GetFormalCharge())
        for atom in rwmol.GetAtoms()
        if not is_metal(atom.GetAtomicNum())
    )
    metal_charges = [
        atom.GetFormalCharge()
        for atom in rwmol.GetAtoms()
        if is_metal(atom.GetAtomicNum())
    ]
    negative_metal_charges = sum(charge for charge in metal_charges if charge < 0)
    total_metal_charge = sum(metal_charges)
    total_num_radical = sum(atom.GetNumRadicalElectrons() for atom in rwmol.GetAtoms())
    num_fragments = len(Chem.GetMolFrags(rwmol))
    metal_atoms = [atom for atom in rwmol.GetAtoms() if is_metal(atom.GetAtomicNum())]
    metal_score = np.std(metal_charges).item() if len(metal_charges) > 1 else 0
    dative_count = 0
    negative_atoms = [atom for atom in rwmol.GetAtoms() if atom.GetFormalCharge() < 0]
    for metal_atom, negative_atom in itertools.product(metal_atoms, negative_atoms):
        if rwmol.GetConformer().GetAtomPosition(negative_atom.GetIdx()).Distance(
            rwmol.GetConformer().GetAtomPosition(metal_atom.GetIdx())
        ) <= ratio * estimate_bond_length(
            rwmol.GetAtomWithIdx(metal_atom.GetIdx()).GetAtomicNum(),
            rwmol.GetAtomWithIdx(negative_atom.GetIdx()).GetAtomicNum(),
            Chem.rdchem.BondType.DATIVE,
        ):
            dative_count += 1

    return (
        10 * total_valence
        - 10 * total_formal_charge
        + 2 * total_metal_charge
        + 50 * negative_metal_charges
        - 10 * total_num_radical
        - 50 * num_fragments
        - 100 * metal_score
        + 2 * dative_count
    )


def clean_carbene_neighbor_unsaturated(omol: pybel.Molecule) -> pybel.Molecule:
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
        if atom1.GetSpinMultiplicity() == 2 and atom3.GetSpinMultiplicity() == 0:
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
    return fresh_omol_charge_radical(omol)


def eliminate_high_positive_charge_atoms(
    omol: pybel.Molecule, given_charge: int
) -> Tuple[pybel.Molecule, int]:
    smarts = pybel.Smarts("[*+1,*+2,*+3]-[Ov1+0,Nv2+0,Sv1+0]")
    while res := smarts.findall(omol):
        idxs = res.pop(0)
        if (
            -sum(
                atom.GetFormalCharge()
                for atom in ob.OBAtomAtomIter(omol.OBMol.GetAtom(idxs[0]))
            )
            >= omol.OBMol.GetAtom(idxs[0]).GetFormalCharge()
        ):
            break
        omol.OBMol.GetAtom(idxs[1]).SetFormalCharge(-1)
        moloplogger.debug(
            f"{DEBUG_TAG} | Eliminate atom {omol.OBMol.GetAtom(idxs[1]).GetAtomicNum()} {idxs[1]} with -1 charge"
        )
        given_charge += 1
    return fresh_omol_charge_radical(omol), given_charge


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
    return fresh_omol_charge_radical(omol), given_charge


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
    return fresh_omol_charge_radical(omol), given_charge


def eliminate_carbene_neighbor_heteroatom(
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
        if atom.OBAtom.GetSpinMultiplicity() == 2:
            for neibhor in ob.OBAtomAtomIter(atom.OBAtom):
                if neibhor.GetSpinMultiplicity():
                    return fresh_omol_charge_radical(omol), given_charge
            for neighbor in ob.OBAtomAtomIter(atom.OBAtom):
                if (
                    neighbor.GetAtomicNum() in HETEROATOM
                    and neighbor.GetFormalCharge() == 0
                    and neibhor.GetSpinMultiplicity() == 0
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
    return fresh_omol_charge_radical(omol), given_charge


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
    return fresh_omol_charge_radical(omol), given_charge


def break_deformed_ene(
    omol: pybel.Molecule, given_charge: int = 0, given_radical: int = 0, tolerance=5
):
    """
    Break deformed ene.

    This function checks if there are any unsatisfied valence states in the molecule,
    and based on the charge and radical conditions, finds and breaks an deformed ene.

    Parameters:
        omol (pybel.Molecule): The molecule to break bonds in.
        given_charge (int, optional): The allowed charge. Defaults to 0.
        given_radical (int, optional): The given radical. Defaults to 0.
        tolerance (int, optional): The tolerance for the torsion angle. Defaults to 5.
    """
    smarts = pybel.Smarts("[*]~[*+0]=,:[*+0]~[*]")
    res = list(smarts.findall(omol))
    while len(res):
        if (
            sum(atom.GetSpinMultiplicity() for atom in ob.OBMolAtomIter(omol.OBMol))
            >= abs(given_charge) + given_radical
        ):
            return fresh_omol_charge_radical(omol)
        idxs = res.pop(0)
        bond = omol.OBMol.GetBond(idxs[0], idxs[1])
        if bond.IsRotor() or bond.GetBondOrder() == 1:
            continue
        torsion_angle = abs(omol.OBMol.GetTorsion(*idxs))
        torsion_angle = min(torsion_angle, 180 - torsion_angle)
        if torsion_angle > tolerance:
            moloplogger.debug(
                f"{DEBUG_TAG} | torsion angle: {torsion_angle}, tolerance: {tolerance}"
            )
            moloplogger.debug(f"{DEBUG_TAG} | break bond {idxs[0]} and {idxs[1]}")
            bond.SetBondOrder(bond.GetBondOrder() - 1)
            bond.GetBeginAtom().SetSpinMultiplicity(
                bond.GetBeginAtom().GetSpinMultiplicity() + 1
            )
            bond.GetEndAtom().SetSpinMultiplicity(
                bond.GetEndAtom().GetSpinMultiplicity() + 1
            )
    smarts = pybel.Smarts("[*]~[*+0](=,:[*+0])~[*]")
    res = list(smarts.findall(omol))
    while len(res):
        if (
            sum(atom.GetSpinMultiplicity() for atom in ob.OBMolAtomIter(omol.OBMol))
            >= abs(given_charge) + given_radical
        ):
            return fresh_omol_charge_radical(omol)
        idxs = res.pop(0)
        bond = omol.OBMol.GetBond(idxs[0], idxs[1])
        if bond.IsRotor() or bond.GetBondOrder() == 1:
            continue
        torsion_angle = abs(omol.OBMol.GetTorsion(*idxs))
        torsion_angle = min(torsion_angle, 180 - torsion_angle)
        if torsion_angle > tolerance:
            moloplogger.debug(
                f"{DEBUG_TAG} | torsion angle: {torsion_angle}, tolerance: {tolerance}"
            )
            moloplogger.debug(f"{DEBUG_TAG} | break bond {idxs[0]} and {idxs[1]}")
            bond.SetBondOrder(bond.GetBondOrder() - 1)
            bond.GetBeginAtom().SetSpinMultiplicity(
                bond.GetBeginAtom().GetSpinMultiplicity() + 1
            )
            bond.GetEndAtom().SetSpinMultiplicity(
                bond.GetEndAtom().GetSpinMultiplicity() + 1
            )
    return fresh_omol_charge_radical(omol)


def break_one_bond(
    omol: pybel.Molecule, given_charge: int = 0, given_radical: int = 0
) -> Tuple[pybel.Molecule, int]:
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
            sum(atom.GetSpinMultiplicity() for atom in ob.OBMolAtomIter(omol.OBMol))
            >= abs(given_charge) + given_radical
        ):
            return fresh_omol_charge_radical(omol), given_charge
        # Get the indices of the first suitable bond's atoms
        idxs = res.pop(0)
        # Log debug information indicating the bond to be broken
        moloplogger.debug(
            f"{DEBUG_TAG} | break bond {idxs[0]} and {idxs[1]}: Double bond or triple bond"
        )
        # Get and reduce the bond order of the found bond
        bond = omol.OBMol.GetBond(idxs[0], idxs[1])
        bond.SetBondOrder(bond.GetBondOrder() - 1)
        bond.GetBeginAtom().SetSpinMultiplicity(
            bond.GetBeginAtom().GetSpinMultiplicity() + 1
        )
        bond.GetEndAtom().SetSpinMultiplicity(
            bond.GetEndAtom().GetSpinMultiplicity() + 1
        )
    smarts = pybel.Smarts("[#7+1,#15+1]=[*+0]")
    # Loop to find suitable bonds
    while res := smarts.findall(omol):
        if (
            sum(atom.GetSpinMultiplicity() for atom in ob.OBMolAtomIter(omol.OBMol))
            >= abs(given_charge) + given_radical
        ):
            return fresh_omol_charge_radical(omol), given_charge
        # Get the indices of the first suitable bond's atoms
        idxs = res.pop(0)
        # Log debug information indicating the bond to be broken
        moloplogger.debug(
            f"{DEBUG_TAG} | break bond {idxs[0]} and {idxs[1]}: *[N+](*)=[O]"
        )
        # Get and reduce the bond order of the found bond
        bond = omol.OBMol.GetBond(idxs[0], idxs[1])
        bond.SetBondOrder(bond.GetBondOrder() - 1)
        bond.GetBeginAtom().SetSpinMultiplicity(
            bond.GetBeginAtom().GetSpinMultiplicity() + 1
        )
        bond.GetEndAtom().SetSpinMultiplicity(
            bond.GetEndAtom().GetSpinMultiplicity() + 1
        )
        omol.OBMol.GetAtom(idxs[0]).SetFormalCharge(
            int(omol.OBMol.GetAtom(idxs[0]).GetFormalCharge() - 1)
        )
        given_charge += 1
    smarts = pybel.Smarts("[*+0]:[*+0]")
    # Loop to find suitable bonds, if only aromatic bonds are present
    while res := smarts.findall(omol):
        if (
            sum(atom.GetSpinMultiplicity() for atom in ob.OBMolAtomIter(omol.OBMol))
            >= abs(given_charge) + given_radical
        ):
            return fresh_omol_charge_radical(omol), given_charge
        # Get the indices of the first suitable bond's atoms
        idxs = res.pop(0)
        # Log debug information indicating the bond to be broken
        moloplogger.debug(f"{DEBUG_TAG} | break bond {idxs[0]} and {idxs[1]}: Aromatic")
        # Get and reduce the bond order of the found bond
        bond = omol.OBMol.GetBond(idxs[0], idxs[1])
        bond.SetBondOrder(bond.GetBondOrder() - 1)
        bond.GetBeginAtom().SetSpinMultiplicity(
            bond.GetBeginAtom().GetSpinMultiplicity() + 1
        )
        bond.GetEndAtom().SetSpinMultiplicity(
            bond.GetEndAtom().GetSpinMultiplicity() + 1
        )

    if all(bond.GetBondOrder() == 1 for bond in ob.OBMolBondIter(omol.OBMol)):
        for bond in ob.OBMolBondIter(omol.OBMol):
            if (
                sum(atom.GetSpinMultiplicity() for atom in ob.OBMolAtomIter(omol.OBMol))
                >= abs(given_charge) + given_radical
            ):
                return fresh_omol_charge_radical(omol), given_charge
            moloplogger.debug(
                f"{DEBUG_TAG} | break bond {bond.GetBeginAtom().GetIdx()} and {bond.GetEndAtom().GetIdx()}:"
                f" Single bond"
            )
            bond.GetBeginAtom().SetSpinMultiplicity(
                bond.GetBeginAtom().GetSpinMultiplicity() + 1
            )
            bond.GetEndAtom().SetSpinMultiplicity(
                bond.GetEndAtom().GetSpinMultiplicity() + 1
            )
            omol.OBMol.DeleteBond(bond)
    return fresh_omol_charge_radical(omol), given_charge


def clean_neighbor_radicals(omol: pybel.Molecule) -> pybel.Molecule:
    """
    Clean the neighbor radicals.
    Parameters:
        omol (pybel.Molecule): The molecule to be transformed.
    """
    for bond in ob.OBMolBondIter(omol.OBMol):
        if (
            bond.GetBeginAtom().GetSpinMultiplicity()
            and bond.GetEndAtom().GetSpinMultiplicity()
        ):
            bond_to_add = min(
                bond.GetBeginAtom().GetSpinMultiplicity(),
                bond.GetEndAtom().GetSpinMultiplicity(),
            )
            bond.SetBondOrder(bond.GetBondOrder() + bond_to_add)
            bond.GetBeginAtom().SetSpinMultiplicity(
                bond.GetBeginAtom().GetSpinMultiplicity() - bond_to_add
            )
            bond.GetEndAtom().SetSpinMultiplicity(
                bond.GetEndAtom().GetSpinMultiplicity() - bond_to_add
            )

            moloplogger.debug(
                f"{DEBUG_TAG} | Cleaning bond {bond.GetBeginAtom().GetIdx()} "
                f"and {bond.GetEndAtom().GetIdx()}"
            )
    return fresh_omol_charge_radical(omol)


def eliminate_charge_spliting(
    omol: pybel.Molecule, given_charge: int
) -> Tuple[pybel.Molecule, int]:
    if (
        all(atom.GetFormalCharge() == 0 for atom in ob.OBMolAtomIter(omol.OBMol))
        and sum(atom.GetSpinMultiplicity() for atom in ob.OBMolAtomIter(omol.OBMol))
        >= 2
    ):
        radical_atoms = [
            atom for atom in ob.OBMolAtomIter(omol.OBMol) if atom.GetSpinMultiplicity()
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
    return fresh_omol_charge_radical(omol), given_charge


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
        if omol.OBMol.GetAtom(idxs[0]).GetSpinMultiplicity() >= 1:
            new_omol = rdmol_to_omol(omol_to_rdmol_by_graph(omol))
            new_omol.OBMol.GetBond(idxs[1], idxs[2]).SetBondOrder(
                new_omol.OBMol.GetBond(idxs[1], idxs[2]).GetBondOrder() - 1
            )
            new_omol.OBMol.GetBond(idxs[0], idxs[1]).SetBondOrder(
                new_omol.OBMol.GetBond(idxs[0], idxs[1]).GetBondOrder() + 1
            )

            result.append(fresh_omol_charge_radical(new_omol))
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
            atom3.GetSpinMultiplicity()
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
    return fresh_omol_charge_radical(omol), given_charge


def eliminate_positive_charges(
    omol: pybel.Molecule, given_charge: int
) -> Tuple[pybel.Molecule, int]:
    smarts = pybel.Smarts("[Nv3+0]=[Nv2+0]")
    while given_charge > 0 and (res := smarts.findall(omol)):
        idxs = res.pop(0)
        omol.OBMol.GetAtom(idxs[1]).SetFormalCharge(1)
        given_charge -= 1
    smarts = pybel.Smarts("[#6v3+0,#6v2+0,#1v0+0]")
    while given_charge > 0 and (res := smarts.findall(omol)):
        idxs = res.pop(0)
        omol.OBMol.GetAtom(idxs[0]).SetFormalCharge(1)
        given_charge -= 1
    for atom in omol.atoms:
        if given_charge <= 0:
            break
        if (
            atom.OBAtom.GetSpinMultiplicity() >= 1
            and atom.OBAtom.GetFormalCharge() == 0
        ):
            to_add = min(atom.OBAtom.GetSpinMultiplicity(), given_charge)
            atom.OBAtom.SetFormalCharge(to_add)
            given_charge -= to_add
    return fresh_omol_charge_radical(omol), given_charge


def eliminate_negative_charges(
    omol: pybel.Molecule, given_charge: int
) -> Tuple[pybel.Molecule, int]:
    possibile_heteroatoms = []
    for atom in omol.atoms:
        if (
            atom.atomicnum in HETEROATOM
            and atom.OBAtom.GetFormalCharge() == 0
            and atom.OBAtom.GetSpinMultiplicity()
        ):
            possibile_heteroatoms.append((atom, HETEROATOM.index(atom.atomicnum)))
    possibile_heteroatoms.sort(key=lambda x: x[1])
    for atom, _ in possibile_heteroatoms:
        if given_charge >= 0:
            break
        to_add = min(atom.OBAtom.GetSpinMultiplicity(), abs(given_charge))
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
        to_add = min(
            omol.OBMol.GetAtom(idxs[0]).GetSpinMultiplicity(), abs(given_charge)
        )
        omol.OBMol.GetAtom(idxs[0]).SetFormalCharge(-to_add)
        given_charge += to_add
        moloplogger.debug(
            f"{DEBUG_TAG} | Eliminate negative charge: {atom.OBAtom.GetIdx()} "
            f"with charge {to_add}"
        )
    smarts = pybel.Smarts("[#1v0+0]")
    while given_charge < 0 and (res := smarts.findall(omol)):
        idxs = res.pop(0)
        to_add = min(
            omol.OBMol.GetAtom(idxs[0]).GetSpinMultiplicity(), abs(given_charge)
        )
        omol.OBMol.GetAtom(idxs[0]).SetFormalCharge(-to_add)
        given_charge += to_add
        moloplogger.debug(
            f"{DEBUG_TAG} | Eliminate negative charge: {atom.OBAtom.GetIdx()} "
            f"with charge {to_add}"
        )
    smarts = pybel.Smarts("[#6v2+0,#6v1+0,#6v0+0]")
    while given_charge < 0 and (res := smarts.findall(omol)):
        idxs = res.pop(0)
        to_add = min(
            omol.OBMol.GetAtom(idxs[0]).GetSpinMultiplicity(), abs(given_charge)
        )
        omol.OBMol.GetAtom(idxs[0]).SetFormalCharge(-to_add)
        given_charge += to_add
        moloplogger.debug(
            f"{DEBUG_TAG} | Eliminate negative charge: {atom.OBAtom.GetIdx()} "
            f"with charge {to_add}"
        )
    return fresh_omol_charge_radical(omol), given_charge


def clean_resonances_0(omol: pybel.Molecule) -> pybel.Molecule:
    smarts = pybel.Smarts("[*-]-[*]=[*]~[*+]")
    res = list(smarts.findall(omol))
    while len(res):
        idxs = res.pop(0)
        if (
            pt.GetDefaultValence(omol.OBMol.GetAtom(idxs[-1]).GetAtomicNum())
            > omol.OBMol.GetAtom(idxs[-1]).GetTotalValence()
        ):
            moloplogger.debug(f"{DEBUG_TAG} | Cleaning resonance 0: {idxs}")
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
    return fresh_omol_charge_radical(omol)


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
        moloplogger.debug(f"{DEBUG_TAG} | Cleaning resonance 1: {idxs}")
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
    return fresh_omol_charge_radical(omol)


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
        moloplogger.debug(f"{DEBUG_TAG} | Cleaning resonance 2: {idxs}")
        omol.OBMol.GetBond(idxs[4], idxs[5]).SetBondOrder(2)
        omol.OBMol.GetBond(idxs[3], idxs[4]).SetBondOrder(1)
        omol.OBMol.GetBond(idxs[1], idxs[3]).SetBondOrder(2)
        omol.OBMol.GetBond(idxs[0], idxs[1]).SetBondOrder(1)
        omol.OBMol.GetAtom(idxs[0]).SetFormalCharge(-1)
        omol.OBMol.GetAtom(idxs[-1]).SetFormalCharge(0)
    return fresh_omol_charge_radical(omol)


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
        moloplogger.debug(f"{DEBUG_TAG} | Cleaning resonance 3: {idxs}")
        omol.OBMol.GetBond(idxs[3], idxs[4]).SetBondOrder(2)
        omol.OBMol.GetBond(idxs[2], idxs[3]).SetBondOrder(1)
        omol.OBMol.GetBond(idxs[1], idxs[2]).SetBondOrder(2)
        omol.OBMol.GetBond(idxs[0], idxs[1]).SetBondOrder(1)
        omol.OBMol.GetAtom(idxs[0]).SetFormalCharge(0)
        omol.OBMol.GetAtom(idxs[-1]).SetFormalCharge(0)
    return fresh_omol_charge_radical(omol)


def clean_resonances_4(omol: pybel.Molecule) -> pybel.Molecule:
    """
    `[#7+,#8+]=[*]-[#6-,#7-,#8-]>>[#7,#8]-[*]=[#6,#7,#8]`

    Parameters:
        omol (pybel.Molecule): The input molecule object

    Returns:
        pybel.Molecule: The cleaned molecule.
    """
    smarts = pybel.Smarts("[#7+,#8+]=[*]-[#6-,#7-,#8-]")
    res = list(smarts.findall(omol))
    while len(res):
        idxs = res.pop(0)
        moloplogger.debug(f"{DEBUG_TAG} | Cleaning resonance 4: {idxs}")
        omol.OBMol.GetBond(idxs[1], idxs[2]).SetBondOrder(
            omol.OBMol.GetBond(idxs[1], idxs[2]).GetBondOrder() + 1
        )
        omol.OBMol.GetBond(idxs[0], idxs[1]).SetBondOrder(
            omol.OBMol.GetBond(idxs[0], idxs[1]).GetBondOrder() - 1
        )
        omol.OBMol.GetAtom(idxs[0]).SetFormalCharge(0)
        omol.OBMol.GetAtom(idxs[-1]).SetFormalCharge(0)
    return fresh_omol_charge_radical(omol)


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
        moloplogger.debug(f"{DEBUG_TAG} | Cleaning resonance 5: {idxs}")
        omol.OBMol.GetBond(idxs[1], idxs[2]).SetBondOrder(
            omol.OBMol.GetBond(idxs[1], idxs[2]).GetBondOrder() + 1
        )
        omol.OBMol.GetBond(idxs[0], idxs[1]).SetBondOrder(
            omol.OBMol.GetBond(idxs[0], idxs[1]).GetBondOrder() - 1
        )
        omol.OBMol.GetAtom(idxs[0]).SetFormalCharge(-1)
        omol.OBMol.GetAtom(idxs[-1]).SetFormalCharge(0)
    return fresh_omol_charge_radical(omol)


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
        moloplogger.debug(f"{DEBUG_TAG} | Cleaning resonance 6: {idxs}")
        omol.OBMol.GetBond(idxs[1], idxs[2]).SetBondOrder(3)
        omol.OBMol.GetBond(idxs[0], idxs[1]).SetBondOrder(1)
        omol.OBMol.GetAtom(idxs[0]).SetFormalCharge(-1)
        omol.OBMol.GetAtom(idxs[-1]).SetFormalCharge(0)
    return fresh_omol_charge_radical(omol)


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
        moloplogger.debug(f"{DEBUG_TAG} | Cleaning resonance 7: {idxs}")
        omol.OBMol.GetBond(idxs[1], idxs[2]).SetBondOrder(1)
        omol.OBMol.GetBond(idxs[0], idxs[1]).SetBondOrder(2)
        omol.OBMol.GetAtom(idxs[0]).SetFormalCharge(0)
        omol.OBMol.GetAtom(idxs[-1]).SetFormalCharge(-1)
    return fresh_omol_charge_radical(omol)


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
        moloplogger.debug(f"{DEBUG_TAG} | Cleaning resonance 8: {idxs}")
        omol.OBMol.GetBond(idxs[3], idxs[4]).SetBondOrder(1)
        omol.OBMol.GetBond(idxs[2], idxs[3]).SetBondOrder(2)
        omol.OBMol.GetBond(idxs[1], idxs[2]).SetBondOrder(1)
        omol.OBMol.GetBond(idxs[0], idxs[1]).SetBondOrder(2)
        omol.OBMol.GetAtom(idxs[0]).SetFormalCharge(0)
        omol.OBMol.GetAtom(idxs[-1]).SetFormalCharge(-1)
    return fresh_omol_charge_radical(omol)


def clean_resonances_9(omol: pybel.Molecule) -> pybel.Molecule:
    """
    Clean splited charges

    `[*+]-[*-]>>[*]=[*]`

    Parameters:
        omol (pybel.Molecule): The input molecule object

    Returns:
        pybel.Molecule: The cleaned molecule.
    """
    smarts = pybel.Smarts("[*+,*+2,*+3]-,=[*-,*-2,*-3]")
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
            moloplogger.debug(f"{DEBUG_TAG} | Cleaning splited charges {idxs}")
            omol.OBMol.GetBond(idxs[0], idxs[1]).SetBondOrder(
                int(omol.OBMol.GetBond(idxs[0], idxs[1]).GetBondOrder() + 1)
            )
            pos_atom.SetFormalCharge(omol.OBMol.GetAtom(idxs[0]).GetFormalCharge() - 1)
            neg_atom.SetFormalCharge(omol.OBMol.GetAtom(idxs[1]).GetFormalCharge() + 1)
    return fresh_omol_charge_radical(omol)


def clean_resonances_10(omol: pybel.Molecule) -> pybel.Molecule:
    """
    Clean double radical

    `[*]-[*]=,#[*]-[*]>>[*]=[*]-,=[*]=[*]`

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
        if atom1.GetSpinMultiplicity() and atom2.GetSpinMultiplicity():
            moloplogger.debug(f"{DEBUG_TAG} | Cleaning double radical")
            omol.OBMol.GetBond(idxs[1], idxs[2]).SetBondOrder(
                int(omol.OBMol.GetBond(idxs[1], idxs[2]).GetBondOrder() - 1)
            )
            omol.OBMol.GetBond(idxs[0], idxs[1]).SetBondOrder(
                int(omol.OBMol.GetBond(idxs[0], idxs[1]).GetBondOrder() + 1)
            )
            omol.OBMol.GetBond(idxs[2], idxs[3]).SetBondOrder(
                int(omol.OBMol.GetBond(idxs[2], idxs[3]).GetBondOrder() + 1)
            )
    return fresh_omol_charge_radical(omol)


def clean_resonances_11(omol: pybel.Molecule) -> pybel.Molecule:
    smarts = pybel.Smarts("[#7v3+0,#8v2+0,#16v2+0]-,=,~[*+1]")
    res = list(smarts.findall(omol))
    while len(res):
        idxs = res.pop(0)
        # atom1 = omol.OBMol.GetAtom(idxs[0])
        atom2 = omol.OBMol.GetAtom(idxs[1])
        if pt.GetDefaultValence(atom2.GetAtomicNum()) - atom2.GetTotalValence() >= 1:
            if (
                int(omol.OBMol.GetBond(idxs[0], idxs[1]).GetBondOrder()) == 1
                or int(omol.OBMol.GetBond(idxs[0], idxs[1]).GetBondOrder()) == 2
            ):
                moloplogger.debug(
                    f"{DEBUG_TAG} | Cleaning N, O, S with neighboring positive charge"
                )
                omol.OBMol.GetBond(idxs[0], idxs[1]).SetBondOrder(
                    int(omol.OBMol.GetBond(idxs[0], idxs[1]).GetBondOrder() + 1)
                )
                omol.OBMol.GetAtom(idxs[0]).SetFormalCharge(1)
                omol.OBMol.GetAtom(idxs[1]).SetFormalCharge(0)
    return fresh_omol_charge_radical(omol)


def clean_resonances_12(omol: pybel.Molecule) -> pybel.Molecule:
    smarts = pybel.Smarts("[#7v3+0,#8v2+0,#16v2+0]-,~[*]=,~[*]-,~[*+1]")
    res = list(smarts.findall(omol))
    while len(res):
        idxs = res.pop(0)
        # atom1 = omol.OBMol.GetAtom(idxs[0])
        atom4 = omol.OBMol.GetAtom(idxs[3])
        if pt.GetDefaultValence(atom4.GetAtomicNum()) - atom4.GetTotalValence() >= 1:
            if (
                int(omol.OBMol.GetBond(idxs[0], idxs[1]).GetBondOrder()) == 1
                and int(omol.OBMol.GetBond(idxs[1], idxs[2]).GetBondOrder()) == 2
                and int(omol.OBMol.GetBond(idxs[2], idxs[3]).GetBondOrder()) == 1
            ):
                moloplogger.debug(
                    f"{DEBUG_TAG} | Cleaning N, O, S with neighboring positive charge"
                )
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
    return fresh_omol_charge_radical(omol)


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
    return fresh_omol_charge_radical(omol)


def valid_metal_valence_radical(
    metal: str | int,
    valence,
    numradical,
) -> bool:
    """
    Check whether the given metal atom with the given valence and radical number is valid or not.

    Parameters:
        metal (str | int): The atomic number or symbol of the metal atom.
        valence (int): The valence of the metal atom.
        numradical (int): The radical number of the metal atom.

    Returns:
        bool: Whether the metal atom with the given valence and radical number is valid or not.
    """
    metal_symbol = pt.GetElementSymbol(metal) if isinstance(metal, int) else metal
    return (
        valence
        in metal_valence_avialable_prior[metal_symbol]
        + metal_valence_avialable_minor[metal_symbol]
    ) and (numradical in get_possible_metal_radicals(metal_symbol, valence))
