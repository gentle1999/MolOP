"""
Author: TMJ
Date: 2025-01-15 23:01:22
LastEditors: TMJ
LastEditTime: 2025-12-24 19:38:50
Description: 请填写简介
"""

import dataclasses
import itertools
from collections.abc import Sequence
from functools import lru_cache

import numpy as np
import timeout_decorator
from openbabel import openbabel as ob
from openbabel import pybel

from molop.config import molopconfig, moloplogger
from molop.utils.consts import (
    HETEROATOM,
    METAL_VALENCE_AVAILABLE_MINOR,
    METAL_VALENCE_AVAILABLE_PRIOR,
    NON_METAL_DICT,
    get_possible_metal_radicals,
)


DEBUG_TAG = "[GRAPH RECONSTRUCTION]"


@dataclasses.dataclass
class MetalAtomPosition:
    """
    The metal atom.
    """

    idx: int
    symbol: str
    element_idx: int
    valence: int
    radical_num: int
    position_x: float
    position_y: float
    position_z: float


def assign_radical_dots(atom: ob.OBAtom) -> int:
    return max(
        0,
        ob.GetTypicalValence(atom.GetAtomicNum(), atom.GetTotalValence(), atom.GetFormalCharge())
        - atom.GetTotalValence(),
    )


def assign_charge_radical_for_atom(atom: ob.OBAtom):
    """
    Assign the charge and radical dots for the overbonded atom.
    """
    if assign_radical_dots(atom):
        atom.SetSpinMultiplicity(assign_radical_dots(atom))
    else:
        if (
            NON_METAL_DICT[atom.GetAtomicNum()].num_outer_electrons == 3
            and atom.GetTotalValence() == 4
        ):
            atom.SetFormalCharge(-1)
        else:
            low_valence_total_elec = (
                NON_METAL_DICT[atom.GetAtomicNum()].num_outer_electrons
                + atom.GetTotalValence()
                + atom.GetSpinMultiplicity()
                - atom.GetFormalCharge()
            ) % 8
            high_valence_total_elec = (
                NON_METAL_DICT[atom.GetAtomicNum()].num_outer_electrons
                - atom.GetTotalValence()
                + atom.GetSpinMultiplicity()
                - atom.GetFormalCharge()
            ) % 2
            if low_valence_total_elec == 0:
                return
            if low_valence_total_elec <= high_valence_total_elec:
                atom.SetFormalCharge(low_valence_total_elec)
            else:
                atom.SetSpinMultiplicity(
                    NON_METAL_DICT[atom.GetAtomicNum()].num_outer_electrons
                    - atom.GetTotalValence()
                    + atom.GetSpinMultiplicity()
                    - atom.GetFormalCharge()
                )


def log_omol_infos(omol: pybel.Molecule, comment: str = ""):
    """
    Log the information of the omol.
    """
    moloplogger.debug(f"{DEBUG_TAG} | {comment}")
    moloplogger.debug(
        f"{DEBUG_TAG} | The omol contains {omol.OBMol.NumAtoms()} atoms "
        f"and {omol.OBMol.NumBonds()} bonds. smiles: {omol.write('smi').strip()}"  # type: ignore
    )
    for atom in omol.atoms:
        obatom: ob.OBAtom = atom.OBAtom
        moloplogger.debug(
            f"{DEBUG_TAG} | Atom {obatom.GetIdx():03d} is {ob.GetSymbol(obatom.GetAtomicNum()):<2s} "
            f"with formal charge {obatom.GetFormalCharge():+d} and "
            f"spin multiplicity {obatom.GetSpinMultiplicity():+d}."
        )


def fresh_omol_charge_radical(omol: pybel.Molecule) -> pybel.Molecule:
    """
    Fresh the omol, set the spin multiplicity of radical atoms to 1.
    """
    for atom in omol.atoms:
        assign_charge_radical_for_atom(atom.OBAtom)
    return omol


def combine_metal_with_omol(
    omol: pybel.Molecule, metal_list: Sequence[MetalAtomPosition]
) -> pybel.Molecule:
    """
    Combine the metal atom with the omol.
    """
    obmol: ob.OBMol = omol.clone.OBMol
    for metal in metal_list:
        atom: ob.OBAtom = obmol.NewAtom()
        atom.SetAtomicNum(metal.element_idx)
        atom.SetFormalCharge(metal.valence)
        atom.SetSpinMultiplicity(metal.radical_num)
        moloplogger.debug(
            f"{DEBUG_TAG} | Set atom {atom.GetIdx():03d} (original index: {metal.idx:03d}) position to "
            f"({metal.position_x:.6f}, {metal.position_y:.6f}, {metal.position_z:.6f})."
        )
        atom.SetVector(metal.position_x, metal.position_y, metal.position_z)
        obmol.RenumberAtoms(
            list(range(1, metal.idx)) + [atom.GetIdx()] + list(range(metal.idx + 1, atom.GetIdx())),
        )
    recovered_omol = pybel.Molecule(obmol)
    log_omol_infos(recovered_omol, "After adding metal atoms")
    return recovered_omol


@timeout_decorator.timeout(molopconfig.max_structure_recovery_time)
def xyz2omol(
    xyz_block: str, total_charge: int = 0, total_radical_electrons: int = 0
) -> pybel.Molecule | None:
    """
    Recover the structure of a molecule from its XYZ block.
    """
    omol = pybel.readstring("xyz", xyz_block)
    removable_metal_atoms: list[ob.OBAtom] = [
        atom.OBAtom for atom in omol.atoms if atom.OBAtom.IsMetal()
    ]
    available_valence_radical_states = [
        [
            MetalAtomPosition(
                idx=obatom.GetIdx(),
                symbol=ob.GetSymbol(obatom.GetAtomicNum()),
                element_idx=obatom.GetAtomicNum(),
                valence=valence,
                radical_num=radical_num,
                position_x=obatom.GetX(),
                position_y=obatom.GetY(),
                position_z=obatom.GetZ(),
            )
            for valence in METAL_VALENCE_AVAILABLE_PRIOR[ob.GetSymbol(obatom.GetAtomicNum())]
            + METAL_VALENCE_AVAILABLE_MINOR[ob.GetSymbol(obatom.GetAtomicNum())]
            for radical_num in get_possible_metal_radicals(
                ob.GetSymbol(obatom.GetAtomicNum()), valence
            )
        ]
        for obatom in removable_metal_atoms
    ]
    for obatom in removable_metal_atoms:
        omol.OBMol.DeleteAtom(obatom)
    log_omol_infos(omol, "After removing metal atoms")
    no_metal_xyz = omol.write("xyz")
    possible_metal_valence_radical_product = [
        product
        for product in itertools.product(*available_valence_radical_states)
        if sum(metal_pos.radical_num for metal_pos in product) <= total_radical_electrons
    ]
    possible_omols: list[pybel.Molecule] = []
    for metal_atom_product in possible_metal_valence_radical_product:
        total_metal_charge = sum(metal_pos.valence for metal_pos in metal_atom_product)
        total_metal_radical_electrons = sum(
            metal_pos.radical_num for metal_pos in metal_atom_product
        )
        moloplogger.debug(
            f"{DEBUG_TAG} | Metal atom product: with total charge {total_metal_charge} "
            f"and total radical {total_metal_radical_electrons}"
        )
        for metal_pos in metal_atom_product:
            moloplogger.debug(
                f"{DEBUG_TAG} | Metal atom {metal_pos.idx} with symbol {metal_pos.symbol}"
                f" and valence {metal_pos.valence} and radical {metal_pos.radical_num}"
            )
        possible_omol = xyz_to_omol_no_metal(
            no_metal_xyz,
            total_charge - total_metal_charge,
            total_radical_electrons - total_metal_radical_electrons,
        )
        if possible_omol is None:
            continue
        possible_omols.append(combine_metal_with_omol(possible_omol, metal_atom_product))
    if not possible_omols:
        return None

    scored_omols: list[tuple[float, pybel.Molecule]] = [
        (omol_score(res), res) for res in possible_omols
    ]
    scored_omols.sort(key=lambda x: x[0])
    scores_records = "\n".join(
        f"{score:.4f} {res.write('smi', opt={'k': True}).strip()}"  # type: ignore
        for score, res in scored_omols
    )
    moloplogger.debug(f"{DEBUG_TAG} | Scores with smiles: \n{scores_records}")
    possible_omols = [res for score, res in scored_omols]
    return possible_omols[0]


def validate_omol(
    omol: pybel.Molecule, total_charge: int = 0, total_radical_electrons: int = 0
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
    if (charge_sum := sum(atom.OBAtom.GetFormalCharge() for atom in omol.atoms)) != total_charge:
        moloplogger.debug(
            f"{DEBUG_TAG} | Charge check failed, total charge: "
            f"{total_charge}, actual charge: {charge_sum}"
        )
        return False

    radical_sum = sum(atom.OBAtom.GetSpinMultiplicity() for atom in omol.atoms)
    radical_sum_singlet = sum(atom.OBAtom.GetSpinMultiplicity() % 2 for atom in omol.atoms)
    if radical_sum_singlet == total_radical_electrons:
        radical_sum = radical_sum_singlet
    if radical_sum != total_radical_electrons:
        moloplogger.debug(
            f"{DEBUG_TAG} | Radical check failed, total radical: "
            f"{total_radical_electrons}, actual radical: {radical_sum}"
        )
        return False
    return True


@lru_cache(maxsize=1024)
def xyz_to_omol_no_metal(
    xyz_block: str, total_charge: int = 0, total_radical_electrons: int = 0
) -> pybel.Molecule | None:
    moloplogger.debug(
        f"{DEBUG_TAG} | Structure recovery started"
        f" with total charge {total_charge} and total radical {total_radical_electrons}"
    )
    if total_radical_electrons < 0:
        moloplogger.debug(f"{DEBUG_TAG} | Cannot recover the structure with negative radical")
        return None
    omol = pybel.readstring("xyz", xyz_block)
    omol = make_connections(omol)  # Connect the additional atoms if possible
    omol = pre_clean(omol)  # Clean the molecule before recovery
    omol = fresh_omol_charge_radical(omol)
    log_omol_infos(omol, "Input smiles: ")

    given_charge = total_charge - sum(atom.OBAtom.GetFormalCharge() for atom in omol.atoms)

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
    omol = fresh_omol_charge_radical(omol)
    log_omol_infos(omol, "After initial elimination: ")

    if validate_omol(omol, total_charge, total_radical_electrons):
        moloplogger.debug(
            f"{DEBUG_TAG} | Final charge: {sum(atom.OBAtom.GetFormalCharge() for atom in omol.atoms)} "
            f"total charge: {total_charge}; Final radical: {omol.OBMol.GetAtom(1).GetSpinMultiplicity()}"
            f" total radical: {total_radical_electrons}; Final smiles: \n{omol.write('smi').strip()}"  # type: ignore
        )
        return clean_resonances(omol)
    possible_resonances = get_radical_resonances(omol)
    moloplogger.debug(
        f"{DEBUG_TAG} | Possible resonance structures number: {len(possible_resonances)}"
    )
    recovered_resonances: list[pybel.Molecule] = []
    for idx, resonance in enumerate(possible_resonances):
        charge = given_charge
        log_omol_infos(
            resonance,
            f"Resonance index: {idx}, charge to be allocated: {charge}, "
            f"radical to be allocated: {total_radical_electrons}",
        )
        resonance, charge = process_resonance(resonance, charge)
        if validate_omol(resonance, total_charge, total_radical_electrons):
            recovered_resonances.append(resonance)
    if len(recovered_resonances) == 0:
        moloplogger.debug(
            f"{DEBUG_TAG} | Cannot recover the structure with given charge and radical"
        )
        return None

    scored_resonances: list[tuple[float, pybel.Molecule]] = [
        (omol_score(res), res) for res in recovered_resonances
    ]
    scored_resonances.sort(key=lambda x: x[0])
    scores_records = "\n".join(
        f"{score:.4f} {res.write('smi', opt={'k': True}).strip()}"  # type: ignore
        for score, res in scored_resonances
    )
    moloplogger.debug(f"{DEBUG_TAG} | Scores with smiles: \n{scores_records}")
    recovered_resonances = [res for score, res in scored_resonances]
    final_omol = recovered_resonances[0]

    moloplogger.debug(
        f"{DEBUG_TAG} | Final charge: {sum(atom.OBAtom.GetFormalCharge() for atom in final_omol.atoms)} "
        f"total charge: {total_charge}; Final radical: {final_omol.OBMol.GetAtom(1).GetSpinMultiplicity()}"
        f" total radical: {total_radical_electrons}; Final smiles: \n{final_omol.write('smi').strip()}"  # type: ignore
    )
    return final_omol


@lru_cache(maxsize=1024)
def process_resonance(resonance: pybel.Molecule, charge: int):
    resonance, charge = eliminate_1_3_dipole(resonance, charge)
    resonance, charge = eliminate_positive_charges(resonance, charge)
    resonance, charge = eliminate_negative_charges(resonance, charge)
    resonance = clean_neighbor_radicals(resonance)
    log_omol_infos(resonance, f"After resonance elimination, charge: {charge}")
    resonance = clean_resonances(resonance)
    return resonance, charge


def make_connections(omol: pybel.Molecule, factor: float = 1.4) -> pybel.Molecule:
    """
    Make connections between atoms in the pybel molecule.
    Parameters:
        omol (pybel.Molecule): The pybel molecule to be connected.
    Returns:
        pybel.Molecule: The connected pybel molecule.
    """
    donate_smarts = pybel.Smarts("[Nv0,Cv1,Nv3,Clv1,Clv2,Clv3,Brv1,Brv2,Brv3,Iv1,Iv2,Iv3]")
    accept_smarts = pybel.Smarts(
        "[Hv0,Bv2,Bv3,Cv0,Cv1,Cv2,Cv3,Nv1,Nv2,Ov0,Ov1,Clv0,Siv3,Pv2,Sv0,Sv1,Brv0,Iv0]"
    )
    donate_atoms: list[int] = list(itertools.chain(*donate_smarts.findall(omol)))
    accept_atoms: list[int] = list(itertools.chain(*accept_smarts.findall(omol)))

    for atom in donate_atoms:
        pairs = [(atom, accept_atom) for accept_atom in accept_atoms]
        pairs = sorted(
            pairs,
            key=lambda x: omol.OBMol.GetAtom(x[0]).GetDistance(omol.OBMol.GetAtom(x[1])),
        )
        for pair_1, pair_2 in pairs:
            donate_atom: ob.OBAtom = omol.OBMol.GetAtom(pair_1)
            accept_atom: ob.OBAtom = omol.OBMol.GetAtom(pair_2)
            distance = donate_atom.GetDistance(accept_atom)
            if (
                distance
                < (
                    ob.GetCovalentRad(donate_atom.GetAtomicNum())
                    + ob.GetCovalentRad(accept_atom.GetAtomicNum())
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
                    moloplogger.debug(f"{DEBUG_TAG} | Set bond order {pair_1} - {pair_2} to 1")
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
    smarts = pybel.Smarts("[Cv5,Nv5,Pv5,Siv5]=,#[*]")
    while res := smarts.findall(omol):
        idxs = res.pop(0)
        obbond: ob.OBBond = omol.OBMol.GetBond(idxs[0], idxs[1])
        obbond.SetBondOrder(obbond.GetBondOrder() - 1)

    smarts = pybel.Smarts("[#6]1([#6]2)([#6]3)[#7]23[#6]1")
    while res := smarts.findall(omol):
        idxs = res.pop(0)
        bcp_n = None
        bcp_c = None
        for idx in idxs:
            indexs = set(idxs) - {idx}
            if all(omol.OBMol.GetBond(idx, idx_2) for idx_2 in indexs):
                if omol.OBMol.GetAtom(idx).GetAtomicNum() == 7:
                    bcp_n = idx
                if omol.OBMol.GetAtom(idx).GetAtomicNum() == 6:
                    bcp_c = idx
        if bcp_n is not None and bcp_c is not None:
            omol.OBMol.DeleteBond(omol.OBMol.GetBond(bcp_n, bcp_c))
            moloplogger.debug(f"{DEBUG_TAG} | Fix N-BCP: {bcp_n} - {bcp_c}")

    smarts = pybel.Smarts("[#6]1([#6]2)[#7]2[#6]1")
    while res := smarts.findall(omol):
        idxs = res.pop(0)
        amine_n = None
        butyl_c = None
        for idx in idxs:
            indexs = set(idxs) - {idx}
            if all(omol.OBMol.GetBond(idx, idx_2) for idx_2 in indexs):
                if omol.OBMol.GetAtom(idx).GetAtomicNum() == 7:
                    amine_n = idx
                if omol.OBMol.GetAtom(idx).GetAtomicNum() == 6:
                    butyl_c = idx
        if amine_n is not None and butyl_c is not None:
            omol.OBMol.DeleteBond(omol.OBMol.GetBond(amine_n, butyl_c))
            moloplogger.debug(f"{DEBUG_TAG} | Fix N-BCP: {amine_n} - {butyl_c}")

    smarts = pybel.Smarts("[Siv5]-[O,F]")
    while res := smarts.findall(omol):
        idxs = res.pop(0)
        omol.OBMol.DeleteBond(omol.OBMol.GetBond(idxs[0], idxs[1]))
        moloplogger.debug(f"{DEBUG_TAG} | Remove Over bonding Si: {idxs[0]} - {idxs[1]}")

    return omol


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
            atom1.SetSpinMultiplicity(atom1.GetSpinMultiplicity() - 1)
            atom3.SetSpinMultiplicity(atom3.GetSpinMultiplicity() + 1)
    return omol


def eliminate_high_positive_charge_atoms(
    omol: pybel.Molecule, given_charge: int
) -> tuple[pybel.Molecule, int]:
    smarts = pybel.Smarts("[*+1,*+2,*+3]-[Ov1+0,Nv2+0,Sv1+0]")
    while res := smarts.findall(omol):
        idxs = res.pop(0)
        atom2: ob.OBAtom = omol.OBMol.GetAtom(idxs[1])
        if (
            -sum(atom.GetFormalCharge() for atom in ob.OBAtomAtomIter(omol.OBMol.GetAtom(idxs[0])))
            >= omol.OBMol.GetAtom(idxs[0]).GetFormalCharge()
        ):
            break
        atom2.SetSpinMultiplicity(atom2.GetSpinMultiplicity() - 1)
        atom2.SetFormalCharge(atom2.GetFormalCharge() - 1)
        moloplogger.debug(
            f"{DEBUG_TAG} | Eliminate atom {omol.OBMol.GetAtom(idxs[1]).GetAtomicNum()} {idxs[1]} with -1 charge"
        )
        given_charge += 1
    return omol, given_charge


def eliminate_CN_in_doubt(omol: pybel.Molecule, given_charge: int) -> tuple[pybel.Molecule, int]:
    smarts = pybel.Smarts("[#6v4+0]=,#[#7v4+1,#15v4+1]")
    doubt_pair = smarts.findall(omol)
    CN_in_doubt = len(doubt_pair)
    if CN_in_doubt > 0:
        moloplogger.debug(
            f"{DEBUG_TAG} | Fixing CN in doubt, number of pairs: {CN_in_doubt}, "
            f"charge to be allocated: {given_charge}"
        )
    if CN_in_doubt % 2 == 0 and CN_in_doubt > 0:
        for atom_1_idx, atom_2_idx in doubt_pair[: CN_in_doubt // 2]:
            atom_1: ob.OBAtom = omol.OBMol.GetAtom(atom_1_idx)
            atom_2: ob.OBAtom = omol.OBMol.GetAtom(atom_2_idx)
            bond: ob.OBBond = omol.OBMol.GetBond(atom_1_idx, atom_2_idx)
            atom_1.SetFormalCharge(-1)
            bond.SetBondOrder(bond.GetBondOrder() - 1)
            atom_2.SetFormalCharge(0)
            given_charge += 2
            moloplogger.debug(
                f"{DEBUG_TAG} | Fix CN in doubt: {atom_1_idx} - {atom_2}, "
                f"charge to be allocated: {given_charge}"
            )
    return omol, given_charge


def eliminate_carboxyl(omol: pybel.Molecule, given_charge: int) -> tuple[pybel.Molecule, int]:
    smarts = pybel.Smarts("[Ov1+0]-C=O")
    res = list(smarts.findall(omol))
    while len(res):
        idxs = res.pop(0)
        atom_1: ob.OBAtom = omol.OBMol.GetAtom(idxs[0])
        atom_1.SetSpinMultiplicity(atom_1.GetSpinMultiplicity() - 1)
        atom_1.SetFormalCharge(atom_1.GetFormalCharge() - 1)
        given_charge += 1
        moloplogger.debug(f"{DEBUG_TAG} | Eliminate atom {idxs[0]} with -1 charge")
    return omol, given_charge


def eliminate_carbene_neighbor_heteroatom(
    omol: pybel.Molecule, given_charge: int
) -> tuple[pybel.Molecule, int]:
    """
    Fix the carbine neighbor heteroatom.
    Parameters:
        omol (pybel.Molecule): The molecule to be transformed.
        given_charge (int): The given charge.
    Returns:
        int: The remaining charge to be allocated.
    """
    for atom in omol.atoms:
        obatom: ob.OBAtom = atom.OBAtom
        if obatom.GetSpinMultiplicity() == 2:
            for neighbor in ob.OBAtomAtomIter(obatom):
                if neighbor.GetSpinMultiplicity():
                    return omol, given_charge
            for neighbor in ob.OBAtomAtomIter(obatom):
                if (
                    neighbor.GetAtomicNum() in HETEROATOM
                    and neighbor.GetFormalCharge() == 0
                    and neighbor.GetSpinMultiplicity() == 0
                ):
                    bond: ob.OBBond = obatom.GetBond(neighbor)
                    bond.SetBondOrder(bond.GetBondOrder() + 1)
                    obatom.SetSpinMultiplicity(0)
                    obatom.SetFormalCharge(obatom.GetFormalCharge() - 1)
                    neighbor.SetFormalCharge(neighbor.GetFormalCharge() + 1)
                    moloplogger.debug(
                        f"{DEBUG_TAG} | eliminating carbine neighbor: {obatom.GetIdx()} and "
                        f"{neighbor.GetIdx()}, charge to be allocated: {given_charge}"
                    )
                    break
    return omol, given_charge


def eliminate_NNN(omol: pybel.Molecule, given_charge: int) -> tuple[pybel.Molecule, int]:
    smarts = pybel.Smarts("[#7v1+0]-[#7v2+0]-[#7v1+0]")
    while res := smarts.findall(omol):
        idxs = res.pop(0)
        atom1: ob.OBAtom = omol.OBMol.GetAtom(idxs[0])
        atom2: ob.OBAtom = omol.OBMol.GetAtom(idxs[1])
        atom3: ob.OBAtom = omol.OBMol.GetAtom(idxs[2])
        bond1: ob.OBBond = omol.OBMol.GetBond(idxs[0], idxs[1])
        bond1.SetBondOrder(bond1.GetBondOrder() + 1)
        bond2: ob.OBBond = omol.OBMol.GetBond(idxs[1], idxs[2])
        bond2.SetBondOrder(bond2.GetBondOrder() + 1)
        atom1.SetSpinMultiplicity(atom1.GetSpinMultiplicity() - 1)
        atom1.SetFormalCharge(atom1.GetFormalCharge() - 1)
        atom2.SetSpinMultiplicity(atom2.GetSpinMultiplicity() - 1)
        atom2.SetFormalCharge(atom2.GetFormalCharge() + 1)
        atom3.SetSpinMultiplicity(atom3.GetSpinMultiplicity() - 1)
        atom3.SetFormalCharge(atom3.GetFormalCharge() - 1)
        given_charge += 1
        moloplogger.debug(f"{DEBUG_TAG} | Fix NNN-: {idxs[0]}, {idxs[1]}, {idxs[2]}")
    return omol, given_charge


def break_deformed_ene(
    omol: pybel.Molecule, given_charge: int = 0, given_radical: int = 0, tolerance=5
) -> pybel.Molecule:
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
            return omol
        idxs = res.pop(0)
        bond2: ob.OBBond = omol.OBMol.GetBond(idxs[0], idxs[1])
        if bond2.IsRotor() or bond2.GetBondOrder() == 1:
            continue
        torsion_angle = abs(omol.OBMol.GetTorsion(*idxs))
        torsion_angle = min(torsion_angle, 180 - torsion_angle)
        if torsion_angle > tolerance:
            moloplogger.debug(
                f"{DEBUG_TAG} | torsion angle: {torsion_angle}, tolerance: {tolerance}"
            )
            moloplogger.debug(f"{DEBUG_TAG} | break bond {idxs[0]} and {idxs[1]}")
            atom1_2: ob.OBAtom = omol.OBMol.GetAtom(idxs[0])
            atom2_2: ob.OBAtom = omol.OBMol.GetAtom(idxs[1])
            bond2.SetBondOrder(bond2.GetBondOrder() - 1)
            atom1_2.SetSpinMultiplicity(atom1_2.GetSpinMultiplicity() + 1)
            atom2_2.SetSpinMultiplicity(atom2_2.GetSpinMultiplicity() + 1)

    smarts = pybel.Smarts("[*]~[*+0](=,:[*+0])~[*]")
    res = list(smarts.findall(omol))
    while len(res):
        if (
            sum(atom.GetSpinMultiplicity() for atom in ob.OBMolAtomIter(omol.OBMol))
            >= abs(given_charge) + given_radical
        ):
            return omol
        idxs = res.pop(0)
        bond: ob.OBBond = omol.OBMol.GetBond(idxs[0], idxs[1])
        if bond.IsRotor() or bond.GetBondOrder() == 1:
            continue
        torsion_angle = abs(omol.OBMol.GetTorsion(*idxs))
        torsion_angle = min(torsion_angle, 180 - torsion_angle)
        if torsion_angle > tolerance:
            moloplogger.debug(
                f"{DEBUG_TAG} | torsion angle: {torsion_angle}, tolerance: {tolerance}"
            )
            moloplogger.debug(f"{DEBUG_TAG} | break bond {idxs[0]} and {idxs[1]}")
            atom1: ob.OBAtom = omol.OBMol.GetAtom(idxs[0])
            atom2: ob.OBAtom = omol.OBMol.GetAtom(idxs[1])
            bond.SetBondOrder(bond.GetBondOrder() - 1)
            atom1.SetSpinMultiplicity(atom1.GetSpinMultiplicity() + 1)
            atom2.SetSpinMultiplicity(atom2.GetSpinMultiplicity() + 1)
    return omol


def break_one_bond(
    omol: pybel.Molecule, given_charge: int = 0, given_radical: int = 0
) -> tuple[pybel.Molecule, int]:
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
    while res := smarts.findall(omol):
        if (
            sum(atom.GetSpinMultiplicity() for atom in ob.OBMolAtomIter(omol.OBMol))
            >= abs(given_charge) + given_radical
        ):
            return omol, given_charge
        # Get the indices of the first suitable bond's atoms
        idxs = res.pop(0)
        # Log debug information indicating the bond to be broken
        moloplogger.debug(
            f"{DEBUG_TAG} | break bond {idxs[0]} and {idxs[1]}: Double bond or triple bond"
        )
        # Get and reduce the bond order of the found bond
        bond: ob.OBBond = omol.OBMol.GetBond(idxs[0], idxs[1])
        bond.SetBondOrder(bond.GetBondOrder() - 1)
        begin_atom: ob.OBAtom = bond.GetBeginAtom()
        end_atom: ob.OBAtom = bond.GetEndAtom()
        begin_atom.SetSpinMultiplicity(begin_atom.GetSpinMultiplicity() + 1)
        end_atom.SetSpinMultiplicity(end_atom.GetSpinMultiplicity() + 1)

    smarts = pybel.Smarts("[#7+1,#15+1]=[*+0]")
    while res := smarts.findall(omol):
        if (
            sum(atom.GetSpinMultiplicity() for atom in ob.OBMolAtomIter(omol.OBMol))
            >= abs(given_charge) + given_radical
        ):
            return omol, given_charge
        # Get the indices of the first suitable bond's atoms
        idxs = res.pop(0)
        # Log debug information indicating the bond to be broken
        moloplogger.debug(f"{DEBUG_TAG} | break bond {idxs[0]} and {idxs[1]}: *[N+](*)=[O]")
        # Get and reduce the bond order of the found bond
        bond2: ob.OBBond = omol.OBMol.GetBond(idxs[0], idxs[1])
        bond2.SetBondOrder(bond2.GetBondOrder() - 1)
        begin_atom2: ob.OBAtom = bond2.GetBeginAtom()
        end_atom2: ob.OBAtom = bond2.GetEndAtom()
        end_atom2.SetSpinMultiplicity(end_atom2.GetSpinMultiplicity() + 1)
        begin_atom2.SetFormalCharge(int(begin_atom2.GetFormalCharge() - 1))
        given_charge += 1
    smarts = pybel.Smarts("[*+0]:[*+0]")
    # Loop to find suitable bonds, if only aromatic bonds are present
    while res := smarts.findall(omol):
        if (
            sum(atom.GetSpinMultiplicity() for atom in ob.OBMolAtomIter(omol.OBMol))
            >= abs(given_charge) + given_radical
        ):
            return omol, given_charge
        # Get the indices of the first suitable bond's atoms
        idxs = res.pop(0)
        # Log debug information indicating the bond to be broken
        moloplogger.debug(f"{DEBUG_TAG} | break bond {idxs[0]} and {idxs[1]}: Aromatic")
        # Get and reduce the bond order of the found bond
        bond3: ob.OBBond = omol.OBMol.GetBond(idxs[0], idxs[1])
        bond3.SetBondOrder(bond3.GetBondOrder() - 1)
        begin_atom3: ob.OBAtom = bond3.GetBeginAtom()
        end_atom3: ob.OBAtom = bond3.GetEndAtom()
        begin_atom3.SetSpinMultiplicity(begin_atom3.GetSpinMultiplicity() + 1)
        end_atom3.SetSpinMultiplicity(end_atom3.GetSpinMultiplicity() + 1)

    if all(bond.GetBondOrder() == 1 for bond in ob.OBMolBondIter(omol.OBMol)):
        for single_bond in ob.OBMolBondIter(omol.OBMol):
            if (
                sum(atom.GetSpinMultiplicity() for atom in ob.OBMolAtomIter(omol.OBMol))
                >= abs(given_charge) + given_radical
            ):
                return omol, given_charge
            moloplogger.debug(
                f"{DEBUG_TAG} | break bond {single_bond.GetBeginAtom().GetIdx()} and "
                f"{single_bond.GetEndAtom().GetIdx()}:"
                f" Single bond"
            )
            single_begin_atom: ob.OBAtom = single_bond.GetBeginAtom()
            single_end_atom: ob.OBAtom = single_bond.GetEndAtom()
            single_begin_atom.SetSpinMultiplicity(single_begin_atom.GetSpinMultiplicity() + 1)
            single_end_atom.SetSpinMultiplicity(single_end_atom.GetSpinMultiplicity() + 1)
            omol.OBMol.DeleteBond(single_bond)
    return omol, given_charge


def clean_neighbor_radicals(omol: pybel.Molecule) -> pybel.Molecule:
    """
    Clean the neighbor radicals.
    Parameters:
        omol (pybel.Molecule): The molecule to be transformed.
    """
    for bond in ob.OBMolBondIter(omol.OBMol):
        begin_atom: ob.OBAtom = bond.GetBeginAtom()
        end_atom: ob.OBAtom = bond.GetEndAtom()
        if begin_atom.GetSpinMultiplicity() and end_atom.GetSpinMultiplicity():
            bond_to_add = min(
                begin_atom.GetSpinMultiplicity(),
                end_atom.GetSpinMultiplicity(),
            )
            bond.SetBondOrder(bond.GetBondOrder() + bond_to_add)
            begin_atom.SetSpinMultiplicity(begin_atom.GetSpinMultiplicity() - bond_to_add)
            end_atom.SetSpinMultiplicity(end_atom.GetSpinMultiplicity() - bond_to_add)
            moloplogger.debug(
                f"{DEBUG_TAG} | Cleaning bond {begin_atom.GetIdx()} and {end_atom.GetIdx()}"
            )
    return omol


def eliminate_charge_spliting(
    omol: pybel.Molecule, given_charge: int
) -> tuple[pybel.Molecule, int]:
    if (
        all(atom.GetFormalCharge() == 0 for atom in ob.OBMolAtomIter(omol.OBMol))
        and sum(atom.GetSpinMultiplicity() for atom in ob.OBMolAtomIter(omol.OBMol)) >= 2
    ):
        radical_atoms: list[ob.OBAtom] = [
            atom for atom in ob.OBMolAtomIter(omol.OBMol) if atom.GetSpinMultiplicity()
        ]
        moloplogger.debug(
            f"{DEBUG_TAG} | Eliminate charge splitting, number of radical atoms: {len(radical_atoms)}"
        )
        while len(radical_atoms) > abs(given_charge) + 1:
            for atom in radical_atoms:
                if atom.GetAtomicNum() == 8:
                    atom.SetSpinMultiplicity(atom.GetSpinMultiplicity() - 1)
                    atom.SetFormalCharge(atom.GetFormalCharge() - 1)
                    given_charge += 1
                    radical_atoms.remove(atom)
                    break
            else:
                break
        while len(radical_atoms) > abs(given_charge) + 1:
            for atom in radical_atoms:
                if atom.GetAtomicNum() == 7:
                    atom.SetSpinMultiplicity(atom.GetSpinMultiplicity() - 1)
                    atom.SetFormalCharge(atom.GetFormalCharge() - 1)
                    given_charge += 1
                    radical_atoms.remove(atom)
                    break
            else:
                break
        while len(radical_atoms) > abs(given_charge) + 1:
            for atom in radical_atoms:
                if atom.GetAtomicNum() == 6 and not any(
                    _atom for _atom in ob.OBAtomAtomIter(atom) if _atom.GetAtomicNum() in HETEROATOM
                ):
                    atom.SetSpinMultiplicity(atom.GetSpinMultiplicity() - 1)
                    atom.SetFormalCharge(atom.GetFormalCharge() - 1)
                    given_charge += 1
                    radical_atoms.remove(atom)
                    break
            else:
                break
        while len(radical_atoms) > abs(given_charge) + 1:
            for atom in radical_atoms:
                if atom.GetAtomicNum() == 6:
                    atom.SetSpinMultiplicity(atom.GetSpinMultiplicity() - 1)
                    atom.SetFormalCharge(atom.GetFormalCharge() - 1)
                    given_charge += 1
                    radical_atoms.remove(atom)
                    break
            else:
                break
    return omol, given_charge


def get_one_step_resonance(omol: pybel.Molecule) -> list[pybel.Molecule]:
    """
    Get one step resonance structures from a molecule.

    Parameters:
        omol (pybel.Molecule): The molecule to be processed.

    Returns:
        List[pybel.Molecule]: A list of one step resonance structures.
    """
    # Get a list of atoms with unpaired electrons (radicals)
    smarts = pybel.Smarts("[*]-,=,:[*]=,#,:[*]")
    res = list(smarts.findall(omol))
    result = []
    for idxs in res:
        atom1: ob.OBAtom = omol.OBMol.GetAtom(idxs[0])
        atom3: ob.OBAtom = omol.OBMol.GetAtom(idxs[2])
        bond1: ob.OBBond = omol.OBMol.GetBond(idxs[0], idxs[1])
        bond2: ob.OBBond = omol.OBMol.GetBond(idxs[1], idxs[2])
        if (
            atom1.GetSpinMultiplicity() == 1
            and bond1.GetBondOrder() <= 2
            and bond2.GetBondOrder() >= 2
        ):
            new_omol = omol.clone
            bond1.SetBondOrder(bond1.GetBondOrder() + 1)
            bond2.SetBondOrder(bond2.GetBondOrder() - 1)
            atom1.SetSpinMultiplicity(atom1.GetSpinMultiplicity() - 1)
            atom3.SetSpinMultiplicity(atom3.GetSpinMultiplicity() + 1)
            result.append(new_omol)
    return result


def get_radical_resonances(omol: pybel.Molecule) -> list[pybel.Molecule]:
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


def eliminate_1_3_dipole(omol: pybel.Molecule, given_charge: int) -> tuple[pybel.Molecule, int]:
    smarts = pybel.Smarts("[*-1]-,=[N+0,O+0]-,=[*]")
    res = list(smarts.findall(omol))
    while len(res):
        idxs = res.pop(0)
        atom1: ob.OBAtom = omol.OBMol.GetAtom(idxs[0])
        atom2: ob.OBAtom = omol.OBMol.GetAtom(idxs[1])
        atom3: ob.OBAtom = omol.OBMol.GetAtom(idxs[2])
        bond2: ob.OBBond = omol.OBMol.GetBond(idxs[1], idxs[2])
        if (
            atom3.GetSpinMultiplicity()
            and NON_METAL_DICT[atom2.GetAtomicNum()].num_outer_electrons + atom2.GetTotalValence()
            == 8
        ):
            moloplogger.debug(
                f"{DEBUG_TAG} | Eliminate 1,3 dipole: {atom1.GetIdx()} "
                f"{atom2.GetIdx()} {atom3.GetIdx()}"
            )
            atom2.SetFormalCharge(atom2.GetFormalCharge() + 1)
            bond2.SetBondOrder(int(bond2.GetBondOrder() + 1))
            atom3.SetSpinMultiplicity(atom3.GetSpinMultiplicity() - 1)
            given_charge -= 1
    return omol, given_charge


def eliminate_positive_charges(
    omol: pybel.Molecule, given_charge: int
) -> tuple[pybel.Molecule, int]:
    smarts = pybel.Smarts("[Nv3+0]=[Nv2+0]")
    while given_charge > 0 and (res := smarts.findall(omol)):
        idxs = res.pop(0)
        abatom: ob.OBAtom = omol.OBMol.GetAtom(idxs[1])
        abatom.SetSpinMultiplicity(abatom.GetSpinMultiplicity() - 1)
        abatom.SetFormalCharge(1)
        given_charge -= 1
    smarts = pybel.Smarts("[#6v3+0,#6v2+0,#1v0+0]")
    while given_charge > 0 and (res := smarts.findall(omol)):
        idxs = res.pop(0)
        abatom2: ob.OBAtom = omol.OBMol.GetAtom(idxs[0])
        abatom2.SetSpinMultiplicity(abatom2.GetSpinMultiplicity() - 1)
        abatom2.SetFormalCharge(1)
        given_charge -= 1
    for atom in omol.atoms:
        obatom: ob.OBAtom = atom.OBAtom
        if given_charge <= 0:
            break
        if obatom.GetSpinMultiplicity() >= 1 and obatom.GetFormalCharge() == 0:
            to_add = min(obatom.GetSpinMultiplicity(), given_charge)
            obatom.SetSpinMultiplicity(obatom.GetSpinMultiplicity() - to_add)
            obatom.SetFormalCharge(to_add)
            given_charge -= to_add
    return omol, given_charge


def eliminate_negative_charges(
    omol: pybel.Molecule, given_charge: int
) -> tuple[pybel.Molecule, int]:
    possible_heteroatoms: list[tuple[ob.OBAtom, int]] = []
    for atom in omol.atoms:
        obatom: ob.OBAtom = atom.OBAtom
        if (
            obatom.GetAtomicNum() in HETEROATOM
            and obatom.GetFormalCharge() == 0
            and obatom.GetSpinMultiplicity() >= 1
        ):
            possible_heteroatoms.append((obatom, HETEROATOM.index(obatom.GetAtomicNum())))
    possible_heteroatoms.sort(key=lambda x: x[1])
    for obatom, _ in possible_heteroatoms:
        if given_charge >= 0:
            break
        to_add = min(obatom.GetSpinMultiplicity(), abs(given_charge))
        obatom.SetSpinMultiplicity(obatom.GetSpinMultiplicity() - to_add)
        obatom.SetFormalCharge(-to_add)
        given_charge += to_add
        moloplogger.debug(
            f"{DEBUG_TAG} | Eliminate negative charge: {obatom.GetIdx()} with charge {to_add}"
        )

    # Step 2.3.2: Try to find the heteroatom with negative charge first.
    smarts = pybel.Smarts("[#6v3+0]")
    while given_charge < 0 and (res := smarts.findall(omol)):
        idxs = res.pop(0)
        obatom1: ob.OBAtom = omol.OBMol.GetAtom(idxs[0])
        to_add = min(obatom1.GetSpinMultiplicity(), abs(given_charge))
        obatom1.SetSpinMultiplicity(obatom1.GetSpinMultiplicity() - to_add)
        obatom1.SetFormalCharge(-to_add)
        given_charge += to_add
        moloplogger.debug(
            f"{DEBUG_TAG} | Eliminate negative charge: {obatom1.GetIdx()} with charge {to_add}"
        )
    smarts = pybel.Smarts("[#1v0+0]")
    while given_charge < 0 and (res := smarts.findall(omol)):
        idxs = res.pop(0)
        obatom2: ob.OBAtom = omol.OBMol.GetAtom(idxs[0])
        to_add = min(obatom2.GetSpinMultiplicity(), abs(given_charge))
        obatom2.SetSpinMultiplicity(obatom2.GetSpinMultiplicity() - to_add)
        obatom2.SetFormalCharge(-to_add)
        given_charge += to_add
        moloplogger.debug(
            f"{DEBUG_TAG} | Eliminate negative charge: {obatom2.GetIdx()} with charge {to_add}"
        )
    smarts = pybel.Smarts("[#6v2+0,#6v1+0,#6v0+0]")
    while given_charge < 0 and (res := smarts.findall(omol)):
        idxs = res.pop(0)
        obatom3: ob.OBAtom = omol.OBMol.GetAtom(idxs[0])
        to_add = min(obatom3.GetSpinMultiplicity(), abs(given_charge))
        obatom3.SetSpinMultiplicity(obatom3.GetSpinMultiplicity() - to_add)
        obatom3.SetFormalCharge(-to_add)
        given_charge += to_add
        moloplogger.debug(
            f"{DEBUG_TAG} | Eliminate negative charge: {obatom3.GetIdx()} with charge {to_add}"
        )
    return omol, given_charge


def clean_resonances_0(omol: pybel.Molecule) -> pybel.Molecule:
    smarts = pybel.Smarts("[*-]-[*]=[*]~[*+]")
    res = list(smarts.findall(omol))
    while len(res):
        idxs = res.pop(0)
        obatom1: ob.OBAtom = omol.OBMol.GetAtom(idxs[0])
        obatom4: ob.OBAtom = omol.OBMol.GetAtom(idxs[-1])
        obbond1: ob.OBBond = omol.OBMol.GetBond(idxs[0], idxs[1])
        obbond2: ob.OBBond = omol.OBMol.GetBond(idxs[1], idxs[2])
        obbond3: ob.OBBond = omol.OBMol.GetBond(idxs[2], idxs[3])
        if (
            NON_METAL_DICT[obatom1.GetAtomicNum()].default_valence > obatom1.GetTotalValence()
            and NON_METAL_DICT[obatom4.GetAtomicNum()].default_valence > obatom4.GetTotalValence()
        ):
            moloplogger.debug(f"{DEBUG_TAG} | Cleaning resonance 0: {idxs}")
            obbond1.SetBondOrder(obbond1.GetBondOrder() + 1)
            obbond2.SetBondOrder(obbond2.GetBondOrder() - 1)
            obbond3.SetBondOrder(obbond3.GetBondOrder() + 1)
            obatom1.SetFormalCharge(obatom1.GetFormalCharge() + 1)
            obatom4.SetFormalCharge(obatom4.GetFormalCharge() - 1)
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
    res = smarts.findall(omol)
    while len(res):
        idxs = res.pop(0)
        obatom1: ob.OBAtom = omol.OBMol.GetAtom(idxs[0])
        obatom3: ob.OBAtom = omol.OBMol.GetAtom(idxs[-1])
        obbond1: ob.OBBond = omol.OBMol.GetBond(idxs[0], idxs[1])
        obbond2: ob.OBBond = omol.OBMol.GetBond(idxs[1], idxs[2])
        moloplogger.debug(f"{DEBUG_TAG} | Cleaning resonance 1: {idxs}")
        obbond1.SetBondOrder(obbond1.GetBondOrder() + 1)
        obbond2.SetBondOrder(obbond2.GetBondOrder() - 1)
        obatom1.SetFormalCharge(obatom1.GetFormalCharge() + 1)
        obatom3.SetFormalCharge(obatom3.GetFormalCharge() - 1)
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
    res = smarts.findall(omol)
    while len(res):
        idxs = res.pop(0)
        obatom1: ob.OBAtom = omol.OBMol.GetAtom(idxs[0])
        obatom5: ob.OBAtom = omol.OBMol.GetAtom(idxs[-1])
        obbond1: ob.OBBond = omol.OBMol.GetBond(idxs[0], idxs[1])
        obbond2: ob.OBBond = omol.OBMol.GetBond(idxs[1], idxs[3])
        obbond3: ob.OBBond = omol.OBMol.GetBond(idxs[3], idxs[4])
        obbond4: ob.OBBond = omol.OBMol.GetBond(idxs[4], idxs[5])
        if (
            obbond4.GetBondOrder() == 1
            and obbond3.GetBondOrder() == 2
            and obbond2.GetBondOrder() == 1
            and obbond1.GetBondOrder() == 2
            and obatom1.GetFormalCharge() == 0
            and obatom5.GetFormalCharge() == -1
        ):
            moloplogger.debug(f"{DEBUG_TAG} | Cleaning resonance 2: {idxs}")
            obbond4.SetBondOrder(2)
            obbond3.SetBondOrder(1)
            obbond2.SetBondOrder(2)
            obbond1.SetBondOrder(1)
            obatom1.SetFormalCharge(-1)
            obatom5.SetFormalCharge(0)
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
    res = smarts.findall(omol)
    while len(res):
        idxs = res.pop(0)
        obatom1: ob.OBAtom = omol.OBMol.GetAtom(idxs[0])
        obatom5: ob.OBAtom = omol.OBMol.GetAtom(idxs[-1])
        obbond1: ob.OBBond = omol.OBMol.GetBond(idxs[0], idxs[1])
        obbond2: ob.OBBond = omol.OBMol.GetBond(idxs[1], idxs[2])
        obbond3: ob.OBBond = omol.OBMol.GetBond(idxs[2], idxs[3])
        obbond4: ob.OBBond = omol.OBMol.GetBond(idxs[3], idxs[4])
        moloplogger.debug(f"{DEBUG_TAG} | Cleaning resonance 3: {idxs}")
        obbond4.SetBondOrder(2)
        obbond3.SetBondOrder(1)
        obbond2.SetBondOrder(2)
        obbond1.SetBondOrder(1)
        obatom1.SetFormalCharge(0)
        obatom5.SetFormalCharge(0)
    return omol


def clean_resonances_4(omol: pybel.Molecule) -> pybel.Molecule:
    """
    `[#7+,#8+]=[*]-[#6-,#7-,#8-]>>[#7,#8]-[*]=[#6,#7,#8]`

    Parameters:
        omol (pybel.Molecule): The input molecule object

    Returns:
        pybel.Molecule: The cleaned molecule.
    """
    smarts = pybel.Smarts("[#7+,#8+]=[*]-[#6-,#7-,#8-]")
    res = smarts.findall(omol)
    while len(res):
        idxs = res.pop(0)
        obatom1: ob.OBAtom = omol.OBMol.GetAtom(idxs[0])
        obatom3: ob.OBAtom = omol.OBMol.GetAtom(idxs[-1])
        obbond1: ob.OBBond = omol.OBMol.GetBond(idxs[0], idxs[1])
        obbond2: ob.OBBond = omol.OBMol.GetBond(idxs[1], idxs[2])
        moloplogger.debug(f"{DEBUG_TAG} | Cleaning resonance 4: {idxs}")
        obbond2.SetBondOrder(obbond2.GetBondOrder() + 1)
        obbond1.SetBondOrder(obbond1.GetBondOrder() - 1)
        obatom1.SetFormalCharge(0)
        obatom3.SetFormalCharge(0)
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
    res = smarts.findall(omol)
    while len(res):
        idxs = res.pop(0)
        obatom1: ob.OBAtom = omol.OBMol.GetAtom(idxs[0])
        obatom3: ob.OBAtom = omol.OBMol.GetAtom(idxs[-1])
        obbond1: ob.OBBond = omol.OBMol.GetBond(idxs[0], idxs[1])
        obbond2: ob.OBBond = omol.OBMol.GetBond(idxs[1], idxs[2])
        if (
            obbond1.GetBondOrder() == 2
            and obbond2.GetBondOrder() == 1
            and obatom3.GetFormalCharge() == -1
            and obatom1.GetFormalCharge() == 0
        ):
            moloplogger.debug(f"{DEBUG_TAG} | Cleaning resonance 5: {idxs}")
            obbond2.SetBondOrder(obbond2.GetBondOrder() + 1)
            obbond1.SetBondOrder(obbond1.GetBondOrder() - 1)
            obatom1.SetFormalCharge(-1)
            obatom3.SetFormalCharge(0)
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
    res = smarts.findall(omol)
    while len(res):
        idxs = res.pop(0)
        obatom1: ob.OBAtom = omol.OBMol.GetAtom(idxs[0])
        obatom3: ob.OBAtom = omol.OBMol.GetAtom(idxs[-1])
        obbond1: ob.OBBond = omol.OBMol.GetBond(idxs[0], idxs[1])
        obbond2: ob.OBBond = omol.OBMol.GetBond(idxs[1], idxs[2])
        moloplogger.debug(f"{DEBUG_TAG} | Cleaning resonance 6: {idxs}")
        obbond2.SetBondOrder(3)
        obbond1.SetBondOrder(1)
        obatom1.SetFormalCharge(-1)
        obatom3.SetFormalCharge(0)
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
    res = smarts.findall(omol)
    while len(res):
        idxs = res.pop(0)
        obatom1: ob.OBAtom = omol.OBMol.GetAtom(idxs[0])
        obatom3: ob.OBAtom = omol.OBMol.GetAtom(idxs[2])
        obbond1: ob.OBBond = omol.OBMol.GetBond(idxs[0], idxs[1])
        obbond2: ob.OBBond = omol.OBMol.GetBond(idxs[1], idxs[2])
        moloplogger.debug(f"{DEBUG_TAG} | Cleaning resonance 7: {idxs}")
        obbond2.SetBondOrder(1)
        obbond1.SetBondOrder(2)
        obatom1.SetFormalCharge(0)
        obatom3.SetFormalCharge(-1)
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
    res = smarts.findall(omol)
    while len(res):
        idxs = res.pop(0)
        moloplogger.debug(f"{DEBUG_TAG} | Cleaning resonance 8: {idxs}")
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
    smarts = pybel.Smarts("[*+,*+2,*+3]-,=[*-,*-2,*-3]")
    res = smarts.findall(omol)
    while len(res):
        idxs = res.pop(0)
        obatom1: ob.OBAtom = omol.OBMol.GetAtom(idxs[0])
        obatom2: ob.OBAtom = omol.OBMol.GetAtom(idxs[1])
        obbond1: ob.OBBond = omol.OBMol.GetBond(idxs[0], idxs[1])
        if (
            NON_METAL_DICT[obatom1.GetAtomicNum()].default_valence - obatom1.GetTotalValence() >= 1
            and NON_METAL_DICT[obatom2.GetAtomicNum()].default_valence - obatom2.GetTotalValence()
            >= 1
        ):
            moloplogger.debug(f"{DEBUG_TAG} | Cleaning splited charges {idxs}")
            bond_to_add = min(
                NON_METAL_DICT[obatom1.GetAtomicNum()].default_valence - obatom1.GetTotalValence(),
                NON_METAL_DICT[obatom2.GetAtomicNum()].default_valence - obatom2.GetTotalValence(),
            )
            obbond1.SetBondOrder(obbond1.GetBondOrder() + bond_to_add)
            obatom1.SetFormalCharge(obatom1.GetFormalCharge() - bond_to_add)
            obatom2.SetFormalCharge(obatom2.GetFormalCharge() + bond_to_add)
    return omol


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
    res = smarts.findall(omol)
    while len(res):
        idxs = res.pop(0)
        obatom1: ob.OBAtom = omol.OBMol.GetAtom(idxs[0])
        obatom3: ob.OBAtom = omol.OBMol.GetAtom(idxs[-1])
        obbond1: ob.OBBond = omol.OBMol.GetBond(idxs[0], idxs[1])
        obbond2: ob.OBBond = omol.OBMol.GetBond(idxs[1], idxs[2])
        obbond3: ob.OBBond = omol.OBMol.GetBond(idxs[2], idxs[3])
        if obatom1.GetSpinMultiplicity() == 1 and obatom3.GetSpinMultiplicity() == 1:
            moloplogger.debug(f"{DEBUG_TAG} | Cleaning double radical")
            obbond2.SetBondOrder(obbond2.GetBondOrder() - 1)
            obbond1.SetBondOrder(obbond1.GetBondOrder() + 1)
            obbond3.SetBondOrder(obbond3.GetBondOrder() + 1)
            obatom1.SetSpinMultiplicity(0)
            obatom3.SetSpinMultiplicity(0)
    return omol


def clean_resonances_11(omol: pybel.Molecule) -> pybel.Molecule:
    smarts = pybel.Smarts("[#7v3+0,#8v2+0,#16v2+0]-,=,:[*+1]")
    res = smarts.findall(omol)
    while len(res):
        idxs = res.pop(0)
        obatom1: ob.OBAtom = omol.OBMol.GetAtom(idxs[0])
        obatom2: ob.OBAtom = omol.OBMol.GetAtom(idxs[1])
        obbond1: ob.OBBond = omol.OBMol.GetBond(idxs[0], idxs[1])
        if NON_METAL_DICT[
            obatom2.GetAtomicNum()
        ].default_valence - obatom2.GetTotalValence() >= 1 and (
            obbond1.GetBondOrder() == 1 or obbond1.GetBondOrder() == 2
        ):
            moloplogger.debug(f"{DEBUG_TAG} | Cleaning N, O, S with neighboring positive charge")
            obbond1.SetBondOrder(obbond1.GetBondOrder() + 1)
            obatom1.SetFormalCharge(1)
            obatom2.SetFormalCharge(0)
    return omol


def clean_resonances_12(omol: pybel.Molecule) -> pybel.Molecule:
    smarts = pybel.Smarts("[#7v3+0,#8v2+0,#16v2+0]-,:[*]=,:[*]-,:[*+1]")
    res = smarts.findall(omol)
    while len(res):
        idxs = res.pop(0)
        obatom1: ob.OBAtom = omol.OBMol.GetAtom(idxs[0])
        obatom4: ob.OBAtom = omol.OBMol.GetAtom(idxs[3])
        obbond1: ob.OBBond = omol.OBMol.GetBond(idxs[0], idxs[1])
        obbond2: ob.OBBond = omol.OBMol.GetBond(idxs[1], idxs[2])
        obbond3: ob.OBBond = omol.OBMol.GetBond(idxs[2], idxs[3])

        if NON_METAL_DICT[
            obatom4.GetAtomicNum()
        ].default_valence - obatom4.GetTotalValence() >= 1 and (
            obbond1.GetBondOrder() == 1
            and obbond2.GetBondOrder() == 2
            and obbond3.GetBondOrder() == 1
        ):
            moloplogger.debug(f"{DEBUG_TAG} | Cleaning N, O, S with neighboring positive charge")
            obbond1.SetBondOrder(obbond1.GetBondOrder() + 1)
            obbond2.SetBondOrder(obbond2.GetBondOrder() - 1)
            obbond3.SetBondOrder(obbond3.GetBondOrder() + 1)
            obatom1.SetFormalCharge(1)
            obatom4.SetFormalCharge(0)
    return omol


def clean_resonances_13(omol: pybel.Molecule) -> pybel.Molecule:
    smarts = pybel.Smarts("[*-]:[*]=[#7+0,#8+0]")
    res = smarts.findall(omol)
    while len(res):
        idxs = res.pop(0)
        obatom1: ob.OBAtom = omol.OBMol.GetAtom(idxs[0])
        obatom3: ob.OBAtom = omol.OBMol.GetAtom(idxs[2])
        obbond1: ob.OBBond = omol.OBMol.GetBond(idxs[0], idxs[1])
        obbond2: ob.OBBond = omol.OBMol.GetBond(idxs[1], idxs[2])
        if (
            NON_METAL_DICT[obatom1.GetAtomicNum()].default_valence - obatom1.GetTotalValence() >= 1
            and obbond1.GetBondOrder() == 1
        ):
            moloplogger.debug(
                f"{DEBUG_TAG} | Cleaning aromatic ring with neighboring negative charge"
            )
            obbond1.SetBondOrder(obbond1.GetBondOrder() + 1)
            obbond2.SetBondOrder(obbond2.GetBondOrder() - 1)
            obatom1.SetFormalCharge(0)
            obatom3.SetFormalCharge(-1)
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
    omol = clean_resonances_13(omol)
    return omol


# ---- Score Functions ----


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


def get_deviation_score(omol: pybel.Molecule, atom_idx: int) -> float:
    def _get_coords(a: ob.OBAtom) -> list[float]:
        vec = a.GetVector()
        return [vec.GetX(), vec.GetY(), vec.GetZ()]

    atom: ob.OBAtom = omol.OBMol.GetAtom(atom_idx)
    neighbor_atoms = list(ob.OBAtomAtomIter(atom))

    if len(neighbor_atoms) == 2:
        angle = omol.OBMol.GetAngle(neighbor_atoms[0], atom, neighbor_atoms[1])
        return abs(angle - 108) / 108.0

    if len(neighbor_atoms) == 3:
        p1 = _get_coords(neighbor_atoms[0])
        p2 = _get_coords(neighbor_atoms[1])
        p3 = _get_coords(neighbor_atoms[2])
        p_atom = _get_coords(atom)
        quality = calculate_shape_quality(p1, p2, p3, p_atom)
        return 1.0 - quality

    return 0.0


def calc_symmetry_penalty(omol: pybel.Molecule) -> float:
    obmol: ob.OBMol = omol.OBMol
    gs = ob.OBGraphSym(obmol)
    symmetry_ids_vec = ob.vectorUnsignedInt()
    gs.GetSymmetry(symmetry_ids_vec)
    return (len(set(symmetry_ids_vec)) - obmol.NumAtoms()) * 2.0


def calculate_charge_penalty(atom: ob.OBAtom) -> float:
    penalty = 0.0
    charge = atom.GetFormalCharge()
    if charge == 0:
        return penalty
    total_electrons = (
        NON_METAL_DICT[atom.GetAtomicNum()].num_outer_electrons
        + atom.GetTotalValence()
        - atom.GetFormalCharge()
        + atom.GetSpinMultiplicity()
    )
    if total_electrons == 8 or total_electrons == 2:
        return penalty
    en = ob.GetElectroNeg(atom.GetAtomicNum())
    if charge > 0:
        penalty += abs(charge) * max(0.0, en - 2) * 3.0
    if charge < 0:
        penalty += abs(charge) * max(0.0, 4 - en) * 3.0
    return penalty


def calculate_radical_penalty(atom: ob.OBAtom) -> float:
    penalty = 0.0
    radical_num = atom.GetSpinMultiplicity()
    if radical_num == 0:
        return penalty
    if atom.GetAtomicNum() in HETEROATOM:
        return penalty + radical_num * 2.0
    return penalty + (3 - atom.GetHvyDegree()) * 1.5


def calculate_coulombic_penalty(bond: ob.OBBond) -> float:
    atom1: ob.OBAtom = bond.GetBeginAtom()
    atom2: ob.OBAtom = bond.GetEndAtom()
    q1 = atom1.GetFormalCharge()
    q2 = atom2.GetFormalCharge()

    if q1 == 0 or q2 == 0:
        return 0.0
    if q1 * q2 > 0:
        return 15.0
    return -0.5


def calculate_physchem_penalty(omol: pybel.Molecule) -> float:
    total_penalty = 0.0
    obmol = omol.OBMol

    for atom in ob.OBMolAtomIter(obmol):
        if atom.IsMetal():
            continue
        total_penalty += calculate_charge_penalty(atom)
        total_penalty += calculate_radical_penalty(atom)

    for bond in ob.OBMolBondIter(obmol):
        total_penalty += calculate_coulombic_penalty(bond)

    return total_penalty


def get_metal_coordination_sphere(
    omol: pybel.Molecule, metal_atom: ob.OBAtom, cutoff: float = 2.8
) -> list[tuple[ob.OBAtom, float]]:
    neighbors = []
    metal_idx = metal_atom.GetIdx()
    for atom in omol.atoms:
        obatom: ob.OBAtom = atom.OBAtom
        if obatom.GetIdx() == metal_idx:
            continue
        dist = metal_atom.GetDistance(obatom)
        if dist <= cutoff:
            neighbors.append((obatom, dist))
    return neighbors


def calculate_metal_penalty(omol: pybel.Molecule) -> float:
    penalty = 0.0

    for atom in omol.atoms:
        metal_atom: ob.OBAtom = atom.OBAtom
        if not metal_atom.IsMetal():
            continue

        symbol = ob.GetSymbol(metal_atom.GetAtomicNum())
        valence = metal_atom.GetFormalCharge()
        if valence <= 0:
            penalty += 10
        prior_list = METAL_VALENCE_AVAILABLE_PRIOR.get(symbol, [])
        minor_list = METAL_VALENCE_AVAILABLE_MINOR.get(symbol, [])

        if valence not in prior_list:
            if valence in minor_list:
                penalty += 2.0
            else:
                penalty += 10.0
        neighbors = get_metal_coordination_sphere(omol, metal_atom, cutoff=2.6)

        for ligand_atom, dist in neighbors:
            ligand_charge = ligand_atom.GetFormalCharge()
            if valence > 0:
                if ligand_charge > 0:
                    penalty += 10.0 * (ligand_charge * valence) / (dist**2)
                elif ligand_charge < 0:
                    penalty -= 2.0 * (abs(ligand_charge) * valence) / dist

    return penalty


def omol_score(omol: pybel.Molecule) -> float:
    """
    Calculate the resonance score of a molecule, the lower the score the better the resonance.

    This function can only be used for comparison between isomers recovered from the same
    set of atomic coordinates.

    Parameters:
        omol (pybel.Molecule): A pybel.Molecule object.

    Returns:
        float: The resonance score of a molecule.
    """
    obmol: ob.OBMol = omol.OBMol
    score: float = 0.0
    # --- 1. Symmetry Score ---
    score += calc_symmetry_penalty(omol)

    # --- 2. Geometry Deviation Penalties ---
    for atom_idx in range(1, obmol.NumAtoms() + 1):
        atom: ob.OBAtom = obmol.GetAtom(atom_idx)
        if atom.IsMetal():
            continue
        if atom.GetSpinMultiplicity() > 0:
            score += get_deviation_score(omol, atom_idx) * 10.0
        if atom.GetFormalCharge() > 0:
            score += get_deviation_score(omol, atom_idx) * 10.0
        if atom.GetFormalCharge() < 0:
            score += (1 - get_deviation_score(omol, atom_idx)) * 10.0

    # --- 3. Physical Chemistry Penalties ---
    score += calculate_physchem_penalty(omol)

    # --- 4. Metal Specific Penalties ---
    score += calculate_metal_penalty(omol)
    return score
