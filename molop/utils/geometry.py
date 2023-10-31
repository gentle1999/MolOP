"""
Author: TMJ
Date: 2023-05-22 10:45:29
LastEditors: TMJ
LastEditTime: 2023-06-02 11:00:54
Description: Including functions related to the three-dimensional structure of molecules
"""

import itertools
from copy import deepcopy
from typing import Iterable, List, Sequence, Optional

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Geometry import Point3D
from scipy.spatial.transform import Rotation as R

from .types import RdMol, RdConformer
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")


def merge_mols(
    mols: List[RdMol], scale=3.0, overlap_threshold=1, random_state=3407
) -> RdMol:
    """
    Merge a list of molecules.

    Parameters:
        mols List[RdMol]:
            List of RDkit molecules.
            All molecules must have conformers.
        scale float:
            intermolecular scale
        overlap_threshold float:
            maximum atomic spacing allowed between molecules
        random_state int:
            random number generator seed

    Returns:
        RdMol:
            merged RDkit molecules
    """
    rng = np.random.RandomState(random_state)
    mols_copy = deepcopy(mols)
    points = fibonacci_sphere(
        len(mols_copy), rng.randint(0, min(random_state * 100, 65535))
    )
    for mol, point in zip(mols_copy, points):
        translate_mol(mol, xyz_to_point(point) * scale)
    merged_mols = mols_copy[0]
    for mol in mols_copy[1:]:
        merged_mols = Chem.CombineMols(merged_mols, mol)
    while is_overlap(merged_mols, threshold=overlap_threshold):
        for mol, point in zip(mols_copy, points):
            translate_mol(mol, xyz_to_point(point) * scale)
        merged_mols = mols_copy[0]
        for mol in mols_copy[1:]:
            merged_mols = Chem.CombineMols(merged_mols, mol)
    return merged_mols


def is_overlap(mol: RdMol, threshold=1.0) -> bool:
    """
    Determine whether the molecular complexes overlap.\n
    Based on whether there is an atomic spacing between molecules that is less than the threshold value.

    Parameters:
        mol RdMol:
            molecule to be determined whether overlap or not
        threshold float:
            maximum atomic spacing allowed between molecules

    Returns:
        Bool:
            True if overlap, False otherwise
    """
    frags = Chem.GetMolFrags(mol)
    for ix1, frag1 in enumerate(frags):
        for ix2, frag2 in enumerate(frags):
            if ix1 != ix2:
                for atid1 in frag1:
                    for atid2 in frag2:
                        if (
                            Point3D.Distance(
                                mol.GetConformer().GetAtomPosition(atid1),
                                mol.GetConformer().GetAtomPosition(atid2),
                            )
                            < threshold
                        ):
                            return True
    return False


def xyz_to_point(xyz: Sequence) -> Point3D:
    return Point3D(float(xyz[0]), float(xyz[1]), float(xyz[2]))


def set_conformer_position(mol: RdMol, positions: np.ndarray, conformer_id=0):
    for atom in mol.GetAtoms():
        mol.GetConformer(conformer_id).SetAtomPosition(
            atom.GetIdx(), xyz_to_point(positions[atom.GetIdx()])
        )


def get_random_vector(scale=5.0, random_state=3407):
    rng = np.random.RandomState(random_state)
    vec = rng.rand(3) - 0.5
    unit_vec = vec / np.linalg.norm(vec)
    return Point3D(unit_vec[0], unit_vec[1], unit_vec[2]) * scale


def fibonacci_sphere(samples=1, random_state=3407):
    """
    Generate points on a Fibonacci sphere surface.

    Parameters:
        samples int:
            number of sample points to generate
        random_state int:
            random number generator seed

    Reutrns:
        List[List[float]]:
            list of points in format [x, y, z]
    """
    rng = np.random.RandomState(random_state)
    rnd = rng.random() * samples

    points = []
    offset = 2.0 / samples
    increment = np.pi * (3.0 - np.sqrt(5.0))
    for i in range(samples):
        y = ((i * offset) - 1) + (offset / 2)
        r = np.sqrt(1 - pow(y, 2))
        phi = ((i + rnd) % samples) * increment
        x = np.cos(phi) * r
        z = np.sin(phi) * r
        points.append([x, y, z])
    return points


def translate_mol(mol: RdMol, vector: Point3D, conformer_id=0):
    for at in mol.GetAtoms():
        mol.GetConformer(conformer_id).SetAtomPosition(
            at.GetIdx(),
            mol.GetConformer(conformer_id).GetAtomPosition(at.GetIdx()) + vector,
        )


def translate_sub_mol(mol: RdMol, vector: Point3D, idx: Iterable[int], conformer_id=0):
    mol_copy = deepcopy(mol)
    for at_idx in idx:
        mol_copy.GetConformer(conformer_id).SetAtomPosition(
            at_idx, mol_copy.GetConformer(conformer_id).GetAtomPosition(at_idx) + vector
        )
    return mol_copy


def translate_anchor(mol: RdMol, idx: int, conformer_id=0):
    """
    Translate the entire molecule such that the atomic coordinates of the specified index are the origin

    Parameters:
        mol RdMol:
            Molecule to be translated
        idx int:
            The index of the anchor atom
        conformer_id int:
            The conformer of id to be translated
    """
    vector = mol.GetConformer(conformer_id).GetAtomPosition(idx) * -1.0
    translate_mol(mol, vector, conformer_id)


def rotate_anchor_to_X(mol: RdMol, idx: int, conformer_id=0):
    """
    Suppose the specified atom is A, and by rotating along the z and y axes, point A is rotated to the positive half of the X axis.

    Parameters:
        mol RdMol:
            Molecule to be rotated
        idx int:
            The index of the anchor atom
        conformer_id int:
            The conformer of id to be rotated
    """
    original_positions = mol.GetConformer(conformer_id).GetPositions()
    r1 = R.from_euler(
        "z", -np.arctan2(original_positions[idx][1], original_positions[idx][0])
    )
    positions = r1.apply(original_positions)
    r2 = R.from_euler("y", np.arctan2(positions[idx][2], positions[idx][0]))
    positions = r2.apply(positions)
    set_conformer_position(mol, positions, conformer_id)


def rotate_anchor_to_XY(mol: RdMol, idx: int, conformer_id=0):
    """
    Suppose the specified atom is A and rotates along the X-axis, rotating point A to quadrant 1 or 2 of the XY plane.

    Paraemters:
        mol RdMol:
            Molecule to be rotated
        idx int:
            The index of the anchor atom
        conformer_id int:
            The conformer of id to be rotated
    """
    original_positions = mol.GetConformer(conformer_id).GetPositions()
    r = R.from_euler(
        "x", -np.arctan2(original_positions[idx][2], original_positions[idx][1])
    )
    positions = r.apply(original_positions)
    set_conformer_position(mol, positions, conformer_id)


def standard_orient(mol: RdMol, idx_list: Sequence[int], conformer_id=0):
    """
    Depending on the input `idx_list`, `translate_anchor`, `rotate_anchor_to_X`, and `rotate_anchor_to_XY` are executed in order to obtain the normalized oriented molecule.\n
    If the length of the input `idx_list` is greater than 3, subsequent atomic numbers are not considered.

    Parameters:
        mol RdMol:
            Molecule to be oriented
        idx_list Sequence[int]:
            A list of indices of the atoms to be translated, rotated, and rotated again
        conformer_id int:
            The conformer of id to be oriented
    """
    methods = [translate_anchor, rotate_anchor_to_X, rotate_anchor_to_XY]
    for method, idx in zip(methods, idx_list):
        method(mol, idx, conformer_id)


def standard_orient_all_conformer(mol: RdMol, idx_list: Sequence[int]):
    """
    All the conformers of the molecule will be oriented.\n
    Depending on the input `idx_list`, `translate_anchor`, `rotate_anchor_to_X`, and `rotate_anchor_to_XY` are executed in order to obtain the normalized oriented molecule.\n
    If the length of the input `idx_list` is greater than 3, subsequent atomic numbers are not considered.

    Parameters:
        mol RdMol:
            Molecule to be oriented
        idx_list Sequence[int]:
            A list of indices of the atoms to be translated, rotated, and rotated again
    """
    conformers = mol.GetConformers()
    if len(conformers) == 0:
        return "No conformers found!"
    for conformer in conformers:
        standard_orient(mol, idx_list, conformer.GetId())


def unigue_conformer(mol: RdMol, idx: int):
    temp_mol = deepcopy(mol)
    conformer_idx = [conf.GetId() for conf in temp_mol.GetConformers()]
    for conformer_id in conformer_idx:
        if conformer_id != idx:
            temp_mol.RemoveConformer(conformer_id)
    temp_mol.GetConformer(idx).SetId(0)
    return temp_mol
