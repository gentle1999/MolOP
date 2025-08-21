"""
Including functions related to the three-dimensional structure of molecules
"""

import itertools
from copy import deepcopy
from typing import Iterable, List, Literal, Sequence, Tuple, Union, overload

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolTransforms
from rdkit.Chem.rdDistGeom import EmbedMolecule
from rdkit.Geometry import Point3D
from scipy.spatial.transform import Rotation as R

from ..utils.types import RdMol, RWMol
from .utils import pt

RDLogger.DisableLog("rdApp.*")  # type: ignore


def merge_mols(
    mols: List[RdMol | RWMol], scale=3.0, overlap_threshold=0.9, random_state=3407
) -> RdMol:
    """
    Merge a list of molecules.

    Parameters:
        mols (List[RdMol | RWMol]):
            List of RDkit molecules.
            All molecules must have conformers.
        scale (float):
            intermolecular scale
        overlap_threshold (float):
            maximum atomic spacing allowed between molecules
        random_state (int):
            random number generator seed

    Returns:
        Chem.rdchem.Mol: merged RDkit molecules
    """
    rng = np.random.RandomState(random_state)
    mols_copy = deepcopy(mols)
    points = fibonacci_sphere(
        len(mols_copy), rng.randint(0, min(random_state * 100, 65535))
    )
    for mol, point in zip(mols_copy, points, strict=True):
        translate_mol(mol, xyz_to_point(point) * scale)
    merged_mols = mols_copy[0]
    for mol in mols_copy[1:]:
        merged_mols = Chem.CombineMols(merged_mols, mol)
    while is_overlap(merged_mols, threshold=overlap_threshold):
        for mol, point in zip(mols_copy, points, strict=True):
            translate_mol(mol, xyz_to_point(point) * scale)
        merged_mols = mols_copy[0]
        for mol in mols_copy[1:]:
            merged_mols = Chem.CombineMols(merged_mols, mol)
    return merged_mols


def is_overlap(mol: RdMol | RWMol, threshold=0.9) -> bool:
    """
    Determine whether the molecular complexes overlap.
    Based on whether there is an atomic spacing between molecules that is less than the threshold value.

    Parameters:
        mol (RdMol | RWMol):
            molecule to be determined whether overlap or not
        threshold (float):
            maximum atomic spacing allowed between molecules

    Returns:
        bool: True if overlap, False otherwise
    """
    frags = Chem.GetMolFrags(mol)
    for frag1, frag2 in itertools.combinations(frags, 2):
        for atid1, atid2 in itertools.product(frag1, frag2):
            if Point3D.Distance(
                mol.GetConformer().GetAtomPosition(atid1),
                mol.GetConformer().GetAtomPosition(atid2),
            ) < threshold * (
                pt.GetRcovalent(mol.GetAtomWithIdx(atid1).GetAtomicNum())
                + pt.GetRcovalent(mol.GetAtomWithIdx(atid2).GetAtomicNum())
            ):
                return True
    return False


def xyz_to_point(xyz: Sequence | np.ndarray) -> Point3D:
    """
    Convert a sequence of three numbers representing an XYZ coordinate to a Point3D object.

    Parameters:
        xyz (Sequence): A sequence of three numbers representing an XYZ coordinate.

    Returns:
        rdkit.Geometry.Point3D: A Point3D object representing the given XYZ coordinate.
    """
    assert len(xyz) == 3, "xyz must have length 3"
    return Point3D(*map(float, xyz))


def set_conformer_position(mol: RdMol, positions: np.ndarray, conformer_id=0):
    """
    Set conformer position of a molecule.

    Parameters:
        mol (RdMol):
            molecule to be set conformer position
        positions (np.ndarray):
            conformer position
        conformer_id (int):
            conformer id
    """
    if mol.GetNumConformers() == 0:
        EmbedMolecule(mol)
    for atom in mol.GetAtoms():
        mol.GetConformer(conformer_id).SetAtomPosition(
            atom.GetIdx(), xyz_to_point(positions[atom.GetIdx()])
        )
    Chem.SanitizeMol(mol)


def get_random_vector(scale=5.0, random_state=3407) -> Point3D:
    """
    Generate a random vector with a given scale and random state.

    Parameters:
        scale (float): The scale factor for the vector. Default is 5.0.
        random_state (int): The random seed for reproducibility. Default is 3407.

    Returns:
        rdkit.Geometry.Point3D: A Point3D object representing the scaled random vector.
    """
    rng = np.random.RandomState(random_state)
    vec = rng.rand(3) - 0.5
    unit_vec = vec / np.linalg.norm(vec)
    return xyz_to_point(unit_vec) * scale


def fibonacci_sphere(samples=1, random_state=3407) -> List[Tuple[float, float, float]]:
    """
    Generate points on a Fibonacci sphere surface.

    Parameters:
        samples (int):
            number of sample points to generate
        random_state (int):
            random number generator seed

    Returns:
        List[Tuple[float, float, float]]:
            list of points in format (x, y, z)
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
        points.append((x, y, z))
    return points


def translate_mol(mol: RdMol, vector: Point3D, conformer_id=0):
    """
    Translates the coordinates of a molecule by a given vector.

    Parameters:
        mol (RdMol):
            The molecule to translate.
        vector (Point3D):
            The vector by which to translate the molecule.
        conformer_id (int):
            The conformer ID to translate. Defaults to 0.
    """
    for at in mol.GetAtoms():
        mol.GetConformer(conformer_id).SetAtomPosition(
            at.GetIdx(),
            mol.GetConformer(conformer_id).GetAtomPosition(at.GetIdx()) + vector,
        )


def translate_sub_mol(
    mol: RdMol, vector: Point3D, idx: Iterable[int], conformer_id=0
) -> RdMol:
    """
    Translates a subset of atoms in a molecule by a given vector.

    Parameters:
        mol (RdMol):
            The molecule to translate.
        vector (Point3D):
            The vector by which to translate the atoms.
        idx (Iterable[int]):
            The indices of the atoms to translate.
        conformer_id (int):
            The conformer ID to use for translation. Defaults to 0.

    Returns:
        Chem.rdchem.Mol: A copy of the molecule with the specified atoms translated.
    """
    mol_copy = deepcopy(mol)
    for at_idx in idx:
        mol_copy.GetConformer(conformer_id).SetAtomPosition(
            at_idx, mol_copy.GetConformer(conformer_id).GetAtomPosition(at_idx) + vector
        )
    return mol_copy


@overload
def translate_anchor(pos: np.ndarray, idx: int) -> np.ndarray: ...
@overload
def translate_anchor(pos: np.ndarray, idx: Sequence[int]) -> np.ndarray: ...
def translate_anchor(pos: np.ndarray, idx: Union[int, Sequence[int]]) -> np.ndarray:
    if isinstance(idx, (int, np.integer)):
        return pos - pos[idx]
    elif isinstance(idx, Sequence):
        for i in idx:
            assert isinstance(
                i, (int, np.integer)
            ), "idx must be int or sequence of int"
        return pos - np.mean(pos[idx], axis=0)
    else:
        raise TypeError("idx must be int or sequence of int")


def translate_mol_anchor(mol: RdMol, idx: Union[int, Sequence[int]], conformer_id=0):
    """
    Translate the entire molecule such that the atomic coordinates of the specified index are the origin

    Parameters:
        mol (RdMol):
            Molecule to be translated
        idx (int | Sequence[int]):
            The index of the anchor atom
        conformer_id (int):
            The conformer of id to be translated
    """
    set_conformer_position(
        mol,
        translate_anchor(mol.GetConformer(conformer_id).GetPositions(), idx),
        conformer_id,
    )


@overload
def rotate_anchor_to_axis(
    pos: np.ndarray, idx: int, axis: Literal["x", "y", "z"]
) -> np.ndarray: ...
@overload
def rotate_anchor_to_axis(
    pos: np.ndarray, idx: Sequence[int], axis: Literal["x", "y", "z"]
) -> np.ndarray: ...
def rotate_anchor_to_axis(
    pos: np.ndarray, idx: Union[int, Sequence[int]], axis: Literal["x", "y", "z"]
) -> np.ndarray:
    if isinstance(idx, (int, np.integer)):
        point = pos[idx]
    elif isinstance(idx, Sequence):
        for i in idx:
            assert isinstance(
                i, (int, np.integer)
            ), "idx must be int or sequence of int"
        point = np.mean(pos[idx], axis=0)
    else:
        raise TypeError("idx must be int or sequence of int")
    if axis == "x":
        r1 = R.from_euler("z", -np.arctan2(point[1], point[0]))
        positions = r1.apply(pos)
        r2 = R.from_euler("y", np.arctan2(positions[idx][2], positions[idx][0]))
        positions = r2.apply(positions)
    elif axis == "y":
        r1 = R.from_euler("x", -np.arctan2(point[2], point[1]))
        positions = r1.apply(pos)
        r2 = R.from_euler("z", np.arctan2(positions[idx][0], positions[idx][1]))
        positions = r2.apply(positions)
    elif axis == "z":
        r1 = R.from_euler("x", np.arctan2(point[1], point[2]))
        positions = r1.apply(pos)
        r2 = R.from_euler("y", -np.arctan2(positions[idx][0], positions[idx][2]))
        positions = r2.apply(positions)
    else:
        raise ValueError("axis must be one of 'x', 'y', or 'z'")
    return positions


def rotate_mol_anchor_to_axis(
    mol: RdMol,
    idx: Union[int, Sequence[int]],
    axis: Literal["x", "y", "z"] = "x",
    conformer_id=0,
):
    """
    Suppose the specified atom is A, and by rotating along two axes, point A is rotated to the positive half of the other axis.

    Parameters:
        mol (RdMol):
            Molecule to be rotated
        idx (int):
            The index of the anchor atom
        axis (Literal["x", "y", "z"]):
            The axis to rotate along
        conformer_id (int):
            The conformer of id to be rotated
    """
    set_conformer_position(
        mol,
        rotate_anchor_to_axis(mol.GetConformer(conformer_id).GetPositions(), idx, axis),
        conformer_id,
    )


@overload
def rotate_anchor_to_plane(
    pos: np.ndarray, idx: int, plane: Literal["xy", "yz", "zx", "yx", "zy", "xz"]
) -> np.ndarray: ...
@overload
def rotate_anchor_to_plane(
    pos: np.ndarray,
    idx: Sequence[int],
    plane: Literal["xy", "yz", "zx", "yx", "zy", "xz"],
) -> np.ndarray: ...
def rotate_anchor_to_plane(
    pos: np.ndarray,
    idx: Union[int, Sequence[int]],
    plane: Literal["xy", "yz", "zx", "yx", "zy", "xz"],
):
    if isinstance(idx, (int, np.integer)):
        point = pos[idx]
    elif isinstance(idx, Sequence):
        for i in idx:
            assert isinstance(
                i, (int, np.integer)
            ), "idx must be int or sequence of int"
        point = np.mean(pos[idx], axis=0)
    else:
        raise TypeError("idx must be int or sequence of int")
    if plane == "xy":
        r = R.from_euler("x", -np.arctan2(point[2], point[1]))
    elif plane == "yx":
        r = R.from_euler("y", -np.arctan2(point[2], point[0]))
    elif plane == "yz":
        r = R.from_euler("y", -np.arctan2(point[0], point[2]))
    elif plane == "zy":
        r = R.from_euler("z", np.arctan2(point[0], point[1]))
    elif plane == "zx":
        r = R.from_euler("z", -np.arctan2(point[1], point[0]))
    elif plane == "xz":
        r = R.from_euler("x", -np.arctan2(point[1], point[2]))
    else:
        raise ValueError("plane must be one of 'xy', 'yz', 'zx', 'yx', 'zy', or 'xz'")
    positions = r.apply(pos)
    return positions


def rotate_mol_anchor_to_plane(
    mol: RdMol,
    idx: int,
    plane: Literal["xy", "yz", "zx", "yx", "zy", "xz"] = "xy",
    conformer_id=0,
):
    """
    Suppose the specified atom is A and rotates along the X-axis, rotating point A to the XY plane.
    And the other two axes are rotated to the positive half of the X-axis and Y-axis, respectively.

    Parameters:
        mol (RdMol):
            Molecule to be rotated
        idx (int):
            The index of the anchor atom
        plane (Literal["xy", "yz", "zx", "yx", "zy", "xz"]):
            The plane to rotate along
        conformer_id (int):
            The conformer of id to be rotated
    """
    set_conformer_position(
        mol,
        rotate_anchor_to_plane(
            mol.GetConformer(conformer_id).GetPositions(), idx, plane
        ),
        conformer_id,
    )


def standard_orient(mol: RdMol, idx_list: Sequence[int], conformer_id=0):
    """
    Depending on the input `idx_list`, `translate_anchor`, `rotate_anchor_to_axis`, and `rotate_anchor_to_plane` are executed in order to obtain the normalized oriented molecule.\n
    If the length of the input `idx_list` is greater than 3, subsequent atomic numbers are not considered.

    Parameters:
        mol (RdMol):
            Molecule to be oriented
        idx_list (Sequence[int]):
            A Sequence of indices of the atoms to be translated, rotated, and rotated again
        conformer_id (int):
            The conformer of id to be oriented
    """
    methods = (
        translate_mol_anchor,
        rotate_mol_anchor_to_axis,
        rotate_mol_anchor_to_plane,
    )
    for method, idx in zip(methods, idx_list, strict=False):
        method(mol, idx, conformer_id=conformer_id)


def standard_orient_all_conformer(mol: RdMol, idx_list: Sequence[int]):
    """
    All the conformers of the molecule will be oriented.\n
    Depending on the input `idx_list`, `translate_anchor`, `rotate_anchor_to_axis`, and `rotate_anchor_to_plane` are executed in order to obtain the normalized oriented molecule.\n
    If the length of the input `idx_list` is greater than 3, subsequent atomic numbers are not considered.

    Parameters:
        mol (RdMol):
            Molecule to be oriented
        idx_list (Sequence[int]):
            A list of indices of the atoms to be translated, rotated, and rotated again
    """
    conformers = mol.GetConformers()
    if len(conformers) == 0:
        raise ValueError("Molecule must have at least one conformer")
    for conformer in conformers:
        standard_orient(mol, idx_list, conformer_id=conformer.GetId())


def unique_conformer(mol: RdMol, idx: int) -> RdMol:
    """
    Remove all conformers except for the specified conformer index.

    Parameters:
        mol (RdMol): The molecule object.
        idx (int): The index of the conformer to keep.

    Returns:
        Chem.rdchem.Mol: The molecule object with only the specified conformer.
    """
    temp_mol = deepcopy(mol)
    conformer_idx = [conf.GetId() for conf in temp_mol.GetConformers()]
    for conformer_id in conformer_idx:
        if conformer_id != idx:
            temp_mol.RemoveConformer(conformer_id)
    temp_mol.GetConformer(idx).SetId(0)
    return temp_mol


def get_geometry_info(mol: RdMol, atom_idxs: Sequence[int]) -> float:
    """
    Get the geometry infos among the atoms

    Parameters:
        mol (RdMol):
            Molecule to be oriented
        atom_idxs (Sequence[int]):
            A Sequence of index of the atoms, starts from 0

    Returns:
        float:
            - If the length of atom_idxs is 2, the bond length with unit Angstrom between the two atoms will be returned.
            - If the length of atom_idxs is 3, the angle with unit degree between  the three atoms will be returned.
            - If the length of atom_idxs is 4, the dihedral angle with unit degree between the four atoms will be returned.
    """
    if len(atom_idxs) == 2:
        return rdMolTransforms.GetBondLength(mol.GetConformer(), *atom_idxs)
    elif len(atom_idxs) == 3:
        return rdMolTransforms.GetAngleRad(mol.GetConformer(), *atom_idxs) / np.pi * 180
    elif len(atom_idxs) == 4:
        return (
            rdMolTransforms.GetDihedralRad(mol.GetConformer(), *atom_idxs) / np.pi * 180
        )
    else:
        raise ValueError("The length of atom_idxs must be 2, 3, or 4")
