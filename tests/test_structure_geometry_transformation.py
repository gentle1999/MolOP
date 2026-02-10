import sys
from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem

from molop.structure.GeometryTransformation import (
    fibonacci_sphere,
    get_geometry_info,
    get_random_vector,
    is_overlap,
    merge_mols,
    rotate_anchor_to_axis,
    rotate_anchor_to_plane,
    rotate_mol_anchor_to_axis,
    rotate_mol_anchor_to_plane,
    standard_orient,
    standard_orient_all_conformer,
    translate_anchor,
    translate_mol,
    translate_mol_anchor,
    translate_sub_mol,
    unique_conformer,
    xyz_to_point,
)


sys.path.append(str(Path(__file__).resolve().parent))

from _helpers_structure import build_rdmol_with_conformer


def test_xyz_to_point_accepts_valid_xyz() -> None:
    point = xyz_to_point([1, 2.5, -3.0])
    assert (point.x, point.y, point.z) == pytest.approx((1.0, 2.5, -3.0), abs=1e-12)


@pytest.mark.parametrize(
    "bad_xyz",
    [
        [1.0, 2.0],
        [1.0, 2.0, 3.0, 4.0],
        np.array([[0.0, 1.0], [2.0, 3.0]]),
    ],
)
def test_xyz_to_point_rejects_invalid_shape(bad_xyz: list[float] | np.ndarray) -> None:
    with pytest.raises(AssertionError, match="xyz must have length 3"):
        xyz_to_point(bad_xyz)


def test_translate_anchor_handles_scalar_and_sequence_indices() -> None:
    pos = np.array(
        [
            [1.0, 2.0, 3.0],
            [2.0, 2.0, 3.0],
            [1.0, 4.0, 3.0],
        ]
    )

    translated_scalar = translate_anchor(pos, 0)
    translated_sequence_single = translate_anchor(pos, [0])
    translated_sequence_mean = translate_anchor(pos, [0, 1])

    assert translated_scalar == pytest.approx(translated_sequence_single, abs=1e-12)
    assert translated_scalar[0] == pytest.approx((0.0, 0.0, 0.0), abs=1e-12)
    assert translated_sequence_mean[0] == pytest.approx((-0.5, 0.0, 0.0), abs=1e-12)
    assert translated_sequence_mean[1] == pytest.approx((0.5, 0.0, 0.0), abs=1e-12)


def test_translate_anchor_rejects_invalid_index_types() -> None:
    pos = np.array([[1.0, 2.0, 3.0], [2.0, 2.0, 3.0]])

    with pytest.raises(TypeError, match="idx must be int or sequence of int"):
        translate_anchor(pos, 1.5)  # type: ignore[arg-type]

    with pytest.raises(AssertionError, match="idx must be int or sequence of int"):
        translate_anchor(pos, [0, "bad"])  # type: ignore[list-item]


@pytest.mark.parametrize("axis", ["x", "y", "z"])
def test_rotate_anchor_to_axis_aligns_anchor_to_requested_axis(axis: str) -> None:
    positions = np.array([[0.0, 0.0, 0.0], [1.0, 2.0, 3.0]])
    rotated = rotate_anchor_to_axis(positions, 1, axis)  # type: ignore[arg-type]

    anchor = rotated[1]
    if axis == "x":
        assert anchor[1] == pytest.approx(0.0, abs=1e-10)
        assert anchor[2] == pytest.approx(0.0, abs=1e-10)
    elif axis == "y":
        assert anchor[0] == pytest.approx(0.0, abs=1e-10)
        assert anchor[2] == pytest.approx(0.0, abs=1e-10)
    else:
        assert anchor[0] == pytest.approx(0.0, abs=1e-10)
        assert anchor[1] == pytest.approx(0.0, abs=1e-10)


@pytest.mark.parametrize("plane", ["xy", "yx", "yz", "zy", "zx", "xz"])
def test_rotate_anchor_to_plane_rotates_deterministically(plane: str) -> None:
    positions = np.array([[0.0, 0.0, 0.0], [1.0, 2.0, 3.0], [2.0, 0.5, 1.5]])
    original_mean = np.mean(positions[[1, 2]], axis=0)
    rotated = rotate_anchor_to_plane(positions, [1, 2], plane)  # type: ignore[arg-type]

    rotated_mean = np.mean(rotated[[1, 2]], axis=0)
    assert np.linalg.norm(rotated_mean) == pytest.approx(np.linalg.norm(original_mean), abs=1e-10)


def test_rotate_anchor_to_axis_rejects_invalid_axis() -> None:
    positions = np.array([[0.0, 0.0, 0.0], [1.0, 2.0, 3.0]])
    with pytest.raises(ValueError, match="axis must be one of"):
        rotate_anchor_to_axis(positions, 1, "invalid")  # type: ignore[arg-type]


def test_rotate_anchor_to_plane_rejects_invalid_plane() -> None:
    positions = np.array([[0.0, 0.0, 0.0], [1.0, 2.0, 3.0]])
    with pytest.raises(ValueError, match="plane must be one of"):
        rotate_anchor_to_plane(positions, 1, "invalid")  # type: ignore[arg-type]


def test_standard_orient_all_conformer_raises_without_conformer() -> None:
    mol = Chem.MolFromSmiles("O")
    assert mol is not None
    assert mol.GetNumConformers() == 0

    with pytest.raises(ValueError, match="at least one conformer"):
        standard_orient_all_conformer(mol, [0])


def test_vector_and_sphere_generators_are_deterministic() -> None:
    vec1 = get_random_vector(scale=2.5, random_state=77)
    vec2 = get_random_vector(scale=2.5, random_state=77)
    assert (vec1.x, vec1.y, vec1.z) == pytest.approx((vec2.x, vec2.y, vec2.z), abs=1e-12)
    assert np.linalg.norm([vec1.x, vec1.y, vec1.z]) == pytest.approx(2.5, abs=1e-10)

    pts1 = fibonacci_sphere(samples=4, random_state=99)
    pts2 = fibonacci_sphere(samples=4, random_state=99)
    assert pts1 == pytest.approx(pts2, abs=1e-12)
    assert len(pts1) == 4
    for point in pts1:
        assert np.linalg.norm(point) == pytest.approx(1.0, abs=1e-10)


def test_molecule_translation_helpers_preserve_expected_geometry() -> None:
    mol = build_rdmol_with_conformer(
        atom_symbols=["O", "H", "H"],
        bonds=[(0, 1, 1), (0, 2, 1)],
        coordinates=[
            (0.0, 0.0, 0.0),
            (1.0, 0.0, 0.0),
            (0.0, 1.0, 0.0),
        ],
    )
    conf = mol.GetConformer()
    before = conf.GetAtomPosition(0)

    translate_mol(mol, xyz_to_point([1.0, -2.0, 3.5]))
    moved = mol.GetConformer().GetAtomPosition(0)
    assert (moved.x - before.x, moved.y - before.y, moved.z - before.z) == pytest.approx(
        (1.0, -2.0, 3.5), abs=1e-12
    )

    sub_mol = translate_sub_mol(mol, xyz_to_point([10.0, 0.0, 0.0]), [1])
    sub_original = mol.GetConformer().GetAtomPosition(1)
    sub_shifted = sub_mol.GetConformer().GetAtomPosition(1)
    assert sub_shifted.x - sub_original.x == pytest.approx(10.0, abs=1e-12)
    atom0_original = mol.GetConformer().GetAtomPosition(0)
    atom0_shifted = sub_mol.GetConformer().GetAtomPosition(0)
    assert (atom0_shifted.x, atom0_shifted.y, atom0_shifted.z) == pytest.approx(
        (atom0_original.x, atom0_original.y, atom0_original.z),
        abs=1e-12,
    )


def test_orientation_helpers_mutate_molecule_without_embedding() -> None:
    mol = build_rdmol_with_conformer(
        atom_symbols=["C", "O", "H"],
        bonds=[(0, 1, 1), (0, 2, 1)],
        coordinates=[
            (1.0, 0.5, 0.2),
            (2.0, 1.5, 1.2),
            (0.8, 1.7, -0.1),
        ],
    )

    translate_mol_anchor(mol, 0)
    anchor_after_translate = mol.GetConformer().GetAtomPosition(0)
    assert (
        anchor_after_translate.x,
        anchor_after_translate.y,
        anchor_after_translate.z,
    ) == pytest.approx((0.0, 0.0, 0.0), abs=1e-10)

    rotate_mol_anchor_to_axis(mol, 1, axis="x")
    rotated_axis_anchor = mol.GetConformer().GetAtomPosition(1)
    assert rotated_axis_anchor.y == pytest.approx(0.0, abs=1e-8)
    assert rotated_axis_anchor.z == pytest.approx(0.0, abs=1e-8)

    rotate_mol_anchor_to_plane(mol, 2, plane="xy")
    rotated_plane_anchor = mol.GetConformer().GetAtomPosition(2)
    assert rotated_plane_anchor.z == pytest.approx(0.0, abs=1e-8)

    standard_orient(mol, [0, 1, 2])
    assert mol.GetNumAtoms() == 3


def test_standard_orient_all_conformer_and_unique_conformer() -> None:
    mol = build_rdmol_with_conformer(
        atom_symbols=["C", "O", "H"],
        bonds=[(0, 1, 1), (0, 2, 1)],
        coordinates=[
            (0.0, 0.0, 0.0),
            (1.0, 0.0, 0.0),
            (0.0, 1.0, 0.0),
        ],
    )
    second_conf = Chem.Conformer(mol.GetConformer())
    second_conf.SetAtomPosition(0, xyz_to_point([0.2, 0.2, 0.2]))
    second_conf.SetAtomPosition(1, xyz_to_point([1.2, 0.2, 0.2]))
    second_conf.SetAtomPosition(2, xyz_to_point([0.2, 1.2, 0.2]))
    mol.AddConformer(second_conf, assignId=True)
    assert mol.GetNumConformers() == 2

    standard_orient_all_conformer(mol, [0, 1, 2])
    for conf in mol.GetConformers():
        origin = conf.GetAtomPosition(0)
        assert (origin.x, origin.y, origin.z) == pytest.approx((0.0, 0.0, 0.0), abs=1e-8)

    only_one = unique_conformer(mol, 1)
    assert only_one.GetNumConformers() == 1
    assert only_one.GetConformer(0).GetId() == 0


def test_get_geometry_info_supports_distance_angle_dihedral_and_validation() -> None:
    mol = build_rdmol_with_conformer(
        atom_symbols=["C", "H", "H", "H"],
        bonds=[(0, 1, 1), (0, 2, 1), (0, 3, 1)],
        coordinates=[
            (0.0, 0.0, 0.0),
            (1.0, 0.0, 0.0),
            (0.0, 1.0, 0.0),
            (0.0, 0.0, 1.0),
        ],
    )

    bond_length = get_geometry_info(mol, [0, 1])
    angle = get_geometry_info(mol, [1, 0, 2])
    dihedral = get_geometry_info(mol, [1, 0, 2, 3])

    assert bond_length == pytest.approx(1.0, abs=1e-12)
    assert angle == pytest.approx(90.0, abs=1e-12)
    assert abs(dihedral) == pytest.approx(90.0, abs=1e-12)

    with pytest.raises(ValueError, match="must be 2, 3, or 4"):
        get_geometry_info(mol, [0])


def test_merge_mols_and_is_overlap_invariants_are_deterministic() -> None:
    mol1 = build_rdmol_with_conformer(
        atom_symbols=["O", "H", "H"],
        bonds=[(0, 1, 1), (0, 2, 1)],
        coordinates=[
            (0.0, 0.0, 0.0),
            (0.96, 0.0, 0.0),
            (-0.24, 0.93, 0.0),
        ],
    )
    mol2 = build_rdmol_with_conformer(
        atom_symbols=["O", "H", "H"],
        bonds=[(0, 1, 1), (0, 2, 1)],
        coordinates=[
            (0.0, 0.0, 0.0),
            (0.96, 0.0, 0.0),
            (-0.24, 0.93, 0.0),
        ],
    )

    overlapped = Chem.CombineMols(mol1, mol2)
    assert is_overlap(overlapped) is True

    merged_1 = merge_mols([mol1, mol2], scale=3.0, overlap_threshold=0.9, random_state=123)
    merged_2 = merge_mols([mol1, mol2], scale=3.0, overlap_threshold=0.9, random_state=123)

    assert merged_1.GetNumAtoms() == mol1.GetNumAtoms() + mol2.GetNumAtoms()
    assert merged_2.GetNumAtoms() == merged_1.GetNumAtoms()

    overlap_1 = is_overlap(merged_1, threshold=0.9)
    overlap_2 = is_overlap(merged_2, threshold=0.9)
    assert overlap_1 == overlap_2
    assert overlap_1 is False
