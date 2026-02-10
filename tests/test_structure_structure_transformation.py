import sys
from pathlib import Path

import pytest
from rdkit import Chem
from rdkit.Chem import BondType
from rdkit.Geometry import Point3D

from molop.structure import StructureTransformation as ST
from molop.structure.utils import bond_list, bond_stereo_list


sys.path.append(str(Path(__file__).resolve().parent))

from _helpers_structure import build_rdmol_with_conformer


def _bond_tuple(begin: int, end: int, bond_type: BondType) -> tuple[int, int, int, int]:
    return (
        begin,
        end,
        bond_list.index(bond_type),
        bond_stereo_list.index(Chem.BondStereo.STEREONONE),
    )


def _radical_mol(radicals: list[int]) -> Chem.Mol:
    rw_mol = Chem.RWMol()
    for _ in radicals:
        rw_mol.AddAtom(Chem.Atom("C"))
    for idx, num_radical in enumerate(radicals):
        rw_mol.GetAtomWithIdx(idx).SetNumRadicalElectrons(num_radical)
    mol = rw_mol.GetMol()
    conf = Chem.Conformer(len(radicals))
    for idx in range(len(radicals)):
        conf.SetAtomPosition(idx, Point3D(float(idx), 0.0, 0.0))
    mol.AddConformer(conf, assignId=True)
    return mol


def _count_bonds_by_type(mol: Chem.Mol, bond_type: BondType) -> int:
    return sum(1 for bond in mol.GetBonds() if bond.GetBondType() == bond_type)


def test_formal_helpers_return_expected_values() -> None:
    mol = Chem.MolFromSmiles("[CH2][O-]")
    assert mol is not None

    charges = ST.get_formal_charges(mol)
    radicals = ST.get_formal_num_radicals(mol)

    assert charges == [0, -1]
    assert radicals == [1, 0]
    assert ST.get_total_charge(mol) == -1
    assert ST.get_total_num_radical(mol) == 1
    assert ST.get_total_multiplicity(mol) == 2


def test_basic_structure_helpers_cover_bonds_resonance_and_equality() -> None:
    mol = Chem.MolFromSmiles("C=C")
    assert mol is not None

    bond_pairs = ST.get_bond_pairs(mol)
    assert len(bond_pairs) == 1
    assert bond_pairs[0][2] == bond_list.index(BondType.DOUBLE)

    resonance = list(ST.get_resonance_structures(Chem.MolFromSmiles("[O-][N+](=O)O")))
    assert len(resonance) >= 1

    assert ST.structure_score(Chem.MolFromSmiles("[CH2][O-]")) == 2
    assert ST.check_mol_equal(Chem.MolFromSmiles("CC"), Chem.MolFromSmiles("CC")) is True


def test_build_mol_from_atoms_and_bonds_without_coords() -> None:
    mol = ST.build_mol_from_atoms_and_bonds(
        atoms=["C", "O"],
        bonds=[_bond_tuple(0, 1, BondType.SINGLE)],
        formal_charges=[0, -1],
        formal_num_radicals=[0, 0],
        coords=None,
    )

    assert mol.GetNumAtoms() == 2
    assert mol.GetNumBonds() == 1
    assert mol.GetNumConformers() == 0
    assert ST.get_total_charge(mol) == -1
    Chem.SanitizeMol(mol)


def test_build_mol_from_atoms_and_bonds_with_valid_coords() -> None:
    mol = ST.build_mol_from_atoms_and_bonds(
        atoms=["C", "O"],
        bonds=[_bond_tuple(0, 1, BondType.DOUBLE)],
        formal_charges=[0, 0],
        formal_num_radicals=[0, 0],
        coords=[(0.0, 0.0, 0.0), (1.2, 0.0, 0.0)],
    )

    assert mol.GetNumAtoms() == 2
    assert mol.GetNumConformers() == 1
    assert mol.GetBondBetweenAtoms(0, 1) is not None
    assert mol.GetBondBetweenAtoms(0, 1).GetBondType() == BondType.DOUBLE
    Chem.SanitizeMol(mol)


def test_build_mol_from_atoms_and_bonds_invalid_coords_raise_value_error(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(ST.Chem, "MolFromXYZBlock", lambda _xyz: None)
    monkeypatch.setattr(ST.Chem, "RWMol", lambda mol: mol)

    with pytest.raises(ValueError, match="Invalid input coordinates\\."):
        ST.build_mol_from_atoms_and_bonds(
            atoms=["C", "O"],
            bonds=[_bond_tuple(0, 1, BondType.SINGLE)],
            coords=[(0.0, 0.0, 0.0), (1.2, 0.0, 0.0)],
        )


@pytest.mark.parametrize(
    ("bond_type", "required_radicals"),
    [
        (BondType.SINGLE, 1),
        (BondType.DOUBLE, 2),
        (BondType.TRIPLE, 3),
    ],
)
def test_transform_replacement_index_radical_selection(
    bond_type: BondType, required_radicals: int
) -> None:
    mol = _radical_mol([0, required_radicals, required_radicals + 1])

    by_relative = ST.transform_replacement_index(mol, bond_tag=bond_type, relative_idx=0)
    by_absolute = ST.transform_replacement_index(mol, bond_tag=bond_type, absolute_idx=2)

    assert by_relative.GetAtomWithIdx(0).GetNumRadicalElectrons() >= required_radicals
    assert by_absolute.GetAtomWithIdx(0).GetNumRadicalElectrons() == required_radicals + 1


@pytest.mark.parametrize("bond_type", [BondType.SINGLE, BondType.DOUBLE, BondType.TRIPLE])
def test_transform_replacement_index_radical_error_cases(bond_type: BondType) -> None:
    mol = _radical_mol([0, 1, 2, 3])

    with pytest.raises(ValueError, match="Relative index is out of range\\."):
        ST.transform_replacement_index(mol, bond_tag=bond_type, relative_idx=99)

    with pytest.raises(ValueError, match="Absolute index is not a radical atom\\."):
        ST.transform_replacement_index(mol, bond_tag=bond_type, absolute_idx=0)


def test_transform_replacement_index_dative_selection_and_errors() -> None:
    mol = _radical_mol([0, 0, 0])

    by_relative = ST.transform_replacement_index(mol, bond_tag=BondType.DATIVE, relative_idx=1)
    by_absolute = ST.transform_replacement_index(mol, bond_tag=BondType.DATIVE, absolute_idx=2)

    assert by_relative.GetAtomWithIdx(0).GetIdx() == 0
    assert by_absolute.GetAtomWithIdx(0).GetIdx() == 0

    with pytest.raises(ValueError, match="Relative index is out of range\\."):
        ST.transform_replacement_index(mol, bond_tag=BondType.DATIVE, relative_idx=10)

    with pytest.raises(ValueError, match="Absolute index is not a legal atom\\."):
        ST.transform_replacement_index(mol, bond_tag=BondType.DATIVE, absolute_idx=10)


def test_transform_replacement_index_rejects_unsupported_bond_type() -> None:
    mol = _radical_mol([1])
    with pytest.raises(ValueError, match="Unsupported bond type"):
        ST.transform_replacement_index(mol, bond_tag=BondType.QUADRUPLE)


def test_get_skeleton_raises_when_query_is_ring_bond() -> None:
    mol = Chem.MolFromSmiles("C1CCCCC1")
    assert mol is not None

    ring_bond = mol.GetBondWithIdx(0)
    start = ring_bond.GetBeginAtomIdx()
    end = ring_bond.GetEndAtomIdx()

    with pytest.raises(RuntimeError, match="query_mol should not be in a ring"):
        ST.get_skeleton(Chem.RWMol(mol), start, end)


def test_attempt_replacement_raises_on_replace_all_endless_loop_guard() -> None:
    mol = Chem.MolFromSmiles("CC")
    assert mol is not None

    with pytest.raises(RuntimeError, match="Endless loop"):
        ST.attempt_replacement(mol, query="C", replacement="CC", replace_all=True)


def test_attempt_replacement_happy_path_returns_sanitizable_molecule() -> None:
    mol = build_rdmol_with_conformer(
        atom_symbols=["C", "C", "O"],
        bonds=[(0, 1, 1), (1, 2, 1)],
        coordinates=[
            (0.0, 0.0, 0.0),
            (1.54, 0.0, 0.0),
            (2.85, 0.0, 0.0),
        ],
    )

    replaced = ST.attempt_replacement(
        mol=mol,
        query="O",
        replacement="[CH3]",
        replace_all=False,
        randomSeed=101,
        start_idx=1,
        end_idx=2,
    )

    assert replaced.GetNumAtoms() >= 3
    assert replaced.GetNumBonds() >= 2
    Chem.SanitizeMol(replaced)


def test_make_dative_bonds_no_metal_returns_molecule_without_crash() -> None:
    water = build_rdmol_with_conformer(
        atom_symbols=["O", "H", "H"],
        bonds=[(0, 1, 1), (0, 2, 1)],
        coordinates=[
            (0.0, 0.0, 0.0),
            (0.96, 0.0, 0.0),
            (-0.24, 0.93, 0.0),
        ],
    )

    result = ST.make_dative_bonds(Chem.RWMol(water))

    assert isinstance(result, Chem.RWMol)
    assert result.GetNumAtoms() == water.GetNumAtoms()
    Chem.SanitizeMol(result)


def test_make_dative_bonds_adds_dative_for_minimal_metal_and_donor_case() -> None:
    rw_mol = Chem.RWMol()
    rw_mol.AddAtom(Chem.Atom("C"))
    rw_mol.AddAtom(Chem.Atom("O"))
    rw_mol.AddAtom(Chem.Atom("Fe"))
    rw_mol.AddBond(0, 1, BondType.DOUBLE)

    mol = rw_mol.GetMol()
    conf = Chem.Conformer(3)
    conf.SetAtomPosition(0, Point3D(0.0, 0.0, 0.0))
    conf.SetAtomPosition(1, Point3D(1.20, 0.0, 0.0))
    conf.SetAtomPosition(2, Point3D(2.80, 0.0, 0.0))
    mol.AddConformer(conf, assignId=True)

    datived = ST.make_dative_bonds(Chem.RWMol(mol))

    assert datived.GetNumAtoms() == 3
    assert _count_bonds_by_type(datived, BondType.DATIVE) >= 1
