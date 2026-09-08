import sys
from pathlib import Path

import pytest
from rdkit import Chem
from rdkit.Chem import BondType

from molop.structure import StructureTransformation as ST
from molop.structure.utils import bond_list, bond_stereo_list
from molop.utils.decorators import ExperimentalWarning


sys.path.append(str(Path(__file__).resolve().parent))

from _helpers_structure import build_rdmol_with_conformer


def _bond_tuple(begin: int, end: int, bond_type: BondType) -> tuple[int, int, int, int]:
    return (
        begin,
        end,
        bond_list.index(bond_type),
        bond_stereo_list.index(Chem.BondStereo.STEREONONE),
    )


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


def test_basic_structure_helpers_cover_bonds() -> None:
    mol = Chem.MolFromSmiles("C=C")
    assert mol is not None

    bond_pairs = ST.get_bond_pairs(mol)
    assert len(bond_pairs) == 1
    assert bond_pairs[0][2] == bond_list.index(BondType.DOUBLE)


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


def test_replace_substituent_rejects_ring_site() -> None:
    mol = Chem.MolFromSmiles("C1CCCCC1")
    assert mol is not None

    ring_bond = mol.GetBondWithIdx(0)
    start = ring_bond.GetBeginAtomIdx()
    end = ring_bond.GetEndAtomIdx()

    with pytest.raises(ValueError, match="ring bonds"):
        ST.replace_substituent(
            mol,
            query="C",
            replacement="[*]C",
            start_idx=start,
            end_idx=end,
        )


def test_multisite_replacement_is_marked_experimental() -> None:
    mol = Chem.MolFromSmiles("CC")
    assert mol is not None
    assert getattr(ST.replace_multisite_substituent, "__experimental__", False) is True

    with (
        pytest.warns(ExperimentalWarning, match="replace_multisite_substituent"),
        pytest.raises(ValueError, match="at least two"),
    ):
        ST.replace_multisite_substituent(
            mol,
            query="C",
            replacement="[*]C",
        )


def test_replace_substituent_requires_unambiguous_single_site() -> None:
    mol = Chem.MolFromSmiles("OCCO")
    assert mol is not None

    with pytest.raises(ValueError, match="multiple pendant substituents"):
        ST.replace_substituent(mol, query="O", replacement="[*]C")


def test_replace_substituent_replace_all_is_finite() -> None:
    mol = Chem.MolFromSmiles("OCCO")
    assert mol is not None

    replaced = ST.replace_substituent(
        mol,
        query="O",
        replacement="[*]C",
        replace_all=True,
        attempt_num=2,
        randomSeed=101,
    )

    assert not replaced.HasSubstructMatch(Chem.MolFromSmarts("O"))
    Chem.SanitizeMol(replaced)


def test_replace_substituent_happy_path_returns_sanitizable_molecule() -> None:
    mol = build_rdmol_with_conformer(
        atom_symbols=["C", "C", "O"],
        bonds=[(0, 1, 1), (1, 2, 1)],
        coordinates=[
            (0.0, 0.0, 0.0),
            (1.54, 0.0, 0.0),
            (2.85, 0.0, 0.0),
        ],
    )

    replaced = ST.replace_substituent(
        mol=mol,
        query="O",
        replacement="[*]C",
        replace_all=False,
        randomSeed=101,
        start_idx=1,
        end_idx=2,
    )

    assert replaced.GetNumAtoms() >= 3
    assert replaced.GetNumBonds() >= 2
    Chem.SanitizeMol(replaced)


@pytest.mark.parametrize("stereo", ["E", "Z"])
def test_replace_substituent_double_bond_records_requested_stereo(stereo: str) -> None:
    mol = Chem.MolFromSmiles("FC=O")
    assert mol is not None

    replaced = ST.replace_substituent(
        mol,
        query="O",
        replacement="[*]=[CH](Cl)",
        start_idx=1,
        end_idx=2,
        attempt_num=2,
        randomSeed=101,
        stereo_policy=stereo,
    )

    bond = replaced.GetBondBetweenAtoms(1, 2)
    assert bond is not None
    assert (
        bond.GetStereo()
        == {
            "E": Chem.BondStereo.STEREOE,
            "Z": Chem.BondStereo.STEREOZ,
        }[stereo]
    )


def test_replace_substituent_accepts_explicit_site_without_query() -> None:
    mol = Chem.MolFromSmiles("CCO")
    assert mol is not None

    replaced = ST.replace_substituent(
        mol,
        query=None,
        replacement="[*]C",
        start_idx=1,
        end_idx=2,
        attempt_num=1,
    )

    assert not replaced.HasSubstructMatch(Chem.MolFromSmarts("O"))
    Chem.SanitizeMol(replaced)
