"""Topology and RDKit field conversion helpers."""

from __future__ import annotations

from collections.abc import Sequence

from rdkit import Chem

from molop.utils.types import RdMol

from .contracts import BondRecord
from .utils import bond_list, bond_stereo_list, pt


def get_bond_pairs(mol: RdMol) -> list[BondRecord]:
    """Return RDKit bonds in MolOP's serialized four-integer representation."""
    return [
        (
            bond.GetBeginAtomIdx(),
            bond.GetEndAtomIdx(),
            bond_list.index(bond.GetBondType()),
            bond_stereo_list.index(bond.GetStereo()),
        )
        for bond in mol.GetBonds()
    ]


def get_formal_charges(mol: RdMol) -> list[int]:
    """Return formal charges in atom-index order."""
    return [atom.GetFormalCharge() for atom in mol.GetAtoms()]


def get_formal_num_radicals(mol: RdMol) -> list[int]:
    """Return radical-electron counts in atom-index order."""
    return [atom.GetNumRadicalElectrons() for atom in mol.GetAtoms()]


def get_total_charge(mol: RdMol) -> int:
    """Return the total formal charge."""
    return sum(get_formal_charges(mol))


def get_total_num_radical(mol: RdMol) -> int:
    """Return the total number of radical electrons."""
    return sum(get_formal_num_radicals(mol))


def get_total_multiplicity(mol: RdMol) -> int:
    """Return MolOP's spin-multiplicity convention for a molecule."""
    return get_total_num_radical(mol) + 1


def build_mol_from_atoms_and_bonds(
    atoms: Sequence[int | str],
    bonds: Sequence[BondRecord],
    formal_charges: Sequence[int] | None = None,
    formal_num_radicals: Sequence[int] | None = None,
    *,
    coords: Sequence[tuple[float, float, float]] | None = None,
) -> RdMol:
    """Build an RDKit molecule from MolOP atom, bond, and optional coordinate fields."""
    atom_count = len(atoms)
    if formal_charges is None:
        formal_charges = [0] * atom_count
    if formal_num_radicals is None:
        formal_num_radicals = [0] * atom_count
    if len(formal_charges) != atom_count or len(formal_num_radicals) != atom_count:
        raise ValueError("Atom properties must have the same length as atoms.")

    if coords is None:
        mol = Chem.RWMol()
        for element, formal_charge, formal_num_radical in zip(
            atoms, formal_charges, formal_num_radicals, strict=True
        ):
            atom = Chem.Atom(element if isinstance(element, str) else pt.GetElementSymbol(element))
            atom.SetFormalCharge(formal_charge)
            atom.SetNumRadicalElectrons(formal_num_radical)
            mol.AddAtom(atom)
    else:
        if len(coords) != atom_count:
            raise ValueError("Coordinates must have the same length as atoms.")
        temp_xyz = (
            f"{atom_count}\n\n"
            + "\n".join(
                f"{element if isinstance(element, str) else pt.GetElementSymbol(element)} "
                f"{x:15.6f} {y:15.6f} {z:15.6f}"
                for element, (x, y, z) in zip(atoms, coords, strict=True)
            )
            + "\n"
        )
        xyz_mol = Chem.MolFromXYZBlock(temp_xyz)
        if xyz_mol is None:
            raise ValueError("Invalid input coordinates.")
        mol = Chem.RWMol(xyz_mol)
        for atom_idx, (formal_charge, formal_num_radical) in enumerate(
            zip(formal_charges, formal_num_radicals, strict=True)
        ):
            atom = mol.GetAtomWithIdx(atom_idx)
            atom.SetFormalCharge(formal_charge)
            atom.SetNumRadicalElectrons(formal_num_radical)

    for begin_atom_idx, end_atom_idx, bond_type_idx, _bond_stereo_idx in bonds:
        if not 0 <= begin_atom_idx < atom_count or not 0 <= end_atom_idx < atom_count:
            raise ValueError("Bond atom index is out of range.")
        mol.AddBond(
            beginAtomIdx=begin_atom_idx,
            endAtomIdx=end_atom_idx,
            order=bond_list[bond_type_idx],
        )

    Chem.SanitizeMol(mol, catchErrors=True)
    Chem.SetAromaticity(mol)
    Chem.DetectBondStereochemistry(mol)
    Chem.SetBondStereoFromDirections(mol)
    Chem.AssignStereochemistryFrom3D(mol)
    Chem.AssignCIPLabels(mol)
    for begin_atom_idx, end_atom_idx, _bond_type_idx, bond_stereo_idx in bonds:
        bond = mol.GetBondBetweenAtoms(begin_atom_idx, end_atom_idx)
        if bond is not None:
            bond.SetStereo(bond_stereo_list[bond_stereo_idx])
    return mol.GetMol()


def reset_atom_index(mol: Chem.rdchem.Mol, mapping: Sequence[int]) -> Chem.rdchem.Mol:
    """Move the selected old atom indices to the front, preserving their order."""
    selected = list(mapping)
    atom_count = mol.GetNumAtoms()
    if len(set(selected)) != len(selected) or any(idx < 0 or idx >= atom_count for idx in selected):
        raise ValueError("mapping must contain unique valid atom indices.")
    new_idx = selected + [idx for idx in range(atom_count) if idx not in selected]
    return Chem.RenumberAtoms(mol, new_idx)
