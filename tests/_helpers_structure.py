from __future__ import annotations

from rdkit import Chem
from rdkit.Chem import rdchem
from rdkit.Geometry import Point3D


def build_xyz_block(atoms: list[tuple[str, float, float, float]], comment: str = "") -> str:
    lines = [str(len(atoms)), comment]
    lines.extend(f"{symbol} {x:.6f} {y:.6f} {z:.6f}" for symbol, x, y, z in atoms)
    return "\n".join(lines) + "\n"


def build_water_xyz_block() -> str:
    return build_xyz_block(
        [
            ("O", 0.000000, 0.000000, 0.000000),
            ("H", 0.957200, 0.000000, 0.000000),
            ("H", -0.239987, 0.927297, 0.000000),
        ],
        comment="water",
    )


def build_methane_xyz_block() -> str:
    return build_xyz_block(
        [
            ("C", 0.000000, 0.000000, 0.000000),
            ("H", 0.629118, 0.629118, 0.629118),
            ("H", -0.629118, -0.629118, 0.629118),
            ("H", -0.629118, 0.629118, -0.629118),
            ("H", 0.629118, -0.629118, -0.629118),
        ],
        comment="methane",
    )


def build_rdmol_with_conformer(
    atom_symbols: list[str],
    bonds: list[tuple[int, int, int]],
    coordinates: list[tuple[float, float, float]],
) -> Chem.Mol:
    if len(atom_symbols) != len(coordinates):
        raise ValueError("atom_symbols and coordinates must have the same length")

    rw_mol = Chem.RWMol()
    for symbol in atom_symbols:
        rw_mol.AddAtom(Chem.Atom(symbol))

    bond_type_by_order = {
        1: rdchem.BondType.SINGLE,
        2: rdchem.BondType.DOUBLE,
        3: rdchem.BondType.TRIPLE,
    }
    for begin_idx, end_idx, bond_order in bonds:
        rw_mol.AddBond(begin_idx, end_idx, bond_type_by_order[bond_order])

    mol = rw_mol.GetMol()
    conformer = Chem.Conformer(len(atom_symbols))
    for atom_idx, (x, y, z) in enumerate(coordinates):
        conformer.SetAtomPosition(atom_idx, Point3D(x, y, z))
    mol.AddConformer(conformer, assignId=True)
    return mol


def build_water_rdmol() -> Chem.Mol:
    return build_rdmol_with_conformer(
        atom_symbols=["O", "H", "H"],
        bonds=[(0, 1, 1), (0, 2, 1)],
        coordinates=[
            (0.000000, 0.000000, 0.000000),
            (0.957200, 0.000000, 0.000000),
            (-0.239987, 0.927297, 0.000000),
        ],
    )


def build_methane_rdmol() -> Chem.Mol:
    return build_rdmol_with_conformer(
        atom_symbols=["C", "H", "H", "H", "H"],
        bonds=[(0, 1, 1), (0, 2, 1), (0, 3, 1), (0, 4, 1)],
        coordinates=[
            (0.000000, 0.000000, 0.000000),
            (0.629118, 0.629118, 0.629118),
            (-0.629118, -0.629118, 0.629118),
            (-0.629118, 0.629118, -0.629118),
            (0.629118, -0.629118, -0.629118),
        ],
    )
