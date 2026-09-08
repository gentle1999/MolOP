"""
Author: TMJ
Date: 2026-07-10 00:00:00
LastEditors: TMJ
LastEditTime: 2026-07-10 00:00:00
Description: Coordinate file extraction helpers
"""

from typing import Any, cast

import numpy as np
from rdkit import Chem
from rdkit.Chem.rdDepictor import Compute2DCoords

from molop.io.codec_exceptions import FormatMismatchError
from molop.io.logic.coords.frame_parsers._xyz_patterns import xyz_patterns
from molop.structure.topology import (
    get_bond_pairs,
    get_formal_charges,
    get_formal_num_radicals,
)
from molop.unit import atom_ureg


pt = Chem.GetPeriodicTable()


def extract_xyz_frame_payload(block: str) -> dict[str, Any]:
    charge = 0
    multiplicity = 1
    lines = block.splitlines()
    if len(lines) < 3:
        raise FormatMismatchError("Not an XYZ frame: missing header or coordinate lines.")
    try:
        atom_num = int(lines[0].strip())
    except ValueError as exc:
        raise FormatMismatchError("Not an XYZ frame: first line is not an atom count.") from exc
    comment = lines[1].strip()
    if matched := xyz_patterns.CHARGE_MULTIPLICITY.match(comment):
        charge = int(matched.group("charge"))
        multiplicity = int(matched.group("multiplicity"))
    if matches := xyz_patterns.ATOMS.find_matches(block):
        if len(matches) != atom_num:
            raise FormatMismatchError("Not an XYZ frame: atom number does not match.")
        atoms = [pt.GetAtomicNumber(matched.group("symbol")) for matched in matches]
        coords = np.array(
            [
                (
                    float(matched.group("x")),
                    float(matched.group("y")),
                    float(matched.group("z")),
                )
                for matched in matches
            ],
            dtype=np.float32,
        )
        return {
            "comment": comment,
            "atoms": atoms,
            "coords": coords * atom_ureg.angstrom,
            "charge": charge,
            "multiplicity": multiplicity,
        }
    raise FormatMismatchError("Not an XYZ frame: no valid atom coordinates found.")


def extract_sdf_frame_payload(block: str) -> dict[str, Any]:
    suppl = Chem.SDMolSupplier()
    suppl.SetData(block, removeHs=False, sanitize=False)
    fake_mol: Chem.Mol = next(suppl)
    if fake_mol is None:
        raise FormatMismatchError("Not an SDF/MOL frame: invalid mol block.")
    formal_charges = get_formal_charges(fake_mol)
    formal_num_radicals = get_formal_num_radicals(fake_mol)
    return {
        "atoms": [cast(Chem.Atom, atom).GetAtomicNum() for atom in fake_mol.GetAtoms()],
        "coords": fake_mol.GetConformer().GetPositions() * atom_ureg.angstrom,
        "charge": sum(formal_charges),
        "multiplicity": sum(formal_num_radicals) + 1,
        "bonds": get_bond_pairs(fake_mol),
        "formal_charges": formal_charges,
        "formal_num_radicals": formal_num_radicals,
    }


def extract_smi_frame_payload(block: str) -> dict[str, Any]:
    smiles = block.strip().split()[0] if block.strip() else ""
    rdmol = Chem.MolFromSmiles(smiles)
    if rdmol is None:
        raise FormatMismatchError("Not a SMILES frame: invalid SMILES token.")
    Compute2DCoords(rdmol)
    coords = rdmol.GetConformer().GetPositions()
    formal_charges = get_formal_charges(rdmol)
    formal_num_radicals = get_formal_num_radicals(rdmol)
    return {
        "atoms": [cast(Chem.Atom, atom).GetAtomicNum() for atom in rdmol.GetAtoms()],
        "coords": coords * atom_ureg.angstrom,
        "charge": sum(formal_charges),
        "multiplicity": sum(formal_num_radicals) + 1,
        "bonds": get_bond_pairs(rdmol),
        "formal_charges": formal_charges,
        "formal_num_radicals": formal_num_radicals,
    }
