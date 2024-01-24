"""
Author: TMJ
Date: 2024-01-11 09:58:27
LastEditors: TMJ
LastEditTime: 2024-01-24 22:11:19
Description: 请填写简介
"""
from openbabel import pybel
from rdkit import Chem

from molop.io.bases.molblock_base import BaseBlockParser
from molop.structure.structure import (
    get_bond_pairs,
    get_formal_charges,
    get_formal_spins,
)
from molop.unit import atom_ureg


class SDFBlockParser(BaseBlockParser):
    """
    Parser for SDF Blocks.
    """

    _block_type = "SDF"

    def __init__(self, block: str, file_path=""):
        self._file_path = file_path
        super().__init__(block)
        self._parse()

    def _parse(self):
        self._omol = pybel.readstring("sdf", self._block)
        self._rdmol = Chem.MolFromMolBlock(self._block, removeHs=False)
        self._atoms = [atom.atomicnum for atom in self._omol.atoms]
        self._coords = [
            (
                atom.coords[0] * atom_ureg.angstrom,
                atom.coords[1] * atom_ureg.angstrom,
                atom.coords[2] * atom_ureg.angstrom,
            )
            for atom in self._omol.atoms
        ]
        self._bonds = get_bond_pairs(self._rdmol)
        self._formal_charges = get_formal_charges(self._rdmol)
        self._formal_spins = get_formal_spins(self._rdmol)
        self._charge = sum(self._formal_charges)
        self._multiplicity = sum(self._formal_spins) * 2 + 1
