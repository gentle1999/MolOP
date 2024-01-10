"""
Author: TMJ
Date: 2023-10-30 18:21:31
LastEditors: TMJ
LastEditTime: 2024-01-09 19:17:29
Description: 请填写简介
"""

import os

from openbabel import pybel
from rdkit import Chem

from molop.io.bases.file_base import BaseFileParser
from molop.io.bases.molblock_base import BaseBlockParser
from molop.structure.structure import (
    get_bond_pairs,
    get_formal_charges,
    get_formal_spins,
)
from molop.unit import atom_ureg


class SDFBlockParser(BaseBlockParser):
    """
    Parser for XYZ Blocks.
    """

    def __init__(self, block: str):
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


class SDFParser(BaseFileParser):
    """
    Parser for SDF files.
    """

    def __init__(self, file_path: str):
        super().__init__(file_path)
        _, file_format = os.path.splitext(file_path)
        if file_format != ".sdf":
            raise ValueError("File format must be .sdf")
        self._parse()

    def _parse(self):
        for mol in pybel.readfile("sdf", self._file_path):
            self.append(SDFBlockParser(mol.write("sdf")))
