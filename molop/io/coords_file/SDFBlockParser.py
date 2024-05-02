"""
Author: TMJ
Date: 2024-02-17 15:17:37
LastEditors: TMJ
LastEditTime: 2024-02-29 12:21:56
Description: 请填写简介
"""
from openbabel import pybel
from rdkit import Chem

from molop.io.bases.molblock_base import BaseBlockParser
from molop.structure.structure import (
    get_bond_pairs,
    get_formal_charges,
    get_formal_num_radicals,
)
from molop.unit import atom_ureg
from molop.utils import is_metal


class SDFBlockParser(BaseBlockParser):
    """
    Parser for SDF Blocks.
    """

    _block_type = "SDF"

    def __init__(self, block: str, file_path=""):
        super().__init__(block)
        self._file_path = file_path
        self._parse()

    def _parse(self):
        self._rdmol = Chem.MolFromMolBlock(self._block, removeHs=False)
        self._atoms = [atom.GetAtomicNum() for atom in self._rdmol.GetAtoms()]
        self._coords = self._rdmol.GetConformer().GetPositions() * atom_ureg.angstrom
        self._bonds = get_bond_pairs(self._rdmol)
        self._formal_charges = get_formal_charges(self._rdmol)
        self._formal_num_radicals = get_formal_num_radicals(self._rdmol)
        self._charge = sum(self._formal_charges)
        if any(is_metal(atom) for atom in self._atoms):
            self._multiplicity = 3
        else:
            self._multiplicity = sum(self._formal_num_radicals) + 1
