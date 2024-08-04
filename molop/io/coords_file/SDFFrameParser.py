"""
Author: TMJ
Date: 2024-06-20 22:46:56
LastEditors: TMJ
LastEditTime: 2024-08-04 21:01:15
Description: 请填写简介
"""

from rdkit import Chem

from molop.io.bases.BaseMolFrameParser import BaseMolFrameParser
from molop.structure.structure import (
    get_bond_pairs,
    get_formal_charges,
    get_formal_num_radicals,
)
from molop.unit import atom_ureg
from molop.utils.functions import is_metal


class SDFFrameParser(BaseMolFrameParser):
    _frame_type = "SDF"

    def _parse(self):
        self._rdmol = Chem.MolFromMolBlock(self.frame_content, removeHs=False)
        self.atoms = [atom.GetAtomicNum() for atom in self._rdmol.GetAtoms()]
        self.coords = self._rdmol.GetConformer().GetPositions() * atom_ureg.angstrom
        self.bonds = get_bond_pairs(self._rdmol)
        self.formal_charges = get_formal_charges(self._rdmol)
        self.formal_num_radicals = get_formal_num_radicals(self._rdmol)
        self.charge = sum(self.formal_charges)
        if any(is_metal(atom) for atom in self.atoms):
            self._multiplicity = 3
        else:
            self._multiplicity = sum(self.formal_num_radicals) + 1
