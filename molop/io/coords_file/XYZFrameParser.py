"""
Author: TMJ
Date: 2024-06-20 21:10:52
LastEditors: TMJ
LastEditTime: 2024-08-04 21:01:45
Description: 请填写简介
"""

import re

import numpy as np
from rdkit import Chem

from molop.io.bases.BaseMolFrameParser import BaseMolFrameParser
from molop.unit import atom_ureg


class XYZFrameParser(BaseMolFrameParser):
    _frame_type = "XYZ"

    def _parse(self):
        """
        Parse the frame.
        """
        lines = self.frame_content.splitlines()
        num_atoms = int(lines[0])
        cm = re.search(
            r"charge\s*([\s-]\d+)\s*multiplicity\s*([\s-]\d+)", self.frame_content
        )
        if cm:
            self.charge = int(cm.group(1))
            self.multiplicity = int(cm.group(2))
        temp_coords = []
        for line in lines[2 : 2 + num_atoms]:
            matches = re.search(
                r"([A-Za-z]+)\s*([\d\.\-]+)\s+([\d\.\-]+)\s+([\d\.\-]+)", line
            )
            atom, x, y, z = matches.groups()
            self.atoms.append(Chem.Atom(atom).GetAtomicNum())
            temp_coords.append((float(x), float(y), float(z)))
        self.coords = np.array(temp_coords) * atom_ureg.angstrom
