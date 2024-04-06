"""
Author: TMJ
Date: 2024-02-17 15:17:37
LastEditors: TMJ
LastEditTime: 2024-03-31 16:05:02
Description: 请填写简介
"""
import re
import numpy as np

from molop.io.bases.molblock_base import BaseBlockParser
from molop.unit import atom_ureg


class XYZBlockParser(BaseBlockParser):
    """
    Parser for XYZ Blocks.
    """

    _block_type = "XYZ"

    def __init__(
        self,
        block: str,
        charge=0,
        multiplicity=1,
        file_path="",
    ):
        super().__init__(block)
        self._file_path = file_path
        self._charge = charge
        self._multiplicity = multiplicity
        self._parse()

    def _parse(self):
        """
        Parse the block.
        """
        lines = self._block.splitlines()
        num_atoms = int(lines[0])
        cm = re.findall(r"charge\s+([\-0-9]+)\s+multiplicity\s+([\-0-9]+)\s+", lines[1])
        if len(cm) == 1:
            self._charge = int(cm[0][0])
            self._multiplicity = int(cm[0][1])
        temp_coords = []
        for line in lines[2 : 2 + num_atoms]:
            if re.search(r"[A-Za-z]+\s+[\d\.\-]+\s+[\d\.\-]+\s+[\d\.\-]+", line):
                atom, x, y, z = line.split()
                self._atoms.append(atom)
                temp_coords.append((float(x), float(y), float(z)))
        self._coords = np.array(temp_coords) * atom_ureg.angstrom
