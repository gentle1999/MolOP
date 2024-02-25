"""
Author: TMJ
Date: 2024-02-17 15:17:37
LastEditors: TMJ
LastEditTime: 2024-02-23 20:53:25
Description: 请填写简介
"""
import re

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
        lines = self._block.split("\n")
        num_atoms = int(lines[0])
        cm = re.findall(r"charge\s+([\-0-9]+)\s+multiplicity\s+([\-0-9]+)\s+", lines[1])
        if len(cm) == 1:
            self._charge = int(cm[0][0])
            self._multiplicity = int(cm[0][1])
        for line in lines[2 : 2 + num_atoms]:
            if re.search(r"[A-Za-z]+\s+[\d\.\-]+\s+[\d\.\-]+\s+[\d\.\-]+", line):
                atom, x, y, z = line.split()
                self._atoms.append(atom)
                self._coords.append(
                    (
                        float(x) * atom_ureg.angstrom,
                        float(y) * atom_ureg.angstrom,
                        float(z) * atom_ureg.angstrom,
                    )
                )
