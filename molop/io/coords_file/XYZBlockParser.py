"""
Author: TMJ
Date: 2024-01-11 09:58:01
LastEditors: TMJ
LastEditTime: 2024-01-27 15:08:01
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
        cm = re.findall("charge\s+([\-0-9]+)\s+multiplicity\s+([\-0-9]+)\s+", lines[1])
        if len(cm) == 1:
            self._charge = int(cm[0][0])
            self._multiplicity = int(cm[0][1])
        for line in lines[2 : 2 + num_atoms]:
            if re.match(r"^\s*[A-Z][a-z]?(\s+\-?\d+(\.\d+)?){3}$", line):
                atom, x, y, z = line.split()
                self._atoms.append(atom)
                self._coords.append(
                    (
                        float(x) * atom_ureg.angstrom,
                        float(y) * atom_ureg.angstrom,
                        float(z) * atom_ureg.angstrom,
                    )
                )
