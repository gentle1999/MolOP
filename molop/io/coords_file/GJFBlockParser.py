'''
Author: TMJ
Date: 2024-01-11 09:58:35
LastEditors: TMJ
LastEditTime: 2024-01-12 14:21:49
Description: 请填写简介
'''
import re

from molop.io.bases.molblock_base import BaseBlockParser
from molop.unit import atom_ureg


class GJFBlockParser(BaseBlockParser):
    """
    Parser for XYZ Blocks.
    """

    def __init__(self, block: str, charge=0, multiplicity=1, parameter_comment=None):
        super().__init__(block)
        self._charge = charge
        self._multiplicity = multiplicity
        self._parameter_comment = parameter_comment
        self._parse()

    def _parse(self):
        """
        Parse the block.
        """
        lines = self._block.split("\n")
        for line in lines[1:]:
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

    @property
    def parameter_comment(self) -> str:
        return self._parameter_comment
