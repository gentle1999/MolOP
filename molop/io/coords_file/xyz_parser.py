"""
Author: TMJ
Date: 2023-10-30 15:41:19
LastEditors: TMJ
LastEditTime: 2023-12-10 16:54:10
Description: 请填写简介
"""
import os
import re

from molop.io.bases.molblock_base import BaseBlockParser
from molop.io.bases.file_base import BaseFileParser
from molop.unit import atom_ureg


class XYZBlockParser(BaseBlockParser):
    """
    Parser for XYZ Blocks.
    """

    def __init__(self, block: str, charge=0, multiplicity=1):
        super().__init__(block)
        self._charge = charge
        self._multiplicity = multiplicity
        self._parse()

    def _parse(self):
        """
        Parse the block.
        """
        lines = self._block.split("\n")
        num_atoms = int(lines[0])
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


class XYZParser(BaseFileParser):
    """
    Parser for XYZ files.
    Supports one or more molecules in one file.
    Make sure molecules in the file have same charge and multiplicity.
    """

    def __init__(self, file_path: str, charge=0, multiplicity=1):
        super().__init__(file_path)
        self._charge = charge
        self._multiplicity = multiplicity
        _, file_format = os.path.splitext(file_path)
        if file_format != ".xyz":
            raise ValueError("File format must be .xyz")
        self._parse()

    def _parse(self):
        """
        Parse the file.
        """
        with open(self.file_path, "r") as fr:
            lines = fr.readlines()
        fr.close()
        anchor = 0
        while anchor < len(lines):
            num_atoms = int(lines[anchor])
            block = XYZBlockParser(
                "".join(lines[anchor : anchor + num_atoms + 2]),
                charge=self._charge,
                multiplicity=self._multiplicity,
            )
            self.append(block)
            anchor += num_atoms + 2
