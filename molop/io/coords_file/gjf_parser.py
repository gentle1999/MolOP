"""
Author: TMJ
Date: 2024-01-09 19:54:01
LastEditors: TMJ
LastEditTime: 2024-01-09 20:01:46
Description: 请填写简介
"""
import os
import re

from molop.io.bases.molblock_base import BaseBlockParser
from molop.io.bases.file_base import BaseFileParser
from molop.unit import atom_ureg


class GJFBlockParser(BaseBlockParser):
    """
    Parser for XYZ Blocks.
    """

    def __init__(self, block: str):
        super().__init__(block)
        self._charge = int(block.split("\n")[0].split()[0])
        self._multiplicity = int(block.split("\n")[0].split()[1])
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


class GJFParser(BaseFileParser):
    """
    Parser for GJF files.
    """

    def __init__(self, file_path: str, charge=0, multiplicity=1):
        super().__init__(file_path)
        self._charge = charge
        self._multiplicity = multiplicity
        _, file_format = os.path.splitext(file_path)
        if file_format != ".gjf":
            raise ValueError("File format must be .gjf")
        self._parse()

    def _parse(self):
        """
        Parse the file.
        """
        with open(self.file_path, "r") as fr:
            lines = fr.readlines()
        fr.close()
        for idx, line in enumerate(lines):
            if re.match(r"^\s*\d+\s+\d+$", line):
                block_start = idx
            if re.match(r"^\s*[A-Z][a-z]?(\s+\-?\d+(\.\d+)?){3}$", line):
                block_end = idx
        self.append(GJFBlockParser("".join(lines[block_start : block_end + 1])))
