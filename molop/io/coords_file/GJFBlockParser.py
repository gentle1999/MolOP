from molop.io.bases.molblock_base import BaseBlockParser
from molop.unit import atom_ureg


import re


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