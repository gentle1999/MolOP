from molop.io.bases.molblock_base import BaseBlockParser
from molop.unit import atom_ureg


import re


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