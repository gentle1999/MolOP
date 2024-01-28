import re
from typing import Literal

from packaging.version import Version

from molop.io.bases.molblock_base import QMBaseBlockParser
from molop.logger.logger import logger
from molop.unit import atom_ureg


class XTBOUTBlockParser(QMBaseBlockParser):
    """
    Parser for xTB out Blocks.
    """

    _block_type = "XTB OUT"

    def __init__(
        self,
        block: str,
        charge=0,
        multiplicity=1,
        file_path="",
        version=None,
        parameter_comment=None,
        only_extract_structure=False,
    ):
        super().__init__(block, only_extract_structure)
        self._file_path = file_path
        self._charge = charge
        self._multiplicity = multiplicity
        self._version = version
        self._parameter_comment = parameter_comment
        try:
            self.__n_atom = int(
                re.findall(
                    "number of atoms\s+\:\s+(\d+)",
                    self._block,
                )[0]
            )
        except:
            self.__n_atom = int(
                re.findall(
                    "final structure:\s+\=+\s+(\d+)",
                    self._block,
                )[0]
            )
        self._parse_coords()
        if not self._only_extract_structure:
            self._parse()

    def _parse(self):
        self._parse_state()
        self._parse_energy()
        self._parse_partial_charges()
        self._parse_orbitals()

    def _parse_coords(self):
        if self._version is None:
            version = Version("0.0.0")
        else:
            version = Version(self._version.split()[1])
        if version <= Version("6.2.2"):
            return self._parse_coords_old()
        else:
            return self._parse_coords_new()

    def _parse_state(self):
        if re.search("GEOMETRY OPTIMIZATION CONVERGED", self._block):
            self._state["geometric_optimization"] = True
        else:
            self._state["geometric_optimization"] = False

    def _parse_energy(self):
        lines = self._block.splitlines()
        for line in reversed(lines):
            if "total E" in line:
                self._energy = (
                    round(float(line.split()[-1]), 6)
                    * atom_ureg.hartree
                    / atom_ureg.particle
                )
                return
            if "TOTAL ENERGY" in line:
                self._energy = (
                    round(float(line.split()[-3]), 6)
                    * atom_ureg.hartree
                    / atom_ureg.particle
                )
                return
        raise ValueError("Energy not found")

    def _parse_partial_charges(self):
        lines = self._block.splitlines()
        charges_sect = False
        charges = []
        for line in lines:
            if "Mol." in line:
                charges_sect = False
            if charges_sect and len(line.split()) == 7:
                charges.append(float(line.split()[4]))
            if "covCN" in line:
                charges_sect = True
        self._partial_charges = charges

    def _parse_orbitals(self):
        self.alpha_energy["homo"] = (
            round(float(re.findall("([\+\-0-9.]+)\s+\(HOMO\)", self._block)[-1]), 6)
            * atom_ureg.eV
            / atom_ureg.particle
        ).to("hartree/particle")
        self.alpha_energy["lumo"] = (
            round(float(re.findall("([\+\-0-9.]+)\s+\(LUMO\)", self._block)[-1]), 6)
            * atom_ureg.eV
            / atom_ureg.particle
        ).to("hartree/particle")
        self.alpha_energy["gap"] = self.alpha_energy["lumo"] - self.alpha_energy["homo"]

    def _parse_coords_old(self):
        """
        e.g.

        ================
         final structure:
        ================
        $coord
            2.52072290250473   -0.04782551206377   -0.50388676977877      C
                    .                 .                    .              .
        """
        lines = self._block.splitlines()
        for i, line in enumerate(lines):
            if "final structure" in line:
                for j in range(i + 3, i + self.__n_atom + 3):
                    line = lines[j]
                    x, y, z, atom = line.split()
                    self._atoms.append(atom)
                    self._coords.append(
                        (
                            (round(float(x), 6) * atom_ureg.bohr).to("angstrom"),
                            (round(float(y), 6) * atom_ureg.bohr).to("angstrom"),
                            (round(float(z), 6) * atom_ureg.bohr).to("angstrom"),
                        )
                    )

    def _parse_coords_new(self):
        """
        e.g.

        ================
         final structure:
        ================
        5
         xtb: 6.2.3 (830e466)
        Cl        1.62694523673790    0.09780349799138   -0.02455489507427
        C        -0.15839164427314   -0.00942638308615    0.00237760557913
        H        -0.46867957388620   -0.59222865914178   -0.85786049981721
        H        -0.44751262498645   -0.49575975568264    0.92748366742968
        H        -0.55236139359212    0.99971129991918   -0.04744587811734
        """
        lines = self._block.splitlines()
        for i, line in enumerate(lines):
            if "final structure" in line:
                for j in range(i + 4, i + self.__n_atom + 4):
                    line = lines[j]
                    atom, x, y, z = line.split()
                    self._atoms.append(atom)
                    self._coords.append(
                        (
                            round(float(x), 6) * atom_ureg.angstrom,
                            round(float(y), 6) * atom_ureg.angstrom,
                            round(float(z), 6) * atom_ureg.angstrom,
                        )
                    )
