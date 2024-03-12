import re
from typing import Literal

import numpy as np
from packaging.version import Version

from molop.io.bases.molblock_base import QMBaseBlockParser
from molop.logger.logger import logger
from molop.unit import atom_ureg
from molop.utils import xtboutpatterns


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
        self._parse_coords()
        if not self._only_extract_structure:
            self._parse()

    def _parse(self):
        self._parse_state()
        self._parse_energy()
        self._parse_partial_charges()
        self._parse_orbitals()
        self._wiberg_bond_order = self._parse_bond_order()
        self._parse_dipole()
        self._parse_thermo()

    def _parse_coords(self):
        coords_match = re.search(xtboutpatterns["coords_start"], self._block)
        if not coords_match:
            logger.error(
                "No coordinates found. If you are using `--hess`, use `--ohess` instead. See https://xtb-docs.readthedocs.io/en/latest/hessian.html for more details."
            )
            raise RuntimeError("No coordinates found.")
        coords_end = re.search(xtboutpatterns["coords_end"], self._block)
        coords_block = self._block[coords_match.start() : coords_end.end()]
        if self._version is None:
            version = Version("0.0.0")
        else:
            version = Version(self._version.split()[1])
        if version <= Version("6.2.2"):
            self._parse_coords_old(coords_block)
        else:
            self._parse_coords_new(coords_block)
        self.__n_atom = len(self.atoms)

    def _parse_state(self):
        if re.search(xtboutpatterns["state"], self._block):
            self._state["geometric_optimization"] = True
        else:
            self._state["geometric_optimization"] = False

    def is_error(self) -> bool:
        if "geometric_optimization" in self._state:
            return self._state["geometric_optimization"] == False
        else:
            return True

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
        if not re.search(xtboutpatterns["homo"], self._block):
            return
        self._alpha_energy["homo"] = (
            round(float(re.findall(xtboutpatterns["homo"], self._block)[-1]), 6)
            * atom_ureg.eV
            / atom_ureg.particle
        ).to("hartree/particle")
        self._alpha_energy["lumo"] = (
            round(float(re.findall(xtboutpatterns["lumo"], self._block)[-1]), 6)
            * atom_ureg.eV
            / atom_ureg.particle
        ).to("hartree/particle")
        self._alpha_energy["gap"] = (
            round((self.alpha_energy["lumo"] - self.alpha_energy["homo"]).m, 6)
            * atom_ureg.eV
            / atom_ureg.particle
        )

    def _parse_coords_old(self, block: str):
        """
        e.g.

        ================
         final structure:
        ================
        $coord
            2.52072290250473   -0.04782551206377   -0.50388676977877      C
                    .                 .                    .              .
        """
        coords = xtboutpatterns["coords_old"].findall(block)
        temp_coords = []
        for x, y, z, atom in coords:
            temp_coords.append((float(x), float(y), float(z)))
            self._atoms.append(atom)
        self._coords = (np.array(temp_coords) * atom_ureg.bohr).to("angstrom")

    def _parse_coords_new(self, block: str):
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
        coords = xtboutpatterns["coords_new"].findall(block)
        temp_coords = []
        for atom, x, y, z in coords:
            temp_coords.append((float(x), float(y), float(z)))
            self._atoms.append(atom)
        self._coords = np.array(temp_coords) * atom_ureg.angstrom

    def _parse_bond_order(self):
        start_matches = xtboutpatterns[f"wiberg_start"].search(self._block)
        end_matches = xtboutpatterns[f"wiberg_end"].search(self._block)
        result = np.zeros((self.__n_atom, self.__n_atom))
        if start_matches:
            temp_block = self._block[start_matches.start() : end_matches.end()]
            for idx, sub_temp_block in enumerate(
                xtboutpatterns["index"].split(temp_block)[1:]
            ):
                bond = list(
                    map(
                        lambda x: int(x) - 1,
                        xtboutpatterns["bond"].findall(sub_temp_block),
                    )
                )
                value = list(
                    map(float, xtboutpatterns["value"].findall(sub_temp_block))
                )[1:]
                for i, v in zip(bond, value):
                    result[idx][i] = v
                    result[i][idx] = v
            return result
        return None

    def _parse_dipole(self):
        dipole_start = xtboutpatterns["dipole_start"].search(self._block)
        dipole_end = xtboutpatterns["dipole_end"].search(self._block)
        if dipole_start and dipole_end:
            temp_block = self._block[dipole_start.start() : dipole_end.end()]
            values = list(map(float, xtboutpatterns["value"].findall(temp_block)))
            self._dipole = (np.array(values[:3]) + np.array(values[3:6])).astype(
                np.float32
            ) * atom_ureg.debye

    def _parse_thermo(self):
        zpc = xtboutpatterns["zp correct"].search(self._block)
        if zpc:
            self._sum_energy["zero-point correction"] = (
                float(zpc.group(1)) * atom_ureg.hartree / atom_ureg.particle
            )
        hc = xtboutpatterns["H correct"].search(self._block)
        if hc:
            self._sum_energy["TCH"] = (
                float(hc.group(1)) * atom_ureg.hartree / atom_ureg.particle
            )
        gc = xtboutpatterns["G correct"].search(self._block)
        if gc:
            self._sum_energy["TCG"] = (
                float(gc.group(1)) * atom_ureg.hartree / atom_ureg.particle
            )
        se = xtboutpatterns["sum E"].search(self._block)
        if se:
            self._sum_energy["zero-point sum"] = (
                float(se.group(1)) * atom_ureg.hartree / atom_ureg.particle
            )
        sh = xtboutpatterns["sum H"].search(self._block)
        if sh:
            self._sum_energy["H sum"] = (
                float(sh.group(1)) * atom_ureg.hartree / atom_ureg.particle
            )
        sg = xtboutpatterns["sum G"].search(self._block)
        if sg:
            self._sum_energy["G sum"] = (
                float(sg.group(1)) * atom_ureg.hartree / atom_ureg.particle
            )
