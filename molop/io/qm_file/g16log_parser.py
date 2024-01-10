import os
import re
from typing import Literal

from molop.io.bases.molblock_base import QMBaseBlockParser
from molop.io.bases.file_base import BaseQMFileParser
from molop.unit import atom_ureg
from molop.logger import logger


class G16LOGBlockParser(QMBaseBlockParser):
    """
    Parser for G16 log Blocks.
    """

    def __init__(
        self, block: str, charge=0, multiplicity=1, n_atom=1, parameter_comment=None
    ):
        super().__init__(block)
        self._charge = charge
        self._multiplicity = multiplicity
        self.__n_atom = n_atom
        self._parameter_comment = parameter_comment
        self._parse()

    def _parse(self):
        self._parse_coords()
        self._parse_energy()
        self._parse_partial_charges()
        self._parse_gradient()
        self._parse_orbitals("Alpha")
        self._parse_orbitals("Beta")
        self._parse_frequencies()
        self._parse_sum_energy()
        self._parse_hessian()
        self._parse_state()

    def _parse_coords(self):
        lines = self._block.splitlines()
        for i, line in enumerate(lines):
            if "Input orientation" in line:
                xyz_lines = lines[i + 5 : i + 5 + self.__n_atom]

                for xyz_line in xyz_lines:
                    _, atom_num, _, x, y, z = xyz_line.split()
                    self._atoms.append(int(atom_num))
                    self._coords.append(
                        (
                            float(x) * atom_ureg.angstrom,
                            float(y) * atom_ureg.angstrom,
                            float(z) * atom_ureg.angstrom,
                        )
                    )

    def _parse_energy(self):
        lines = self._block.splitlines()
        for line in reversed(lines):
            if "SCF Done" in line or "E(CIS)" in line:
                self._energy = (
                    float(line.split()[4]) * atom_ureg.hartree / atom_ureg.particle
                )
                self._state["SCF Done"] = True
                return
            if "E(CORR)" in line or "E(CI)" in line:
                self._energy = (
                    float(line.split()[3]) * atom_ureg.hartree / atom_ureg.particle
                )
                return
            if "E(CIS(D))" in line:
                self._energy = (
                    float(line.split()[5]) * atom_ureg.hartree / atom_ureg.particle
                )
                return
            if line.startswith(" Energy=") and "NIter=" in line:
                self._energy = (
                    float(line.split()[1]) * atom_ureg.hartree / atom_ureg.particle
                )
                return

    def _parse_partial_charges(self):
        lines = self._block.splitlines()
        charges_section = False
        charges = []
        for line in reversed(lines):
            if len(charges) == self.__n_atom:
                self._partial_charges = list(reversed(charges))
                break
            if charges_section:
                charges.append(float(line.split()[2]))
            if "sum of mulliken charges" in line.lower():
                charges_section = True

    def _parse_gradient(self):
        lines = self._block.splitlines()
        raw_gradient = []

        for i, line in enumerate(lines):
            if "Forces (Hartrees/Bohr)" not in line:
                continue
            raw_gradient = []  # NOTE: possibly multiple gradients in a file
            for force_line in lines[i + 3 : i + 3 + self.__n_atom]:
                try:
                    _, _, fx, fy, fz = force_line.split()
                    raw_gradient.append(
                        (
                            float(fx) * atom_ureg.hartree / atom_ureg.bohr,
                            float(fy) * atom_ureg.hartree / atom_ureg.bohr,
                            float(fz) * atom_ureg.hartree / atom_ureg.bohr,
                        )
                    )
                except ValueError:
                    logger.warning("Failed to set gradient line")
        self._gradient = raw_gradient

    def _parse_orbitals(self, orbital: Literal["Alpha", "Beta"]):
        occ_patern = (
            r"\s+Alpha\s+occ.\s+eigenvalues\s--\s+([-0-9.\s]+)"
            if orbital == "Alpha"
            else r"\s+Beta\s+occ.\s+eigenvalues\s--\s+([-0-9.\s]+)"
        )
        virt_patern = (
            r"\s+Alpha\s+virt.\s+eigenvalues\s--\s+([-0-9.\s]+)"
            if orbital == "Alpha"
            else r"\s+Beta\s+virt.\s+eigenvalues\s--\s+([-0-9.\s]+)"
        )
        lines = self._block.splitlines()
        orbitals = []
        homo_idx = 0
        for i, line in enumerate(lines):
            match_occ = re.findall(occ_patern, line)
            match_virt = re.findall(virt_patern, line)
            if match_occ:
                orbitals.extend(
                    (
                        e * atom_ureg.hartree / atom_ureg.particle
                        for e in map(float, match_occ[0].split())
                    )
                )
                homo_idx += len(match_occ[0].split())
            if match_virt:
                orbitals.extend(
                    (
                        e * atom_ureg.hartree / atom_ureg.particle
                        for e in map(float, match_virt[0].split())
                    )
                )
        if len(orbitals) > 0:
            if orbital == "Alpha":
                self._alpha_FMO_orbits = orbitals
                self._alpha_energy["homo"] = self._alpha_FMO_orbits[homo_idx - 1]
                self._alpha_energy["lumo"] = self._alpha_FMO_orbits[homo_idx]
                self._alpha_energy["gap"] = (
                    self._alpha_energy["lumo"] - self._alpha_energy["homo"]
                )
            else:
                self._beta_FMO_orbits = orbitals
                self._beta_energy["homo"] = self._beta_FMO_orbits[homo_idx - 1]
                self._beta_energy["lumo"] = self._beta_FMO_orbits[homo_idx]
                self._beta_energy["gap"] = (
                    self._beta_energy["lumo"] - self._beta_energy["homo"]
                )

    def _parse_frequencies(self):
        lines = self._block.splitlines()
        freq_list = []
        for idx, line in enumerate(lines):
            if " and normal coordinates:" in line:
                while lines[idx + 8 + self.__n_atom] != "":
                    num = len(lines[idx + 1].split())
                    freqs = [{} for i in range(num)]
                    for i, freq in enumerate(freqs):
                        freq["freq"] = (
                            float(lines[idx + 3].split()[2 + i]) * atom_ureg.cm_1
                        )
                        freq["is imaginary"] = freq["freq"] < 0
                        freq["Reduced masses"] = (
                            float(lines[idx + 4].split()[3 + i]) * atom_ureg.amu
                        )
                        freq["force constants"] = (
                            float(lines[idx + 5].split()[3 + i])
                            * atom_ureg.mdyne
                            / atom_ureg.angstrom
                        )
                        freq["IR intensities"] = (
                            float(lines[idx + 6].split()[3 + i])
                            * atom_ureg.kmol
                            / atom_ureg.mol
                        )
                        freq["normal coordinates"] = []
                        for j in range(idx + 8, idx + 8 + self.__n_atom):
                            freq["normal coordinates"].append(
                                (
                                    float(lines[j].split()[i * 3 + 2])
                                    * atom_ureg.angstrom,
                                    float(lines[j].split()[i * 3 + 3])
                                    * atom_ureg.angstrom,
                                    float(lines[j].split()[i * 3 + 4])
                                    * atom_ureg.angstrom,
                                )
                            )
                    freq_list.extend(freqs)
                    idx += 7 + self.__n_atom
                self._frequencies = freq_list
                break

    def _parse_sum_energy(self):
        lines = self._block.splitlines()
        for line in reversed(lines):
            if "Sum of electronic and thermal Free Energies=" in line:
                self._sum_energy["thermal gibbs free energy"] = (
                    float(line.split()[-1]) * atom_ureg.hartree / atom_ureg.particle
                )
            elif "Sum of electronic and thermal Enthalpies=" in line:
                self._sum_energy["thermal enthalpy"] = (
                    float(line.split()[-1]) * atom_ureg.hartree / atom_ureg.particle
                )
            elif "Sum of electronic and thermal Energies=" in line:
                self._sum_energy["thermal energy"] = (
                    float(line.split()[-1]) * atom_ureg.hartree / atom_ureg.particle
                )
            elif "Sum of electronic and zero-point Energies=" in line:
                self._sum_energy["zero-point"] = (
                    float(line.split()[-1]) * atom_ureg.hartree / atom_ureg.particle
                )
            elif "Thermal correction to Gibbs Free Energy=" in line:
                self._correction["thermal gibbs free energy"] = (
                    float(line.split()[-1]) * atom_ureg.hartree / atom_ureg.particle
                )
            elif "Thermal correction to Enthalpy=" in line:
                self._correction["thermal enthalpy"] = (
                    float(line.split()[-1]) * atom_ureg.hartree / atom_ureg.particle
                )
            elif "Thermal correction to Energy=" in line:
                self._correction["thermal energy"] = (
                    float(line.split()[-1]) * atom_ureg.hartree / atom_ureg.particle
                )
            elif "Zero-point correction=" in line:
                self._correction["zero-point"] = (
                    float(line.split()[-2]) * atom_ureg.hartree / atom_ureg.particle
                )
                break

    def _parse_hessian(self):
        lines = self._block.splitlines()
        hess_lines = []
        append_line = False
        for line in reversed(lines):
            if r"\\@" in line or line.startswith(" @") or line.startswith(r" \@"):
                append_line = True
            if append_line:
                hess_lines.append(line.strip("\n").strip(" "))
            if "NImag" in line:
                append_line = "end"
                break
        if append_line == "end":
            hess_str = "".join(hess_lines[::-1]).split(r"\\")[-3]
            hess_val = [float(val) for val in hess_str.split(",")]
            n = self.__n_atom * 3
            if len(hess_val) != n * (n + 1) // 2:
                logger.warning(
                    f"The number of elements {len(hess_val)} in the Hessian matrix is not consistent with the number of atoms {self.__n_atom}."
                )
            else:
                self._hessian = hess_val

    def _parse_state(self):
        lines = self._block.splitlines()
        for line in reversed(lines):
            if "Normal termination" in line:
                self._state["Normal termination"] = True
            elif "RMS     Displacement" in line:
                self._state["RMS Displacement"] = True if "YES" in line else False
            elif "Maximum Displacement" in line:
                self._state["Maximum Displacement"] = True if "YES" in line else False
            elif "RMS     Force" in line:
                self._state["RMS Force"] = True if "YES" in line else False
            elif "Maximum Force" in line:
                self._state["Maximum Force"] = True if "YES" in line else False


class G16LOGParser(BaseQMFileParser):
    def __init__(self, file_path: str, charge=None, multiplicity=None, show_progress=False):
        super().__init__(file_path, show_progress)
        self.__force_charge = charge
        self.__force_multiplicity = multiplicity
        _, file_format = os.path.splitext(file_path)
        if file_format != ".log":
            raise ValueError("File format must be .log")
        self._parse()

    def _parse(self):
        """
        Parse the file.
        """
        with open(self.file_path, "r") as fr:
            lines = fr.readlines()
        fr.close()
        full_text = "".join(lines)
        charge, multi = map(
            int,
            re.findall(r"Charge\s*=\s*(\d+)\s+Multiplicity\s*=\s*(\d+)", full_text)[0],
        )
        if self.__force_charge is not None:
            charge = self.__force_charge
        if self.__force_multiplicity is not None:
            multi = self.__force_multiplicity
        pattern = r"""\s+\*+
\s+Gaussian\s+\d+\:\s+[A-Za-z0-9-.]+\s+\d+-[A-Za-z]{3}-\d{4}
\s+\d+-[A-Za-z]{3}-\d{4}\s+
\s+\*+
([a-zA-Z%0-9.=\s\_\\\/\*\+\-]+)
\s+-+
([a-zA-Z%0-9.\=\s\-\+#(),\*\/\\^\n]+)
\s+-+"""
        self._parameter_comment = "\n".join(re.findall(pattern, full_text)[0])
        n_atom = int(re.findall(r"NAtoms=\s*(\d+)", full_text)[0])
        block_starts = [
            idx for idx, line in enumerate(lines) if "Input orientation:" in line
        ] + [len(lines)]
        for idx, start in enumerate(block_starts[:-1]):
            self.append(
                G16LOGBlockParser(
                    "".join(lines[start : block_starts[idx + 1]]),
                    charge=charge,
                    multiplicity=multi,
                    n_atom=n_atom,
                    parameter_comment=self._parameter_comment,
                ),
            )
