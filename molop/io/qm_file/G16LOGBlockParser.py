import os
import re
from typing import Literal

from molop.io.bases.molblock_base import QMBaseBlockParser
from molop.logger.logger import logger
from molop.unit import atom_ureg


class G16LOGBlockParser(QMBaseBlockParser):
    """
    Parser for G16 log Blocks.
    """

    _block_type = "G16 LOG"

    def __init__(
        self,
        block: str,
        charge=0,
        multiplicity=1,
        n_atom=1,
        file_path="",
        version=None,
        parameter_comment=None,
        only_extract_structure=False,
    ):
        super().__init__(block, only_extract_structure)
        self._file_path = file_path
        self._charge = charge
        self._multiplicity = multiplicity
        self.__n_atom = n_atom
        self._version = version
        self._parameter_comment = parameter_comment
        self._parse_coords()
        if not self._only_extract_structure:
            self._parse()

    def _parse(self):
        self._parse_energy()
        self._parse_partial_charges_and_spins()
        self._parse_gradient()
        self._parse_orbitals("Alpha")
        self._parse_orbitals("Beta")
        self._parse_frequencies()
        self._parse_sum_energy()
        # self._parse_hessian()
        self._parse_state()
        # self._parse_nbo()

    def _parse_coords(self):
        lines = self._block.splitlines()
        for i, line in enumerate(lines):
            if "Input orientation:" in line:
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
                break
            if "Standard orientation:" in line:
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
                break

    def _parse_energy(self):
        lines = self._block.splitlines()
        for line in reversed(lines):
            if "SCF Done" in line or "E(CIS)" in line:
                self._energy = (
                    round(float(line.split()[4]), 6)
                    * atom_ureg.hartree
                    / atom_ureg.particle
                )
                self._state["SCF Done"] = True
                return
            if "E(CORR)" in line or "E(CI)" in line:
                self._energy = (
                    round(float(line.split()[3]), 6)
                    * atom_ureg.hartree
                    / atom_ureg.particle
                )
                return
            if "E(CIS(D))" in line:
                self._energy = (
                    round(float(line.split()[5]), 6)
                    * atom_ureg.hartree
                    / atom_ureg.particle
                )
                return
            if line.startswith(" Energy=") and "NIter=" in line:
                self._energy = (
                    round(float(line.split()[1]), 6)
                    * atom_ureg.hartree
                    / atom_ureg.particle
                )
                return
        try:
            self._energy = (
                round(
                    float(
                        "".join(
                            re.findall(r"HF=([\-0-9\s+\.]+)\\", self._block)[0].split(
                                "\n "
                            )
                        )
                    ),
                    6,
                )
                * atom_ureg.hartree
                / atom_ureg.particle
            )
        except:
            self._state["SCF Done"] = False

    def _parse_partial_charges_and_spins(self):
        lines = self._block.splitlines()
        charges_section = False
        charges = []
        spins = []
        for line in reversed(lines):
            if len(charges) == self.__n_atom:
                self._partial_charges = list(reversed(charges))
                self._spin_densities = list(reversed(spins))
                break
            if charges_section:
                charges.append(float(line.split()[2]))
                if len(line.split()) == 4:
                    spins.append(float(line.split()[3]))
            if "sum of mulliken charges" in line.lower():
                charges_section = True
        s = re.findall(
            r"<S\*\*2>=\s+([0-9\.\-]+)\s+S=\s+([0-9\.\-]+)",
            self._block,
        )
        if s:
            self._spin_multiplicity = float(s[-1][1])
            self._spin_eigenvalue = float(s[-1][0])

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
                    logger.warning(
                        f"Failed to set gradient line {i} in {self._file_path}"
                    )
        self._gradients = raw_gradient

    def _parse_orbitals(self, orbital: Literal["Alpha", "Beta"]):
        occ_patern = (
            r"\s+Alpha\s+occ.\s+eigenvalues\s--\s([-0-9.\s]+)"
            if orbital == "Alpha"
            else r"\s+Beta\s+occ.\s+eigenvalues\s--\s([-0-9.\s]+)"
        )
        virt_patern = (
            r"\s+Alpha\s+virt.\s+eigenvalues\s--\s([-0-9.\s]+)"
            if orbital == "Alpha"
            else r"\s+Beta\s+virt.\s+eigenvalues\s--\s([-0-9.\s]+)"
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
                        round(e, 6) * atom_ureg.hartree / atom_ureg.particle
                        for e in map(
                            float,
                            [
                                match_occ[0][j : j + 10]
                                for j in range(0, len(match_occ[0]), 10)
                            ],
                        )
                    )
                )
                homo_idx += len(match_occ[0].split())
            if match_virt:
                orbitals.extend(
                    (
                        round(e, 6) * atom_ureg.hartree / atom_ureg.particle
                        for e in map(
                            float,
                            [
                                match_virt[0][j : j + 10]
                                for j in range(0, len(match_virt[0]), 10)
                            ],
                        )
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
                        freq["reduced masses"] = (
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
        def get_energy(pattern: str):
            match = re.findall(
                pattern,
                self._block,
            )
            if match:
                return (
                    round(float(match[0]), 6) * atom_ureg.hartree / atom_ureg.particle
                )

        self._sum_energy["G gas"] = get_energy(
            r"Sum of electronic and thermal Free Energies=\s+([\-0-9.]+)"
        )
        self._sum_energy["H gas"] = get_energy(
            r"Sum of electronic and thermal Enthalpies=\s+([\-0-9.]+)"
        )
        self._sum_energy["E gas"] = get_energy(
            r"Sum of electronic and thermal Energies=\s+([\-0-9.]+)"
        )
        self._sum_energy["zero-point gas"] = get_energy(
            r"Sum of electronic and zero-point Energies=\s+([\-0-9.]+)"
        )
        self._sum_energy["TCG"] = get_energy(
            r"Thermal correction to Gibbs Free Energy=\s+([\-0-9.]+)"
        )
        self._sum_energy["TCH"] = get_energy(
            r"Thermal correction to Enthalpy=\s+([\-0-9.]+)"
        )
        self._sum_energy["TCE"] = get_energy(
            r"Thermal correction to Energy=\s+([\-0-9.]+)"
        )
        self._sum_energy["zero-point correction"] = get_energy(
            r"Zero-point correction=\s+([\-0-9.]+)"
        )

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
                    f"The number of elements {len(hess_val)} in the Hessian matrix is not consistent with the number of atoms {self.__n_atom} in {self._file_path}"
                )
            else:
                self._hessian = hess_val

    def _parse_state(self):
        lines = self._block.splitlines()
        self._state["Normal termination"] = False
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
                return

    def _parse_nbo(self):
        lines = self._block.splitlines()
        for idx, line in enumerate(lines):
            if "Summary of Natural Population Analysis" in line:
                for i in range(idx + 6, idx + 6 + self.__n_atom):
                    _, _, natural_charge, core, valence, rydberg, _ = lines[i].split()
                    self._nbo_analysis.append(
                        {
                            "natural_charge": float(natural_charge),
                            "core": float(core),
                            "valence": float(valence),
                            "rydberg": float(rydberg),
                        }
                    )
                break
