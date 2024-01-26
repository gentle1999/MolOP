"""
Author: TMJ
Date: 2024-01-24 13:04:53
LastEditors: TMJ
LastEditTime: 2024-01-24 16:25:48
Description: 请填写简介
"""
import itertools
import math
import os
import re
from typing import Literal

from molop.io.bases.molblock_base import QMBaseBlockParser
from molop.logger.logger import logger
from molop.unit import atom_ureg


class G16FCHKBlockParser(QMBaseBlockParser):
    """
    Parser for G16 fchk Blocks.
    """

    _block_type = "G16 FCHK"

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

    def _parse_coords(self):
        lines = self._block.splitlines()
        coords = []
        for i, line in enumerate(lines):
            if "Atomic numbers" in line:
                for j in range(i + 1, i + 1 + math.ceil(self.__n_atom / 6.0)):
                    self._atoms.extend(list(map(int, lines[j].split())))
                assert (
                    len(self.atoms) == self.__n_atom
                ), "Number of atoms is not consistent."
            if "Current cartesian coordinates" in line:
                for j in range(i + 1, i + 1 + math.ceil(self.__n_atom * 3 / 5.0)):
                    coords.extend(list(map(float, lines[j].split())))
                assert (
                    len(coords) == self.__n_atom * 3
                ), "Number of coordinates is not consistent."
                break
        for i in range(0, len(coords), 3):
            self._coords.append(
                (
                    (coords[i] * atom_ureg.bohr).to("angstrom"),
                    (coords[i + 1] * atom_ureg.bohr).to("angstrom"),
                    (coords[i + 2] * atom_ureg.bohr).to("angstrom"),
                )
            )

    def _parse(self):
        self._parse_energy()
        self._parse_partial_charges()
        self._parse_gradient()
        self._parse_orbitals("Alpha")
        self._parse_orbitals("Beta")
        self._parse_frequencies()
        # self._parse_hessian()
        self._parse_state()
        # self._parse_nbo()

    def _parse_energy(self):
        try:
            self._energy = (
                round(
                    float(
                        re.findall(
                            "Total Energy\s+[A-Z]+\s+([\-\+0-9\.E]+)", self._block
                        )[0]
                    ),
                    6,
                )
                * atom_ureg.hartree
                / atom_ureg.particle
            )
        except:
            self._energy = (
                round(
                    float(
                        re.findall(
                            "SCF Energy\s+[A-Z]+\s+([\-\+0-9\.E]+)", self._block
                        )[0]
                    ),
                    6,
                )
                * atom_ureg.hartree
                / atom_ureg.particle
            )
        try:
            self._sum_energy["thermal energy"] = (
                round(
                    float(
                        re.findall(
                            "Thermal Energy\s+[A-Z]+\s+([\-\+0-9\.E]+)", self._block
                        )[0]
                    ),
                    6,
                )
                * atom_ureg.hartree
                / atom_ureg.particle
            )
        except:
            pass
        try:
            self._sum_energy["thermal enthalpy"] = (
                round(
                    float(
                        re.findall(
                            "Thermal Enthalpy\s+[A-Z]+\s+([\-\+0-9\.E]+)", self._block
                        )[0]
                    ),
                    6,
                )
                * atom_ureg.hartree
                / atom_ureg.particle
            )
        except:
            pass
        try:
            self._sum_energy["thermal gibbs free energy"] = (
                round(
                    float(
                        re.findall(
                            "Thermal Free Energy\s+[A-Z]+\s+([\-\+0-9\.E]+)",
                            self._block,
                        )[0]
                    ),
                    6,
                )
                * atom_ureg.hartree
                / atom_ureg.particle
            )
        except:
            pass

    def _parse_partial_charges(self):
        lines = self._block.splitlines()
        for i, line in enumerate(lines):
            if "Mulliken Charges" in line:
                for j in range(i + 1, i + 1 + math.ceil(self.__n_atom / 5.0)):
                    self._partial_charges.extend(list(map(float, lines[j].split())))
                assert (
                    len(self._partial_charges) == self.__n_atom
                ), "Number of charges is not consistent."
                break

    def _parse_gradient(self):
        lines = self._block.splitlines()
        raw_gradient = []
        for i, line in enumerate(lines):
            if "Cartesian Gradient" in line:
                for j in range(i + 1, i + 1 + math.ceil(self.__n_atom * 3 / 5.0)):
                    raw_gradient.extend(list(map(float, lines[j].split())))
                assert (
                    len(raw_gradient) == self.__n_atom * 3
                ), "Number of gradient is not consistent."
                break
        for i in range(0, len(raw_gradient), 3):
            self._gradients.append(
                (
                    raw_gradient[i] * atom_ureg.hartree / atom_ureg.bohr,
                    raw_gradient[i + 1] * atom_ureg.hartree / atom_ureg.bohr,
                    raw_gradient[i + 2] * atom_ureg.hartree / atom_ureg.bohr,
                )
            )

    def _parse_orbitals(self, orbitals: Literal["Alpha", "Beta"]):
        lines = self._block.splitlines()
        orbitals_energy = []
        for i, line in enumerate(lines):
            if f"{orbitals} Orbital Energies" in line:
                for j in range(i + 1, i + 1 + math.ceil(int(line.split()[-1]) / 5.0)):
                    orbitals_energy.extend(list(map(float, lines[j].split())))
                break
        occ = math.ceil(self.total_electrons / 2.0)
        if len(orbitals_energy) == 0:
            return
        if orbitals == "Alpha":
            self._alpha_FMO_orbits = [
                (round(energy, 6) * atom_ureg.hartree / atom_ureg.particle)
                for energy in orbitals_energy
            ]
            self._alpha_energy["homo"] = self._alpha_FMO_orbits[occ - 1]
            self._alpha_energy["lumo"] = self._alpha_FMO_orbits[occ]
            self._alpha_energy["gap"] = (
                self._alpha_energy["lumo"] - self._alpha_energy["homo"]
            )
        else:
            self._beta_FMO_orbits = [
                (round(energy, 6) * atom_ureg.hartree / atom_ureg.particle)
                for energy in orbitals_energy
            ]
            self._beta_energy["homo"] = self._beta_FMO_orbits[occ - 1]
            self._beta_energy["lumo"] = self._beta_FMO_orbits[occ]
            self._beta_energy["gap"] = (
                self._beta_energy["lumo"] - self._beta_energy["homo"]
            )

    def _parse_frequencies(self):
        try:
            num_freqs = int(
                re.findall(
                    "Number of Normal Modes\s+[A-Z]+\s+([\-\+0-9\.E]+)", self._block
                )[0]
            )
            freqs = []
            freq_modes = []
            lines = self._block.splitlines()
            for idx, line in enumerate(lines):
                if "Vib-E2" in line:
                    for j in range(
                        idx + 1, idx + 1 + math.ceil(int(line.split()[-1]) / 5.0)
                    ):
                        freqs.extend(list(map(float, lines[j].split())))
                if "Vib-Modes" in line:
                    for j in range(
                        idx + 1, idx + 1 + math.ceil(int(line.split()[-1]) / 5.0)
                    ):
                        freq_modes.extend(list(map(float, lines[j].split())))
            for idx in range(num_freqs):
                self._frequencies.append(
                    {
                        "is imaginary": freqs[idx] < 0,
                        "freq": freqs[idx] * atom_ureg.cm_1,
                        "Reduced masses": freqs[idx + num_freqs] * atom_ureg.amu,
                        "force constants": freqs[idx + num_freqs * 2]
                        * atom_ureg.mdyne
                        / atom_ureg.angstrom,
                        "IR intensities": freqs[idx + num_freqs * 3]
                        * atom_ureg.kmol
                        / atom_ureg.mol,
                        "normal coordinates": [
                            (
                                x * atom_ureg.angstrom,
                                y * atom_ureg.angstrom,
                                z * atom_ureg.angstrom,
                            )
                            for x, y, z in [
                                freq_modes[i : i + 3]
                                for i in range(
                                    idx * 3 * self.__n_atom,
                                    (idx + 1) * 3 * self.__n_atom,
                                    3,
                                )
                            ]
                        ],
                    }
                )
        except:
            logger.info(f"Frequencies not found in {self._file_path}")

    def _parse_state(self):
        try:
            self._state["Job Status"] = re.findall(
                "Job Status\s+[A-Z]+\s+([\-\+0-9\.E]+)", self._block
            )[0]
        except:
            self._state["Job Status"] = False
