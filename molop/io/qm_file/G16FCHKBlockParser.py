"""
Author: TMJ
Date: 2024-01-24 13:04:53
LastEditors: TMJ
LastEditTime: 2024-01-24 16:25:48
Description: 请填写简介
"""
import math
import re
from typing import Literal

import numpy as np

from molop.io.bases.molblock_base import QMBaseBlockParser
from molop.logger.logger import logger
from molop.unit import atom_ureg
from molop.utils import g16fchkpatterns, parameter_comment_parser


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

        (
            _,
            self._route_params,
            self._dieze_tag,
            self._functional,
            self._basis_set,
        ) = parameter_comment_parser("\n" + self._parameter_comment)
        self._parse_coords()
        if not self._only_extract_structure:
            self._parse()

    @property
    def route_params(self) -> dict:
        return self._route_params

    @property
    def dieze_tag(self) -> Literal["#N", "#P", "#T"]:
        return self._dieze_tag

    @property
    def functional(self) -> str:
        return self._functional

    @property
    def basis_set(self) -> str:
        return self._basis_set

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
        temp_coords = []
        for i in range(0, len(coords), 3):
            temp_coords.append((coords[i], coords[i + 1], coords[i + 2]))
        self._coords = (np.array(temp_coords) * atom_ureg.bohr).to("angstrom")

    def _parse(self):
        self._parse_energy()
        self._parse_partial_charges()
        self._parse_gradient()
        self._parse_orbitals()
        self._parse_frequencies()
        self._parse_spin()
        # self._parse_hessian()
        self._parse_dipole()
        self._parse_state()
        # self._parse_nbo()

    def _parse_energy(self):
        total_energy = re.search(g16fchkpatterns["total energy"], self._block)
        if total_energy:
            self._energy = (
                round(float(total_energy.group(1)), 6)
                * atom_ureg.hartree
                / atom_ureg.particle
            )
        else:
            scf_energy = re.search(g16fchkpatterns["scf energy"], self._block)
            if scf_energy:
                self._energy = (
                    round(float(scf_energy.group(1)), 6)
                    * atom_ureg.hartree
                    / atom_ureg.particle
                )
        thermal_energy = re.search(g16fchkpatterns["thermal energy"], self._block)
        if thermal_energy:
            self._sum_energy["E gas"] = (
                round(float(thermal_energy.group(1)), 6)
                * atom_ureg.hartree
                / atom_ureg.particle
            )
        thermal_enthalpy = re.search(g16fchkpatterns["thermal enthalpy"], self._block)
        if thermal_enthalpy:
            self._sum_energy["H gas"] = (
                round(float(thermal_enthalpy.group(1)), 6)
                * atom_ureg.hartree
                / atom_ureg.particle
            )
        thermal_free_energy = re.search(
            g16fchkpatterns["thermal free energy"], self._block
        )
        if thermal_enthalpy:
            self._sum_energy["G gas"] = (
                round(float(thermal_free_energy.group(1)), 6)
                * atom_ureg.hartree
                / atom_ureg.particle
            )

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
        self._gradients = (
            np.array(raw_gradient).reshape(-1, 3) * atom_ureg.hartree / atom_ureg.bohr
        )

    def _parse_orbitals(self):
        lines = self._block.splitlines()
        orbitals_energy = []
        for i, line in enumerate(lines):
            matches = re.search(g16fchkpatterns["orbital"], line)
            if matches:
                for j in range(i + 1, i + 1 + math.ceil(int(line.split()[-1]) / 5.0)):
                    orbitals_energy.extend(list(map(float, lines[j].split())))
                break
        occ = math.ceil(self.total_electrons / 2.0)
        if len(orbitals_energy) == 0:
            return
        if matches.group(1) == "Alpha":
            self._alpha_FMO_orbits = (
                np.array(orbitals_energy, dtype=np.float32)
                * atom_ureg.hartree
                / atom_ureg.particle
            )
            self._alpha_energy["homo"] = self._alpha_FMO_orbits[occ - 1]
            self._alpha_energy["lumo"] = self._alpha_FMO_orbits[occ]
            self._alpha_energy["gap"] = (
                round((self._alpha_energy["lumo"] - self._alpha_energy["homo"]).m, 6)
                * atom_ureg.hartree
                / atom_ureg.particle
            )
        else:
            self._beta_FMO_orbits = (
                np.array(orbitals_energy, dtype=np.float32)
                * atom_ureg.hartree
                / atom_ureg.particle
            )
            self._beta_energy["homo"] = self._beta_FMO_orbits[occ - 1]
            self._beta_energy["lumo"] = self._beta_FMO_orbits[occ]
            self._beta_energy["gap"] = (
                round((self._beta_energy["lumo"] - self._beta_energy["homo"]).m, 6)
                * atom_ureg.hartree
                / atom_ureg.particle
            )

    def _parse_frequencies(self):
        freq_num = re.search(g16fchkpatterns["freq num"], self._block)
        if freq_num:
            num_freqs = int(freq_num.group(1))
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
                        "reduced masses": freqs[idx + num_freqs] * atom_ureg.amu,
                        "force constants": freqs[idx + num_freqs * 2]
                        * atom_ureg.mdyne
                        / atom_ureg.angstrom,
                        "IR intensities": freqs[idx + num_freqs * 3]
                        * atom_ureg.kmol
                        / atom_ureg.mol,
                        "normal coordinates": np.array(
                            freq_modes[
                                idx * 3 * self.__n_atom : (idx + 1) * 3 * self.__n_atom
                            ]
                        ).reshape(-1, 3)
                        * atom_ureg.angstrom,
                    }
                )

    def _parse_spin(self):
        matches = re.search(g16fchkpatterns["spin"], self._block)
        if matches:
            self._spin_eigenvalue = round(float(matches.group(1)), 2)
            self._spin_multiplicity = round(
                math.sqrt(self._spin_eigenvalue + 0.25) - 0.5, 2
            )

    def _parse_dipole(self):
        matches = re.search(g16fchkpatterns["dipole"], self._block)
        if matches:
            self._dipole = np.array(
                list(map(float, matches.groups()))
            ) * atom_ureg.debye

    #TODO NBO section
    # route section: pop=saveNBO or pop=saveNLMO
    # ref: http://sobereva.com/134

    def _parse_state(self):
        matches = re.search(g16fchkpatterns["job status"], self._block)
        if matches:
            self._state["Job Status"] = matches.group(1) == "1"
        else:
            self._state["Job Status"] = False

    def is_error(self) -> bool:
        if "Job Status" in self._state:
            return self._state["Job Status"] == False
        else:
            return True
