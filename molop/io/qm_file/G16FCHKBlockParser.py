"""
Author: TMJ
Date: 2024-01-24 13:04:53
LastEditors: TMJ
LastEditTime: 2024-01-24 16:25:48
Description: 请填写简介
"""

import math
import re
from typing import Literal, List

import numpy as np

from molop.config import molopconfig
from molop.io.bases.molblock_base import QMBaseBlockParser
from molop.logger.logger import logger
from molop.unit import atom_ureg
from molop.utils.g16patterns import (
    g16fchkpatterns,
    get_solvent,
    get_solvent_model,
    parameter_comment_parser,
    shell_to_orbitals,
)


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
        self.qm_software = "Gaussian"
        self._file_path = file_path
        self._charge = charge
        self._multiplicity = multiplicity
        self.__n_atom = n_atom
        self.version = version
        self.parameter_comment = parameter_comment

        (
            self._route_params,
            self._dieze_tag,
        ) = parameter_comment_parser(self.parameter_comment)
        self.solvent_model = get_solvent_model(self.route_params)
        self.solvent = get_solvent(self.route_params)
        self._parse_functional_basis()
        self._parse_coords()
        if not self._only_extract_structure:
            self._parse()

    @property
    def route_params(self) -> dict:
        return self._route_params

    @property
    def dieze_tag(self) -> Literal["#N", "#P", "#T"]:
        return self._dieze_tag

    def _parse_functional_basis(self):
        for idx, line in enumerate(self._block.splitlines()):
            if idx == 1:
                self.functional = line.split()[1].lower()
                self.basis = line.split()[2]
                break

    def _parse_coords(self):
        atomic_numbers = self._parse_block(
            "atomic_number_start", "int_digits", end_tag="atomic_number_end"
        )
        if len(atomic_numbers) > 0:
            self._atoms = list(map(int, atomic_numbers))
        else:
            raise RuntimeError("No atomic number found in fchk file.")
        coords = self._parse_block("coords_start", "float_digits", end_tag="coords_end")
        if len(coords) > 0:
            coords = np.array(list(map(float, coords))).reshape(-1, 3)
        else:
            raise RuntimeError("No coords found in fchk file.")
        self._coords = (coords * atom_ureg.bohr).to("angstrom")

    def _parse(self):
        self._parse_energy()
        self._parse_spin()
        self._parse_state()
        self._parse_thermal_energy()
        self._parse_orbitals()
        self._parse_mulliken_charges()
        self._parse_gradient()
        self._parse_frequencies()
        # self._parse_hessian()
        self._parse_dipole()
        # self._parse_nbo()

    def _parse_energy(self):
        total_energy = g16fchkpatterns["total energy"].search(self._block)
        if total_energy:
            self._total_energy = (
                round(float(total_energy.group(1)), 6)
                * atom_ureg.hartree
                / atom_ureg.particle
            )
        scf_energy = g16fchkpatterns["scf energy"].search(self._block)
        if scf_energy:
            self.scf_energy = (
                round(float(scf_energy.group(1)), 6)
                * atom_ureg.hartree
                / atom_ureg.particle
            )
        for matches in g16fchkpatterns["mp2-4"].finditer(self._block):
            setattr(
                self,
                f"{matches.group(1).lower()}_energy",
                round(float(matches.group(2)), 6)
                * atom_ureg.hartree
                / atom_ureg.particle,
            )
        cluster_energy = g16fchkpatterns["cluster energy"].search(self._block)
        if cluster_energy:
            self.ccsd_energy = (
                round(float(cluster_energy.group(1)), 6)
                * atom_ureg.hartree
                / atom_ureg.particle
            )

    def _parse_thermal_energy(self):
        thermal_energy = g16fchkpatterns["thermal energy"].search(self._block)
        if thermal_energy:
            self.thermal_energy["U_T"] = (
                round(float(thermal_energy.group(1)), 6)
                * atom_ureg.hartree
                / atom_ureg.particle
            )
        thermal_enthalpy = g16fchkpatterns["thermal enthalpy"].search(self._block)
        if thermal_enthalpy:
            self.thermal_energy["H_T"] = (
                round(float(thermal_enthalpy.group(1)), 6)
                * atom_ureg.hartree
                / atom_ureg.particle
            )
        thermal_free_energy = g16fchkpatterns["thermal free energy"].search(self._block)
        if thermal_enthalpy:
            self.thermal_energy["G_T"] = (
                round(float(thermal_free_energy.group(1)), 6)
                * atom_ureg.hartree
                / atom_ureg.particle
            )

    def _parse_orbitals(self):
        alpha_orbitals_energy = self._parse_block(
            "alpha_start", "float_digits", "alpha_end"
        )
        beta_orbitals_energy = self._parse_block(
            "beta_start", "float_digits", "beta_end"
        )
        if len(alpha_orbitals_energy) > 0:
            elec_num = int(g16fchkpatterns["alpha_elec"].search(self._block).group(1))
            alpha_orbitals_energy = list(
                map(
                    float,
                    alpha_orbitals_energy,
                )
            )
            occ = math.ceil(elec_num / 2.0)
            self.alpha_FMO_orbits = (
                np.array(alpha_orbitals_energy, dtype=np.float32)
                * atom_ureg.hartree
                / atom_ureg.particle
            )
            self.alpha_energy["homo"] = self.alpha_FMO_orbits[occ - 1]
            self.alpha_energy["lumo"] = self.alpha_FMO_orbits[occ]
            self.alpha_energy["gap"] = (
                round((self.alpha_energy["lumo"] - self.alpha_energy["homo"]).m, 6)
                * atom_ureg.hartree
                / atom_ureg.particle
            )

        if len(beta_orbitals_energy) > 0:
            elec_num = int(g16fchkpatterns["beta_elec"].search(self._block).group(1))
            beta_orbitals_energy = list(
                map(
                    float,
                    beta_orbitals_energy,
                )
            )
            occ = math.ceil(elec_num / 2.0)
            self.beta_FMO_orbits = (
                np.array(beta_orbitals_energy, dtype=np.float32)
                * atom_ureg.hartree
                / atom_ureg.particle
            )
            self.beta_energy["homo"] = self.beta_FMO_orbits[occ - 1]
            self.beta_energy["lumo"] = self.alpha_FMO_orbits[occ]
            self.beta_energy["gap"] = (
                round((self.beta_energy["lumo"] - self.beta_energy["homo"]).m, 6)
                * atom_ureg.hartree
                / atom_ureg.particle
            )

    def _parse_mulliken_charges(self):
        mulliken_match = self._parse_block(
            "mulliken_start", "float_digits", "gradient_start"
        )
        if len(mulliken_match) > 0:
            self.mulliken_charges = list(map(float, mulliken_match))

    def _parse_gradient(self):
        gradients = self._parse_block("gradient_start", "float_digits", "gradient_end")
        if len(gradients) > 0:
            gradients = list(map(float, gradients))
        self.forces = (
            np.array(gradients).reshape(-1, 3) * atom_ureg.hartree / atom_ureg.bohr
        )

    def _parse_frequencies(self):
        freq_num = g16fchkpatterns["freq num"].search(self._block)
        if freq_num:
            num_freqs = int(freq_num.group(1))
            freqs = list(
                map(
                    float,
                    self._parse_block("vib_e2_start", "float_digits", "vib_mode_start"),
                )
            )
            freq_modes = list(
                map(
                    float,
                    self._parse_block("vib_mode_start", "float_digits", "vib_mode_end"),
                )
            )
            self.frequencies = [
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
                for idx in range(num_freqs)
            ]

    def _parse_spin(self):
        matches = g16fchkpatterns["spin"].search(self._block)
        if matches:
            self.spin_square = round(float(matches.group(1)), 2)
            self.spin_quantum_number = round(
                math.sqrt(self.spin_square + 0.25) - 0.5, 2
            )
            if molopconfig.allow_spin_change:
                self._multiplicity = int(round(2 * self.spin_quantum_number + 1, 0))

    def _parse_dipole(self):
        matches = re.search(g16fchkpatterns["dipole"], self._block)
        if matches:
            self.dipole = np.array(list(map(float, matches.groups()))) * atom_ureg.debye

    # TODO NBO section
    # route section: pop=saveNBO or pop=saveNLMO
    # ref: http://sobereva.com/134

    def _parse_state(self):
        matches = re.search(g16fchkpatterns["job status"], self._block)
        if matches:
            self.status["Job Status"] = matches.group(1) == "1"
        else:
            self.status["Job Status"] = False

    def is_error(self) -> bool:
        if self.total_energy is None:
            return True
        if "Job Status" in self.status:
            return self.status["Job Status"] == False
        else:
            return True
        
    def is_optimized(self) -> bool:
        if "opt" in self.route_params and not self.is_error():
            return True
        else:
            return False

    def _parse_block(
        self, start_tag: str, item_map: str, end_tag: str = None
    ) -> List[str]:
        start_match = g16fchkpatterns[start_tag].search(self._block)
        temp_list = []
        if start_match:
            num = int(start_match.group(1))
            if end_tag:
                end_match = g16fchkpatterns[end_tag].search(self._block)
                temp_list = [
                    matches.group(0)
                    for i, matches in enumerate(
                        g16fchkpatterns[item_map].finditer(
                            self._block[start_match.end() : end_match.start()]
                        )
                    )
                    if i < num
                ]
            else:
                temp_list = [
                    matches.group(0)
                    for i, matches in enumerate(
                        g16fchkpatterns[item_map].finditer(
                            self._block[start_match.end() :]
                        )
                    )
                    if i < num
                ]
        return temp_list
