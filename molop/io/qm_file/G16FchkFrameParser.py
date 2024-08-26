"""
Author: TMJ
Date: 2024-06-21 11:03:25
LastEditors: TMJ
LastEditTime: 2024-06-24 20:45:32
Description: 请填写简介
"""

import math
import re
from typing import List, Literal

import numpy as np
from pydantic import Field, computed_field

from molop.config import molopconfig
from molop.io.bases.BaseMolFrameParser import BaseQMMolFrameParser
from molop.io.bases.DataClasses import (
    ChargeSpinPopulations,
    Energies,
    MolecularOrbitals,
    Polarizability,
    TotalSpin,
    Vibrations,
)
from molop.unit import atom_ureg
from molop.utils.g16patterns import (
    g16fchkpatterns,
    get_solvent,
    get_solvent_model,
    link0_parser,
    parameter_comment_parser,
)


class G16FchkFrameParser(BaseQMMolFrameParser):
    _frame_type: str = "G16 Fchk"
    qm_software: str = Field(default="Gaussian")
    n_atom: int = Field(default=0, exclude=True, repr=False)
    only_extract_structure: bool = Field(default=False, exclude=True, repr=False)

    @property
    def link(self) -> str:
        """
        Get the link.
        """
        return self.keywords.split("#")[0]

    @property
    def route(self):
        return "#" + self.keywords.split("#")[1].replace("\n", " ")

    @property
    def link0(self) -> dict:
        """
        Get the link0.
        """
        return link0_parser(self.link)

    @property
    def route_params(self) -> dict:
        """
        Get the route parameters.
        """
        return parameter_comment_parser(self.route)[0]

    @property
    def dieze_tag(self) -> Literal["#N", "#P", "#T"]:
        return parameter_comment_parser(self.route)[1]

    def _parse(self):
        self.__block = self.frame_content
        self._parse_functional_basis()
        self.solvent_model = get_solvent_model(self.route_params)
        self.solvent = get_solvent(self.route_params)
        self._parse_coords()
        if self.only_extract_structure:
            return
        self._parse_energy()
        self._parse_spin()
        self._parse_gradient()
        self._parse_orbitals()
        self.charge_spin_populations = ChargeSpinPopulations(
            mulliken_charges=self._parse_mulliken_charges()
        )
        self._parse_polarizability()
        self._parse_vibrations()
        self._parse_state()

    def _parse_functional_basis(self):
        for idx, line in enumerate(self.__block.splitlines()):
            if idx == 1:
                self.functional = line.split()[1].lower()
                self.basis = line.split()[2].lower()
                break

    def _parse_coords(self):
        atomic_numbers = self.__parse_block__(
            "atomic_number_start", "int_digits", end_tag="atomic_number_end"
        )
        if len(atomic_numbers) > 0:
            self.atoms = list(map(int, atomic_numbers))
        else:
            raise RuntimeError("No atomic number found in fchk file.")
        coords = self.__parse_block__(
            "coords_start", "float_digits", end_tag="coords_end"
        )
        if len(coords) > 0:
            coords = np.array(list(map(float, coords))).reshape(-1, 3)
        else:
            raise RuntimeError("No coords found in fchk file.")
        self.coords = (coords * atom_ureg.bohr).to("angstrom")
        self.standard_coords = self.coords

    def _parse_energy(self) -> Energies:
        energies = {}
        total_energy = g16fchkpatterns["total energy"].search(self.__block)
        if total_energy:
            energies["total_energy"] = (
                float(total_energy.group(1)) * atom_ureg.hartree / atom_ureg.particle
            )
        scf_energy = g16fchkpatterns["scf energy"].search(self.__block)
        if scf_energy:
            energies["scf_energy"] = (
                float(scf_energy.group(1)) * atom_ureg.hartree / atom_ureg.particle
            )
        for matches in g16fchkpatterns["mp2-4"].finditer(self.__block):
            energies[f"{matches.group(1).lower()}_energy"] = (
                float(matches.group(2)) * atom_ureg.hartree / atom_ureg.particle
            )
        cluster_energy = g16fchkpatterns["cluster energy"].search(self.__block)
        if cluster_energy:
            energies["ccsd_energy"] = (
                float(cluster_energy.group(1)) * atom_ureg.hartree / atom_ureg.particle
            )
        if len(energies):
            self.energies = Energies.model_validate(energies)

    def _parse_spin(self) -> TotalSpin:
        matches = g16fchkpatterns["spin"].search(self.__block)
        spins = {}
        if matches:
            spins["spin_square"] = float(matches.group(1))
            spins["spin_quantum_number"] = math.sqrt(spins["spin_square"] + 0.25) - 0.5
            if molopconfig.allow_spin_change:
                self.multiplicity = int(round(2 * spins["spin_quantum_number"] + 1, 0))
            self.total_spin = TotalSpin.model_validate(spins)

    def _parse_mulliken_charges(self):
        mulliken_match = self.__parse_block__(
            "mulliken_start", "float_digits", "gradient_start"
        )
        if len(mulliken_match) > 0:
            return list(map(float, mulliken_match))
        return []

    def _parse_orbitals(self):
        alpha_elec_num = int(
            g16fchkpatterns["alpha_elec"].search(self.__block).group(1)
        )
        beta_elec_num = int(g16fchkpatterns["beta_elec"].search(self.__block).group(1))
        alpha_orbitals_energy = self.__parse_block__(
            "alpha_start", "float_digits", "alpha_end"
        )
        beta_orbitals_energy = self.__parse_block__(
            "beta_start", "float_digits", "beta_end"
        )
        if len(alpha_orbitals_energy) > 0:
            alpha_orbitals_energy = list(
                map(
                    float,
                    alpha_orbitals_energy,
                )
            )
            alpha_occ = [
                True if idx < alpha_elec_num else False
                for idx, orbital in enumerate(alpha_orbitals_energy)
            ]
        if len(beta_orbitals_energy) > 0:
            beta_orbitals_energy = list(
                map(
                    float,
                    beta_orbitals_energy,
                )
            )
            beta_occ = [
                True if idx < beta_elec_num else False
                for idx, orbital in enumerate(beta_orbitals_energy)
            ]
        else:
            beta_orbitals_energy = alpha_orbitals_energy
            beta_occ = alpha_occ
        if len(alpha_orbitals_energy):
            self.molecular_orbitals = MolecularOrbitals(
                alpha_energies=np.array(alpha_orbitals_energy)
                * atom_ureg.hartree
                / atom_ureg.particle,
                beta_energies=np.array(beta_orbitals_energy)
                * atom_ureg.hartree
                / atom_ureg.particle,
                alpha_occupancies=alpha_occ,
                beta_occupancies=beta_occ,
            )

    def _parse_gradient(self):
        gradients = self.__parse_block__(
            "gradient_start", "float_digits", "gradient_end"
        )
        if len(gradients) > 0:
            gradients = list(map(float, gradients))
        self.forces = (
            -1 * np.array(gradients).reshape(-1, 3) * atom_ureg.hartree / atom_ureg.bohr
        )

    def _parse_vibrations(self):
        freq_num = g16fchkpatterns["freq num"].search(self.__block)
        if freq_num:
            num_freqs = int(freq_num.group(1))
            freqs = list(
                map(
                    float,
                    self.__parse_block__(
                        "vib_e2_start", "float_digits", "vib_mode_start"
                    ),
                )
            )
            freq_modes = list(
                map(
                    float,
                    self.__parse_block__(
                        "vib_mode_start", "float_digits", "vib_mode_end"
                    ),
                )
            )
            self.vibrations = Vibrations(
                frequencies=np.array([freqs[idx] for idx in range(num_freqs)])
                * atom_ureg.cm_1,
                reduced_masses=np.array(
                    [freqs[idx + num_freqs] for idx in range(num_freqs)]
                )
                * atom_ureg.amu,
                force_constants=np.array(
                    [freqs[idx + num_freqs * 2] for idx in range(num_freqs)]
                )
                * atom_ureg.mdyne
                / atom_ureg.angstrom,
                IR_intensities=np.array(
                    [freqs[idx + num_freqs * 3] for idx in range(num_freqs)]
                )
                * atom_ureg.kmol
                / atom_ureg.mol,
                vibration_mode=[
                    np.array(
                        freq_modes[idx * 3 * self.n_atom : (idx + 1) * 3 * self.n_atom]
                    ).reshape(-1, 3)
                    * atom_ureg.angstrom
                    for idx in range(num_freqs)
                ],
            )

    def _parse_polarizability(self):
        polar = {}
        dipole = self.__parse_block__("dipole_start", "float_digits")
        if len(dipole) > 0:
            polar["dipole"] = np.array(list(map(float, dipole))) * atom_ureg.debye
        polarizability = self.__parse_block__("polarizability_start", "float_digits")
        if len(polarizability) > 0:
            polar["polarizability"] = list(map(float, polarizability))
        quadrupole = self.__parse_block__("quadrupole_start", "float_digits")
        if len(quadrupole) > 0:
            polar["quadrupole"] = (
                np.array(list(map(float, quadrupole)))
                * atom_ureg.debye
                * atom_ureg.angstrom
            )
        if len(polar):
            self.polarizability = Polarizability.model_validate(polar)

    def _parse_state(self):
        if matches := re.search(g16fchkpatterns["job status"], self.__block):
            self.status.normal_terminated = matches.group(1) == "1"
        else:
            self.status.normal_terminated = False
        self.status.scf_converged = self.energies.total_energy is not None
        if "opt" in self.route_params and not self.is_error:
            self.geometry_optimization_status.geometry_optimized = True

    @computed_field
    @property
    def is_error(self) -> bool:
        if self.energies.total_energy is None:
            return True
        if not self.status.normal_terminated:
            return True
        return False

    @property
    def is_optimized(self) -> bool:
        return self.geometry_optimization_status.geometry_optimized

    def __parse_block__(
        self, start_tag: str, item_map: str, end_tag: str = None
    ) -> List[str]:
        start_match = g16fchkpatterns[start_tag].search(self.__block)
        temp_list = []
        if start_match:
            num = int(start_match.group(1))
            if end_tag:
                end_match = g16fchkpatterns[end_tag].search(self.__block)
                temp_list = [
                    matches.group(0)
                    for i, matches in enumerate(
                        g16fchkpatterns[item_map].finditer(
                            self.__block[start_match.end() : end_match.start()]
                        )
                    )
                    if i < num
                ]
            else:
                temp_list = [
                    matches.group(0)
                    for i, matches in enumerate(
                        g16fchkpatterns[item_map].finditer(
                            self.__block[start_match.end() :]
                        )
                    )
                    if i < num
                ]
        return temp_list
