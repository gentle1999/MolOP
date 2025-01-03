"""
Author: TMJ
Date: 2024-06-21 11:02:52
LastEditors: TMJ
LastEditTime: 2024-06-21 11:39:32
Description: 请填写简介
"""

import math
from itertools import chain
from typing import Literal, Tuple

import numpy as np
from pint.facets.plain import PlainQuantity
from pydantic import Field, computed_field

from molop.config import molopconfig
from molop.io.bases.BaseMolFrameParser import BaseQMMolFrameParser
from molop.io.bases.DataClasses import (
    BondOrders,
    ChargeSpinPopulations,
    Energies,
    GeometryOptimizationStatus,
    MolecularOrbitals,
    Polarizability,
    ThermalEnergies,
    TotalSpin,
    Vibrations,
)
from molop.logger.logger import moloplogger
from molop.unit import atom_ureg
from molop.utils.functions import fill_symmetric_matrix
from molop.utils.g16patterns import (
    g16logpatterns,
    link0_parser,
    parameter_comment_parser,
    semi_empirical_methods,
)


class G16LogFrameParser(BaseQMMolFrameParser):
    _frame_type: str = "G16 LOG"
    qm_software: str = Field(default="Gaussian")
    n_atom: int = Field(default=0, exclude=True, repr=False)
    only_extract_structure: bool = Field(default=False, exclude=True, repr=False)

    @property
    def link(self) -> str:
        """
        Get the link.
        """
        return self.keywords.split("\n\n")[0].replace("\n", " ")

    @property
    def route(self):
        return self.keywords.split("\n\n")[1].replace("\n", " ")

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

    @property
    def _tail(self) -> str:
        if (start := g16logpatterns["tail_start"].search(self.__block)) and (
            end := g16logpatterns["tail_end"].search(self.__block)
        ):
            tail = self.__block[start.end() : end.start()]
            tail = tail.replace("\n ", "").replace("|", "\\")
            return tail
        return ""

    def _parse(self):
        self.__block = self.frame_content
        self._parse_time()
        self._parse_coords()
        if self.only_extract_structure:
            return
        self._parse_rotation_consts()
        self._parse_functional_basis()
        self._parse_solvent()
        energies = self._parse_energy()
        self._parse_total_spin()
        self._parse_orbitals()
        mulliken_charges = self._parse_mulliken_charges()
        mulliken_spins = self._parse_mulliken_spin_density()
        apt_charges = self._parse_apt_charges()
        lowdin_charges = self._parse_lowdin_charges()
        self.polarizability = self._parse_polarizability()
        hirshfeld_charges, hirshfeld_spins, hirshfeld_q_cm5 = (
            self._parse_hirshfeld_spin_charges()
        )
        wiberg_bond_order = self._parse_bond_order("wiberg")
        atom_atom_overlap_bond_order = self._parse_bond_order("atom_atom_overlap")
        mo_bond_order = self._parse_bond_order("mo")
        self.charge_spin_populations = ChargeSpinPopulations(
            mulliken_charges=mulliken_charges,
            mulliken_spins=mulliken_spins,
            apt_charges=apt_charges,
            lowdin_charges=lowdin_charges,
            hirshfeld_charges=hirshfeld_charges,
            hirshfeld_spins=hirshfeld_spins,
            hirshfeld_q_cm5=hirshfeld_q_cm5,
            npa_charges=self._parse_npa_charges(),
        )
        self.bond_orders = BondOrders(
            wiberg_bond_order=wiberg_bond_order,
            mo_bond_order=mo_bond_order,
            atom_atom_overlap_bond_order=atom_atom_overlap_bond_order,
            nbo_bond_order=self._parse_nbo_bond_order(),
        )
        self._parse_vibrations()
        self._parse_temperature()
        thermal_energies = self._parse_thermal_energies()
        self._parse_forces()
        self._parse_status()

        self.energies = Energies.model_validate(
            {
                **energies.model_dump(
                    exclude_unset=True, exclude_none=True, exclude_defaults=True
                ),
                **self._parse_tail_energies().model_dump(
                    exclude_unset=True, exclude_none=True, exclude_defaults=True
                ),
            }
        )
        self.thermal_energies = ThermalEnergies.model_validate(
            {
                **thermal_energies.model_dump(
                    exclude_unset=True, exclude_none=True, exclude_defaults=True
                ),
                **self._parse_tail_thermal_energies().model_dump(
                    exclude_unset=True, exclude_none=True, exclude_defaults=True
                ),
            }
        )
        self._add_em()

    def _parse_time(self):
        if time_match := g16logpatterns["time"].findall(self.__block):
            self.running_time = (
                sum(
                    float(cpu_time) + float(elapsed_time)
                    for mem, cpu_time, elapsed_time in time_match
                )
                * atom_ureg.second
            )

    def _parse_coords(self):
        if corrds_end := g16logpatterns["coords_end"].search(self.__block):
            block = self.__block[: corrds_end.end()]
            if input_coords := g16logpatterns["input_coords_start"].search(block):
                coords = g16logpatterns["coords"].findall(
                    block,
                )[: self.n_atom]
                atoms = []
                temp_coords = []
                for atom_num, x, y, z in coords:
                    atoms.append(int(atom_num))
                    temp_coords.append((float(x), float(y), float(z)))
                self.coords = np.array(temp_coords) * atom_ureg.angstrom
                self.atoms = atoms
            if standard_coords := g16logpatterns["standard_coords_start"].search(block):
                coords = g16logpatterns["coords"].findall(block)[-self.n_atom :]
                atoms = []
                temp_coords = []
                for atom_num, x, y, z in coords:
                    atoms.append(int(atom_num))
                    temp_coords.append((float(x), float(y), float(z)))
                self.standard_coords = np.array(temp_coords) * atom_ureg.angstrom
                self.atoms = atoms
            if self.coords.shape == (1, 0):
                self.coords = self.standard_coords

    def _parse_rotation_consts(self):
        if rotational_matches := g16logpatterns["rotation_consts"].search(self.__block):
            rots = []
            for i in range(1, 4):
                try:
                    temp_rot = float(rotational_matches.group(i))
                except:
                    temp_rot = 0
                rots.append(temp_rot)
            self.rotation_constants = np.array(rots) * atom_ureg.Ghertz

    def _parse_functional_basis(self):
        if functional_match := g16logpatterns["functional"].search(self.__block):
            self.functional = functional_match.group(1).lower()
        if basis_match := g16logpatterns["basis"].search(self.__block):
            self.basis = basis_match.group(1)
        if g16logpatterns["Pseudopotential"].search(self.__block):
            self.basis = "pseudopotential"

    def _add_em(self):
        if "em" in self.route_params.keys():
            self.functional = f"{self.functional}-{self.route_params['em']}"
        if "empiricaldispersion" in self.route_params.keys():
            self.functional = (
                f"{self.functional}-{self.route_params['empiricaldispersion']}"
            )

    def _parse_solvent(self):
        solvent_start = g16logpatterns["solvent_start"].search(self.__block)
        solvent_end = g16logpatterns["solvent_end"].search(self.__block)
        if solvent_start and solvent_end:
            solvent_block = self.__block[solvent_start.start() : solvent_end.end()]
            self.__block = self.__block[solvent_end.start() :]
            solvent_model = (
                g16logpatterns["solvent_model"].search(solvent_block).group(1)
            )
            solvent_atom_radii = (
                g16logpatterns["solvent_atomic_radii"].search(solvent_block).group(1)
            )
            if solvent_model == "PCM":
                self.solvent_model = "IEFPCM"
            elif solvent_model == "C-PCM":
                self.solvent_model = "C-PCM"
            if solvent_atom_radii == "SMD-Coulomb":
                self.solvent_model = "SMD-IEFPCM"
            else:
                self.solvent_model = solvent_model
            self.solvent = g16logpatterns["solvent_type"].search(solvent_block).group(1)
            self.solvent_epsilon = float(
                g16logpatterns["solvent_eps"].search(solvent_block).group(1)
            )
            self.solvent_epsilon_infinite = float(
                g16logpatterns["solvent_eps_inf"].search(solvent_block).group(2)
            )
        else:
            self.solvent_model = ""
            self.solvent = ""
            self.solvent_epsilon = None
            self.solvent_epsilon_infinite = None

    def _parse_energy(self) -> Energies:
        energies = {}
        if scf_energy_match := g16logpatterns["scf_energy"].search(self.__block):
            self.__block = self.__block[scf_energy_match.start() :]
            energies["scf_energy"] = (
                float(scf_energy_match.group(1)) * atom_ureg.hartree
            )
            self.status.scf_converged = True
        mp2_4_energy_match = g16logpatterns["mp2-4"].findall(self.__block)
        for mp_energy in mp2_4_energy_match:
            energies[f"{mp_energy[0].lower()}_energy"] = (
                float(mp_energy[1].replace("D", "E")) * atom_ureg.hartree
            )
        if ccsd_energy_match := g16logpatterns["ccsd"].search(self.__block):
            energies["ccsd_energy"] = (
                float(ccsd_energy_match.group(1)) * atom_ureg.hartree
            )
        if ccsd_energy_match := g16logpatterns["ccsd(t)"].search(self.__block):
            energies["ccsd_energy"] = (
                float(ccsd_energy_match.group(1).replace("D", "E")) * atom_ureg.hartree
            )
        return Energies.model_validate(energies)

    def _parse_tail_energies(self) -> Energies:
        energies = {}
        tail = self._tail
        if tail == "":
            return Energies.model_validate(energies)
        groups = tail.split("\\")
        functional = groups[2].lower()
        basis = groups[3].lower()
        if functional[1:] in semi_empirical_methods:
            self.method = "SEMI-EMPIRICAL"
        elif self.functional.endswith("hf"):
            if functional.endswith("hf"):
                self.method = "HF"
            else:
                self.method = "POST-HF"
        else:
            self.method = "DFT"
        self.functional = functional
        if basis != "gen":
            self.basis = basis
        energies_match = g16logpatterns["tail_match"].findall(tail)
        for e, v in energies_match:
            energies_value = float(v) * atom_ureg.hartree
            if "HF" in e:
                energies["scf_energy"] = energies_value
            if "MP2" in e:
                energies["mp2_energy"] = energies_value
            if "MP3" in e:
                energies["mp3_energy"] = energies_value
            if "MP4" in e:
                energies["mp4_energy"] = energies_value
            if "CCSD" in e:
                energies["ccsd_energy"] = energies_value
        return Energies.model_validate(energies)

    def _parse_thermal_energies(self) -> ThermalEnergies:
        thermal_energies = {}
        mappings = {
            "zero-point Energies": "U_0",
            "thermal Energies": "U_T",
            "thermal Enthalpies": "H_T",
            "thermal Free Energies": "G_T",
            ("Zero-point", ""): "ZPVE",
            ("Thermal", " to Energy"): "TCE",
            ("Thermal", " to Enthalpy"): "TCH",
            ("Thermal", " to Gibbs Free Energy"): "TCG",
        }
        matches = g16logpatterns["sum energies"].findall(self.__block)
        if matches:
            for tag, val in matches:
                thermal_energies[mappings[tag]] = (
                    float(val) * atom_ureg.hartree / atom_ureg.particle
                )
        matches = g16logpatterns["corrections"].findall(self.__block)
        if matches:
            for tag1, tag2, val in matches:
                thermal_energies[mappings[(tag1, tag2)]] = (
                    float(val) * atom_ureg.hartree / atom_ureg.particle
                )
        thermal_Cv_S_matches = g16logpatterns["thermal_Cv_S_start"].search(self.__block)
        if thermal_Cv_S_matches:
            block = self.__block[thermal_Cv_S_matches.end() :]
            for line in block.splitlines():
                if "Total" in line:
                    thermal_energies["C_V"] = (
                        float(line.split()[-2])
                        * atom_ureg.calorie
                        / atom_ureg.mol
                        / atom_ureg.kelvin
                    )
                    thermal_energies["S"] = (
                        float(line.split()[-1])
                        * atom_ureg.calorie
                        / atom_ureg.mol
                        / atom_ureg.kelvin
                    )
                    break
        return ThermalEnergies.model_validate(thermal_energies)

    def _parse_tail_thermal_energies(self) -> ThermalEnergies:
        thermal_energies = {}
        thermal_match = g16logpatterns["tail_thermal_match"].findall(self._tail)
        for e, v in thermal_match:
            if "ZeroPoint" in e:
                thermal_energies["ZPVE"] = (
                    float(v) * atom_ureg.hartree / atom_ureg.particle
                )
            if "Thermal" in e:
                thermal_energies["TCE"] = (
                    float(v) * atom_ureg.hartree / atom_ureg.particle
                )
            if "ETot" in e:
                thermal_energies["U_T"] = (
                    float(v) * atom_ureg.hartree / atom_ureg.particle
                )
            if "HTot" in e:
                thermal_energies["H_T"] = (
                    float(v) * atom_ureg.hartree / atom_ureg.particle
                )
            if "GTot" in e:
                thermal_energies["G_T"] = (
                    float(v) * atom_ureg.hartree / atom_ureg.particle
                )
        return ThermalEnergies.model_validate(thermal_energies)

    def _parse_polarizability(self) -> Polarizability:
        polar = {}
        if isotropic_polarizability := g16logpatterns[
            "isotropic_polarizability"
        ].search(self.__block):
            polar["isotropic_polarizability"] = (
                float(isotropic_polarizability.group(1)) * atom_ureg.bohr**3
            )
        if polarizability := g16logpatterns["polarizability"].search(self.__block):
            polar["polarizability"] = (
                np.array(list(map(float, polarizability.groups()))) * atom_ureg.bohr**3
            )
        if electric_dipole_moment := g16logpatterns["electric_dipole_moment"].search(
            self.__block
        ):
            temp_block = self.__block[electric_dipole_moment.start() :].splitlines()
            polar["electric_dipole_moment"] = (
                np.array(
                    [
                        float(temp_block[i].split()[2].replace("D", "E"))
                        for i in range(4, 7)
                    ]
                )
                * atom_ureg.debye
            )
        if polarizability_alter := g16logpatterns["polarizability_alter"].search(
            self.__block
        ):
            temp_block = self.__block[polarizability_alter.start() :].splitlines()
            polar["isotropic_polarizability"] = (
                float(temp_block[4].split()[1].replace("D", "E")) * atom_ureg.bohr**3
            )
            polar["anisotropic_polarizability"] = (
                float(temp_block[5].split()[1].replace("D", "E")) * atom_ureg.bohr**3
            )
            polar["polarizability_tensor"] = (
                np.array(
                    [
                        float(temp_block[i].split()[1].replace("D", "E"))
                        for i in range(6, 12)
                    ]
                )
                * atom_ureg.bohr**3
            )
        polar["dipole"] = self._parse_dipole()
        polar["quadrupole"] = self._parse_quadrupole()
        polar["octapole"] = self._parse_octapole()
        polar["hexadecapole"] = self._parse_hexadecapole()
        polar["electronic_spatial_extent"] = self._parse_electronic_spatial_extent()
        return Polarizability.model_validate(polar)

    def _parse_electronic_spatial_extent(self) -> PlainQuantity:
        ese_match = g16logpatterns["electronic_spatial_extent"].search(self.__block)
        if ese_match:
            return float(ese_match.group(1)) * atom_ureg.bohr**2

    def _parse_dipole(self) -> PlainQuantity:
        start = g16logpatterns["dipole_start"].search(self.__block)
        end = g16logpatterns["quadrupole_start"].search(self.__block)
        if start and end:
            dipole = np.array(
                [
                    float(l)
                    for l in g16logpatterns["dipole"].findall(
                        self.__block[start.start() : end.end()]
                    )
                ]
            )
            return dipole * atom_ureg.debye

    def _parse_quadrupole(self) -> PlainQuantity:
        start = g16logpatterns["quadrupole_start"].search(self.__block)
        end = g16logpatterns["octapole_start"].search(self.__block)
        if start and end:
            quadrupole = np.array(
                [
                    float(l)
                    for l in g16logpatterns["quadrupole"].findall(
                        self.__block[start.start() : end.end()]
                    )
                ]
            )
            return quadrupole * atom_ureg.debye * atom_ureg.angstrom

    def _parse_octapole(self) -> PlainQuantity:
        start = g16logpatterns["octapole_start"].search(self.__block)
        end = g16logpatterns["hexadecapole_start"].search(self.__block)
        if start and end:
            octapole = np.array(
                [
                    float(l)
                    for l in g16logpatterns["octapole"].findall(
                        self.__block[start.start() : end.end()]
                    )
                ]
            )
            return octapole * atom_ureg.debye * atom_ureg.angstrom**2

    def _parse_hexadecapole(self) -> PlainQuantity:
        start = g16logpatterns["hexadecapole_start"].search(self.__block)
        end = g16logpatterns["hexadecapole_end"].search(self.__block)
        if start and end:
            hexadecapole = np.array(
                [
                    float(l)
                    for l in g16logpatterns["hexadecapole"].findall(
                        self.__block[start.start() : end.end()]
                    )
                ]
            )
            return hexadecapole * atom_ureg.debye * atom_ureg.angstrom**3

    def _parse_total_spin(self) -> TotalSpin:
        total_spin = {}
        s = g16logpatterns["spins"].search(self.__block)
        if s:
            total_spin["spin_square"] = float(s.group(1))
            total_spin["spin_quantum_number"] = float(s.group(2))
            if molopconfig.allow_spin_change:
                self.multiplicity = int(
                    round(2 * total_spin["spin_quantum_number"] + 1, 0)
                )
            self.total_spin = TotalSpin.model_validate(total_spin)

    def _parse_orbitals(self):
        orbital_start = g16logpatterns["orbital_start"].search(self.__block)
        orbital_end = g16logpatterns["orbital_end"].search(self.__block)
        if orbital_start and orbital_end:
            block = self.__block[orbital_start.start() : orbital_end.end()]
            self.__block = self.__block[orbital_end.start() :]
            matches = g16logpatterns["orbital"].findall(block)
            temp_alpha_orbitals = []
            temp_alpha_occupancy = []
            temp_beta_orbitals = []
            temp_beta_occupancy = []
            for orbital_type, occ_stat, energies in matches:
                if orbital_type == "Alpha":
                    energies = [
                        e
                        for e in map(
                            float,
                            [energies[j : j + 10] for j in range(0, len(energies), 10)],
                        )
                    ]
                    if occ_stat == "occ.":
                        temp_alpha_occupancy.extend([True for e in energies])
                    else:
                        temp_alpha_occupancy.extend([False for e in energies])
                    temp_alpha_orbitals.extend(energies)
                elif orbital_type == "Beta":
                    energies = [
                        e
                        for e in map(
                            float,
                            [energies[j : j + 10] for j in range(0, len(energies), 10)],
                        )
                    ]
                    if occ_stat == "occ.":
                        temp_beta_occupancy.extend([True for e in energies])
                    else:
                        temp_beta_occupancy.extend([False for e in energies])
                    temp_beta_orbitals.extend(energies)
            if len(temp_beta_orbitals) == 0:
                temp_beta_orbitals = temp_alpha_orbitals
                temp_beta_occupancy = temp_alpha_occupancy
            elif len(temp_beta_orbitals) != len(temp_alpha_orbitals):
                raise ValueError(
                    "Number of beta orbitals does not match number of alpha orbitals"
                )
            self.molecular_orbitals = MolecularOrbitals(
                alpha_energies=np.array(temp_alpha_orbitals) * atom_ureg.hartree,
                beta_energies=np.array(temp_beta_orbitals) * atom_ureg.hartree,
                alpha_occupancies=temp_alpha_occupancy,
                beta_occupancies=temp_beta_occupancy,
            )

    def _parse_mulliken_charges(self):
        mulliken_start = g16logpatterns["mulliken start"].search(self.__block)
        mulliken_end = g16logpatterns["mulliken end"].search(self.__block)
        if mulliken_start and mulliken_end:
            charges = g16logpatterns["mulliken match"].findall(
                self.__block[mulliken_start.start() : mulliken_end.end()]
            )
            return [float(charge) for charge in charges]
        return []

    def _parse_mulliken_spin_density(self):
        spin_densities_start = g16logpatterns["mulliken spin density start"].search(
            self.__block
        )
        spin_densities_end = g16logpatterns["mulliken spin density end"].search(
            self.__block
        )
        if spin_densities_start and spin_densities_end:
            spin_densities = g16logpatterns["mulliken spin density match"].findall(
                self.__block[spin_densities_start.start() : spin_densities_end.end()]
            )
            self.__block = self.__block[spin_densities_end.end() :]
            return [float(spin_density) for spin_density in spin_densities]
        return []

    def _parse_apt_charges(self):
        apt_charges_start = g16logpatterns["apt_start"].search(self.__block)
        apt_charges_end = g16logpatterns["apt_end"].search(self.__block)
        if apt_charges_start and apt_charges_end:
            apt_charges = g16logpatterns["apt_match"].findall(
                self.__block[apt_charges_start.start() : apt_charges_end.end()]
            )
            return [float(charge) for charge in apt_charges]
        return []

    def _parse_lowdin_charges(self):
        lowdin_charges_start = g16logpatterns["lowdin_start"].search(self.__block)
        lowdin_charges_end = g16logpatterns["lowdin_end"].search(self.__block)
        if lowdin_charges_start and lowdin_charges_end:
            lowdin_charges = g16logpatterns["lowdin_match"].findall(
                self.__block[lowdin_charges_start.start() : lowdin_charges_end.end()]
            )
            return [float(charge) for charge in lowdin_charges]
        return []

    def _parse_hirshfeld_spin_charges(self):
        start = g16logpatterns["hirshfeld_start"].search(self.__block)
        end = g16logpatterns["hirshfeld_end"].search(self.__block)
        if start and end:
            hirshfeld = g16logpatterns["hirshfeld"].findall(
                self.__block[start.start() : end.end()]
            )
            return (
                [float(charge) for charge, _, _ in hirshfeld],
                [float(spin) for _, spin, _ in hirshfeld],
                [float(q_cm5) for _, _, q_cm5 in hirshfeld],
            )
        return [], [], []

    def _parse_npa_charges(self):
        start_matches = g16logpatterns[f"npa charge start"].search(self.__block)
        end_matches = g16logpatterns[f"npa charge end"].search(self.__block)
        if start_matches and end_matches:
            temp_block = self.__block[start_matches.start() : end_matches.end()]
            return list(
                map(float, g16logpatterns["npa charge match"].findall(temp_block))
            )
        return []

    def _parse_vibrations(self):
        start = g16logpatterns["freq start"].search(self.__block)
        end = g16logpatterns["freq end"].search(self.__block)
        if start and end:
            block = self.__block[start.end() :]
            self.__block = self.__block[end.start() :]
            freqs = np.array(
                g16logpatterns["freq"].findall(block), dtype=np.float32
            ).flatten()
            red_masses = np.array(
                g16logpatterns["freq red. masses"].findall(block), dtype=np.float32
            ).flatten()
            frc_consts = np.array(
                g16logpatterns["freq frc consts"].findall(block), dtype=np.float32
            ).flatten()
            ir_intens = np.array(
                g16logpatterns["freq IR Inten"].findall(block), dtype=np.float32
            ).flatten()
            freq_modes = np.array(
                g16logpatterns["freq mode"].findall(block), dtype=np.float32
            )
            freq_modes_reindex = []
            for i in range(0, len(freq_modes), 3 * self.n_atom):
                freq_modes_reindex.append(freq_modes[i : i + 3 * self.n_atom : 3])
                freq_modes_reindex.append(
                    freq_modes[i + 1 : i + 1 + 3 * self.n_atom : 3]
                )
                freq_modes_reindex.append(
                    freq_modes[i + 2 : i + 2 + 3 * self.n_atom : 3]
                )
            self.vibrations = Vibrations(
                frequencies=freqs * atom_ureg.cm_1,
                reduced_masses=red_masses * atom_ureg.amu,
                force_constants=frc_consts * atom_ureg.mdyne / atom_ureg.angstrom,
                IR_intensities=ir_intens * atom_ureg.km / atom_ureg.mol,
                vibration_modes=[
                    vib * atom_ureg.angstrom for vib in freq_modes_reindex
                ],
            )

    def _parse_forces(self):
        start = g16logpatterns["forces start"].search(self.__block)
        end = g16logpatterns["forces end"].search(self.__block)
        if start and end:
            forces = g16logpatterns["forces match"].findall(
                self.__block[start.start() : end.end()]
            )
            temp_forces = [(float(fx), float(fy), float(fz)) for fx, fy, fz in forces]
            self.forces = (
                np.array(temp_forces, np.float32) * atom_ureg.hartree / atom_ureg.bohr
            )

    def _parse_temperature(self):
        temperature_match = g16logpatterns["Temperature"].search(self.__block)
        if temperature_match:
            self.temperature = float(temperature_match.group(1)) * atom_ureg.kelvin

    def _parse_status(self):
        feature_mapping = {
            "Maximum Force": "max_force",
            "RMS     Force": "rms_force",
            "Maximum Displacement": "max_displacement",
            "RMS     Displacement": "rms_displacement",
        }
        if "opt" in self.route_params:
            if matches := g16logpatterns["opt stat"].findall(self.__block):
                geometric_metric = {}
                for key, val, threshold, converged in matches:
                    geometric_metric[feature_mapping[key]] = float(val)
                    geometric_metric[f"{feature_mapping[key]}_threshold"] = float(
                        threshold
                    )
                self.geometry_optimization_status = (
                    GeometryOptimizationStatus.model_validate(geometric_metric)
                )
        if matches := g16logpatterns["termination"].findall(self.__block):
            self.status.normal_terminated = matches[0] == "Normal"
            self.status.scf_converged = True
            if "opt" in self.route_params:
                self.geometry_optimization_status.geometry_optimized = True

    def _parse_bond_order(self, _type: Literal["wiberg", "atom_atom_overlap", "mo"]):
        if (
            start_matches := g16logpatterns[f"{_type}_start"].search(self.__block)
        ) and (end_matches := g16logpatterns[f"{_type}_end"].search(self.__block)):
            temp_block = self.__block[start_matches.start() : end_matches.end()]
            self.__block = self.__block[end_matches.start() :]
            digits = list(map(float, g16logpatterns["digit"].findall(temp_block)))
            blks = []
            for idx in range(math.ceil(len(digits) / (self.n_atom * 9))):
                block_length = (
                    min(self.n_atom * 9, len(digits) - idx * self.n_atom * 9)
                    // self.n_atom
                )
                blks.append(
                    np.array(
                        digits[
                            idx * self.n_atom * 9 : idx * self.n_atom * 9
                            + min(self.n_atom * 9, len(digits) - idx * self.n_atom * 9)
                        ]
                    ).reshape(-1, block_length)
                )
            return np.concatenate(blks, axis=1)
        return np.array([[]])

    def _parse_nbo_bond_order(self):
        start_matches = g16logpatterns[f"nbo summary start"].search(self.__block)
        end_matches = g16logpatterns[f"nbo summary end"].search(self.__block)
        if start_matches and end_matches:
            temp_block = self.__block[start_matches.start() : end_matches.end()]
            self.__block = self.__block[end_matches.start() :]
            nbo_bond_order = np.zeros((self.n_atom, self.n_atom))
            for idx1, idx2, occ, energy in g16logpatterns["nbo summary match"].findall(
                temp_block
            ):
                _idx1, _idx2, _occ = (
                    int(idx1),
                    int(idx2),
                    float(occ),
                )
                nbo_bond_order[_idx1 - 1, _idx2 - 1] += _occ
                nbo_bond_order[_idx2 - 1, _idx1 - 1] += _occ
            return nbo_bond_order
        return np.array([[]])

    @property
    def hessian(self) -> np.ndarray:
        """
        Get the hessian matrix (Force constants) from the frame
        """
        hessian_in_body_start = g16logpatterns["hessian_in_body_start"].search(
            self.frame_content
        )
        hessian_in_body_end = g16logpatterns["hessian_in_body_end"].search(
            self.frame_content
        )
        if hessian_in_body_start and hessian_in_body_end:
            hessians = g16logpatterns["hessian_in_body_match"].findall(
                self.frame_content[
                    hessian_in_body_start.start() : hessian_in_body_end.end()
                ]
            )
            return fill_symmetric_matrix(
                np.array(list(map(lambda x: float(x.replace("D", "E")), hessians)))
            )

        tail = self._tail
        if len(tail):
            hessian_in_tail_start = g16logpatterns["hessian_in_tail_start"].search(tail)
            if hessian_in_tail_start:
                hessian_in_tail_end = g16logpatterns["hessian_in_tail_end"].search(
                    tail[hessian_in_tail_start.end() :]
                )
                if hessian_in_tail_end:
                    hessians = g16logpatterns["hessian_in_tail_match"].findall(
                        tail[
                            hessian_in_tail_start.end() : hessian_in_tail_start.end()
                            + hessian_in_tail_end.end()
                        ]
                    )
                    return fill_symmetric_matrix(
                        np.array(list(map(lambda x: float(x), hessians)))
                    )

    @computed_field
    @property
    def task_type(self) -> Literal["sp", "opt", "freq"]:
        if self.route_params:
            if "opt" in self.route_params:
                return "opt"
            if "freq" in self.route_params:
                return "freq"
        return "sp"

    @computed_field
    @property
    def is_error(self) -> bool:
        """
        Check if the current frame is an error frame
        """
        if self.energies.total_energy is None:
            return True
        if not self.status.scf_converged:
            return True
        if self.next is None and not self.status.normal_terminated:
            return True
        return False

    @property
    def is_optimized(self) -> bool:
        return self.geometry_optimization_status.geometry_optimized
