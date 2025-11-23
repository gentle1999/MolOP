from typing import TYPE_CHECKING, Any, Mapping, Optional, Protocol

import numpy as np
from pint.facets.numpy.quantity import NumpyQuantity
from pint.facets.plain import PlainQuantity, PlainUnit
from rdkit import Chem

from molop.config import moloplogger
from molop.io.base_models.DataClasses import (
    ChargeSpinPopulations,
    Energies,
    GeometryOptimizationStatus,
    MolecularOrbitals,
    Polarizability,
    ThermalInformations,
    TotalSpin,
    Vibrations,
)
from molop.io.base_models.FrameParser import BaseFrameParser
from molop.io.base_models.SearchPattern import MolOPPattern
from molop.io.patterns.G16Patterns import g16_log_patterns
from molop.io.QM_models.G16LogFileFrame import (
    G16LogFileFrameDisk,
    G16LogFileFrameMemory,
)
from molop.unit import atom_ureg
from molop.utils.functions import fill_symmetric_matrix, merge_models

pt = Chem.GetPeriodicTable()


def extract_coords(
    coords_match: list[tuple[str | Any, ...]],
) -> tuple[list[int], NumpyQuantity] | tuple[None, None]:
    try:
        atoms: list[int] = []
        coords: list[tuple[float, float, float]] = []
        for row in coords_match:
            atom_num, x, y, z = row
            atoms.append(
                int(atom_num) if atom_num.isdigit() else pt.GetAtomicNumber(atom_num)
            )
            coords.append((float(x), float(y), float(z)))
        coords_ = np.array(coords) * atom_ureg.angstrom
        return atoms, coords_
    except (ValueError, IndexError) as e:
        moloplogger.warning(f"Error extracting coordinates: {e}")
        return None, None


class G16FileFrameParserProtocol(Protocol):
    _block: str
    only_extract_structure: bool
    _file_frame_class_: type[G16LogFileFrameDisk] | type[G16LogFileFrameMemory]


if TYPE_CHECKING:

    class _G16LogFileFrameParserProtocol(G16FileFrameParserProtocol): ...
else:

    class _G16LogFileFrameParserProtocol(object): ...


class G16LogFileFrameParserMixin(_G16LogFileFrameParserProtocol):
    def _parse_frame(self) -> Mapping[str, Any]:
        infos: dict[str, Any] = {"qm_software": "Gaussian"}
        if time := self._parse_time():
            infos["running_time"] = time
        coords_match = self._parse_coords()
        if (atoms := coords_match[0]) and (coords := coords_match[1]) is not None:
            infos["coords"] = coords
            infos["atoms"] = atoms
        coords_match = self._parse_coords_standard_orientation()
        if (atoms := coords_match[0]) and (coords := coords_match[1]) is not None:
            atoms, standard_orientation_coords = coords_match
            infos["standard_coords"] = standard_orientation_coords
            infos["atoms"] = atoms
        if self.only_extract_structure:
            return infos
        if (rotation_consts := self._parse_rotation_consts()) is not None:
            infos["rotation_constants"] = rotation_consts
        energies, total_spin = self._parse_energies_and_total_spin()
        if energies:
            infos["energies"] = energies
        if total_spin:
            infos["total_spin"] = total_spin
        if (polarizability := self._parse_polarizability()) is not None:
            infos["polarizability"] = polarizability
        if pops := self._parse_populations():
            infos.update(pops)
        if (vibrations := self._parse_vibrations()) is not None:
            infos["vibrations"] = vibrations
        if (thermal_info := self._parse_thermal_infos()) is not None:
            infos["thermal_informations"] = thermal_info
        if (forces := self._parse_forces()) is not None:
            infos["forces"] = forces
        if (hessian := self._parse_hessian()) is not None:
            infos["hessian"] = hessian
        if (berny := self._parse_berny()) is not None:
            infos["geometry_optimization_status"] = berny
        if (
            polarizability := self._parse_electric_dipole_moment_and_polarizability()
        ) is not None:
            if "polarizability" in infos:
                infos["polarizability"] = merge_models(
                    infos["polarizability"], polarizability, force_update=True
                )
            else:
                infos["polarizability"] = polarizability
        if (tail := self._parse_tail()) is not None:
            for key, value in tail.items():
                if key not in infos:
                    infos[key] = value
            if (tail_energies := self._parse_tail_energies()) is not None:
                if "energies" in infos:
                    infos["energies"] = merge_models(infos["energies"], tail_energies)
                else:
                    infos["energies"] = tail_energies
            if (tail_thermal_info := self._parse_tail_thermal_infos()) is not None:
                if "thermal_informations" in infos:
                    infos["thermal_informations"] = merge_models(
                        infos["thermal_informations"], tail_thermal_info
                    )
                else:
                    infos["thermal_informations"] = tail_thermal_info
            if (tail_polarizability := self._parse_tail_polarizability()) is not None:
                if "polarizability" in infos:
                    infos["polarizability"] = merge_models(
                        infos["polarizability"], tail_polarizability
                    )
                else:
                    infos["polarizability"] = tail_polarizability
            if (tail_hessian := self._parse_tail_hessian()) is not None:
                if "hessian" not in infos:
                    infos["hessian"] = tail_hessian
        return infos

    def _parse_time(self) -> Optional[PlainQuantity]:
        if time_match := g16_log_patterns.PROCEDURE_TIME.match_content(self._block):
            return (
                sum(
                    float(time_part[1]) + float(time_part[2])
                    for time_part in time_match
                )
                * atom_ureg.second
            )
        return None

    def _parse_coords(self) -> tuple[list[int], NumpyQuantity] | tuple[None, None]:
        focus_content, continued_content = g16_log_patterns.INPUT_COORDS.split_content(
            self._block
        )
        if focus_content == "":
            return None, None
        if coords_match := g16_log_patterns.INPUT_COORDS.get_matches(focus_content):
            self._block = continued_content
            return extract_coords(coords_match)
        return None, None

    def _parse_coords_standard_orientation(
        self,
    ) -> tuple[Optional[list[int]], Optional[NumpyQuantity]]:
        focus_content, self._block = g16_log_patterns.STANDARD_COORDS.split_content(
            self._block
        )
        if focus_content == "":
            return None, None
        if coords_match := g16_log_patterns.STANDARD_COORDS.get_matches(focus_content):
            return extract_coords(coords_match)
        return None, None

    def _parse_rotation_consts(self) -> Optional[NumpyQuantity]:
        if matches := g16_log_patterns.ROTATIONAL_CONST.get_matches(self._block):
            return np.array(list(map(float, matches[0]))) * atom_ureg.gigahertz
        return None

    def _parse_energies_and_total_spin(
        self,
    ) -> tuple[Optional[Energies], Optional[TotalSpin]]:
        scf_energies_dict: dict[str, Optional[PlainQuantity]] = {}
        total_spin_dict: dict[str, Optional[float]] = {}
        focus_content, self._block = g16_log_patterns.SCF_ENERGIES.split_content(
            self._block
        )
        if focus_content == "":
            return (
                (
                    Energies.model_validate(scf_energies_dict)
                    if scf_energies_dict
                    else None
                ),
                TotalSpin.model_validate(total_spin_dict) if total_spin_dict else None,
            )
        if matches := g16_log_patterns.SCF_ENERGY_AND_FUNCTIONAL.get_matches(
            focus_content
        ):
            scf_energies_dict["scf_energy"] = float(matches[0][1]) * atom_ureg.hartree
        if matches := g16_log_patterns.SPIN_SPIN_SQUERE.get_matches(focus_content):
            total_spin_dict["spin_square"] = float(matches[0][0])
            total_spin_dict["spin_quantum_number"] = float(matches[0][1])
        if matches := g16_log_patterns.ENERGY_MP2_4.get_matches(focus_content):
            for match in matches:
                scf_energies_dict[f"{match[0].lower()}_energy"] = (
                    float(match[1].replace("D", "E")) * atom_ureg.hartree
                )
        if matches := g16_log_patterns.ENERGY_MP5.get_matches(focus_content):
            scf_energies_dict["mp5_energy"] = (
                float(matches[0][1].replace("D", "E")) * atom_ureg.hartree
            )
        if matches := g16_log_patterns.ENERGY_CCSD.get_matches(focus_content):
            scf_energies_dict["ccsd_energy"] = (
                float(matches[0][0].replace("D", "E")) * atom_ureg.hartree
            )
        if matches := g16_log_patterns.ENERGY_CCSD_T.get_matches(focus_content):
            scf_energies_dict["ccsd_energy"] = (
                float(matches[0][0].replace("D", "E")) * atom_ureg.hartree
            )
        return (
            Energies.model_validate(scf_energies_dict) if scf_energies_dict else None
        ), (TotalSpin.model_validate(total_spin_dict) if total_spin_dict else None)

    def _parse_polarizability(self) -> Optional[Polarizability]:
        focus_content, self._block = (
            g16_log_patterns.ISOTROPIC_POLARIZABILITY.split_content(self._block)
        )
        if focus_content == "":
            return None
        if matches := g16_log_patterns.ISOTROPIC_POLARIZABILITY.get_matches(
            focus_content
        ):
            return Polarizability(
                isotropic_polarizability=float(matches[0][0]) * atom_ureg.bohr**3
            )
        return None

    def _parse_populations(self) -> dict[str, Any]:
        infos: dict[str, Any] = {}
        mo: dict[str, Any] = {}
        pops: dict[str, Any] = {}
        polars: dict[str, Any] = {}
        try:
            focus_content, self._block = (
                g16_log_patterns.POPULATION_ANALYSIS.split_content(self._block)
            )
            if focus_content == "":
                if mo:
                    infos["molecular_orbitals"] = MolecularOrbitals.model_validate(mo)
                if pops:
                    infos["charge_spin_populations"] = (
                        ChargeSpinPopulations.model_validate(pops)
                    )
                if polars:
                    infos["polarizability"] = Polarizability.model_validate(polars)
                return infos
            patterns_and_keys_1: list[tuple[MolOPPattern, str]] = [
                (
                    g16_log_patterns.MOLECULAR_ORBITALS_SYMMETRY_ALPHA,
                    "alpha_symmetries",
                ),
                (g16_log_patterns.MOLECULAR_ORBITALS_SYMMETRY_BETA, "beta_symmetries"),
                (g16_log_patterns.MOLECULAR_ORBITALS_SYMMETRY, "alpha_symmetries"),
            ]
            for pattern, key in patterns_and_keys_1:
                sub_focus_content, focus_content = pattern.split_content(focus_content)
                if matches := pattern.get_matches(sub_focus_content):
                    mo[key] = [sym for match in matches for sym in match[1].split()]
            if matches := g16_log_patterns.ELECTRONIC_STATE.get_matches(focus_content):
                mo["electronic_state"] = matches[0][0]
            if matches := g16_log_patterns.MOLECULAR_ORBITALS.get_matches(
                focus_content
            ):
                temp_alpha_orbitals = []
                temp_alpha_occupancies = []
                temp_beta_orbitals = []
                temp_beta_occupancies = []
                for orbital_type, occ_stat, energies in matches:
                    energies = [
                        float(energies[j : j + 10]) for j in range(0, len(energies), 10)
                    ]
                    if occ_stat not in ["occ.", "virt."]:
                        raise ValueError(
                            f"Invalid {orbital_type} occupancy status: {occ_stat}"
                        )
                    if orbital_type == "Alpha":
                        temp_alpha_orbitals.extend(energies)
                        temp_alpha_occupancies.extend(
                            [occ_stat == "occ."] * len(energies)
                        )
                    elif orbital_type == "Beta":
                        temp_beta_orbitals.extend(energies)
                        temp_beta_occupancies.extend(
                            [occ_stat == "occ."] * len(energies)
                        )
                    else:
                        raise ValueError(f"Invalid orbital type: {orbital_type}")
                mo["alpha_energies"] = np.array(temp_alpha_orbitals) * atom_ureg.hartree
                mo["alpha_occupancies"] = temp_alpha_occupancies
                mo["beta_energies"] = np.array(temp_beta_orbitals) * atom_ureg.hartree
                mo["beta_occupancies"] = temp_beta_occupancies
            patterns_and_keys_2: list[tuple[MolOPPattern, str]] = [
                (g16_log_patterns.MULLIKEN_SPIN_DENSITY, "mulliken_spins"),
                (g16_log_patterns.MULLIKEN_POPULATION, "mulliken_charges"),
                (g16_log_patterns.APT_POPULATION, "apt_charges"),
                (g16_log_patterns.LOWDIN_POPULATION, "lowdin_charges"),
            ]
            for pattern, key in patterns_and_keys_2:
                sub_focus_content, focus_content = pattern.split_content(focus_content)
                if matches := pattern.get_matches(sub_focus_content):
                    pops[key] = [
                        float(match[0] if key != "mulliken_spins" else match[1])
                        for match in matches
                    ]
            sub_focus_content, focus_content = (
                g16_log_patterns.ELECTRONIC_SPATIAL_EXTENT.split_content(focus_content)
            )
            if matches := g16_log_patterns.ELECTRONIC_SPATIAL_EXTENT.get_matches(
                sub_focus_content
            ):
                polars["electronic_spatial_extent"] = (
                    float(matches[0][0]) * atom_ureg.bohr**2
                )
            patterns_and_keys_3: list[tuple[MolOPPattern, str, PlainUnit]] = [
                (g16_log_patterns.DIPOLE_MOMENT, "dipole", atom_ureg.debye),
                (
                    g16_log_patterns.QUADRUPOLE_MOMENT,
                    "quadrupole",
                    atom_ureg.debye * atom_ureg.angstrom,
                ),
                (
                    g16_log_patterns.TRACELESS_QUADRUPOLE_MOMENT,
                    "traceless_quadrupole",
                    atom_ureg.debye * atom_ureg.angstrom,
                ),
                (
                    g16_log_patterns.OCTAPOLE_MOMENT,
                    "octapole",
                    atom_ureg.debye * atom_ureg.angstrom**2,
                ),
                (
                    g16_log_patterns.HEXADECAPOLE_MOMENT,
                    "hexadecapole",
                    atom_ureg.debye * atom_ureg.angstrom**3,
                ),
            ]
            for pattern, key, unit in patterns_and_keys_3:
                sub_focus_content, focus_content = pattern.split_content(focus_content)
                if matches := pattern.get_matches(sub_focus_content):
                    polars[key] = (
                        np.array([float(match[0]) for match in matches]) * unit
                    )

            if matches := g16_log_patterns.EXACT_POLARIZABILITY.get_matches(
                self._block
            ):
                polars["polarizability_tensor"] = (
                    np.array([float(match) for match in matches[0]]) * atom_ureg.bohr**3
                )

            sub_focus_content, self._block = (
                g16_log_patterns.HIRSHFELD_POPULATION.split_content(self._block)
            )
            if matches := g16_log_patterns.HIRSHFELD_POPULATION.get_matches(
                sub_focus_content
            ):
                pops["hirshfeld_charges"] = [float(match[0]) for match in matches]
                pops["hirshfeld_spins"] = [float(match[1]) for match in matches]
                pops["hirshfeld_q_cm5"] = [float(match[5]) for match in matches]
            if matches := g16_log_patterns.DIPOLE_BEFORE_FORCE.match_content(
                self._block
            ):
                polars["dipole"] = (
                    np.array([float(match[0].replace("D", "E")) for match in matches])
                    * atom_ureg.atomic_unit_of_current
                    * atom_ureg.atomic_unit_of_time
                    * atom_ureg.bohr
                )
            if matches := g16_log_patterns.POLARIZIABILITIES_BEFORE_FORCE.match_content(
                self._block
            ):
                polars["polarizability_tensor"] = (
                    np.array([float(match[0].replace("D", "E")) for match in matches])
                    * atom_ureg.bohr**3
                )
            if mo:
                infos["molecular_orbitals"] = MolecularOrbitals.model_validate(mo)
            if pops:
                infos["charge_spin_populations"] = ChargeSpinPopulations.model_validate(
                    pops
                )
            if polars:
                infos["polarizability"] = Polarizability.model_validate(polars)
        except (ValueError, IndexError) as e:
            moloplogger.warning(
                f"Error parsing populations: {e}, now block: {self._block}"
            )
        except Exception as e:
            moloplogger.error(
                f"Unexpected error occurred while parsing populations: {e}"
            )
        return infos

    def _parse_vibrations(self) -> Optional[Vibrations]:
        vib_dict: dict[str, Any] = {}
        start_index = self._block.find(
            "Harmonic frequencies (cm**-1), IR intensities (KM/Mole), Raman scattering"
        )
        end_index = self._block.find("-------------------", start_index)
        if start_index == -1 or end_index == -1:
            focus_content, self._block = (
                g16_log_patterns.FREQUENCY_ANALYSIS.split_content(self._block)
            )
        else:
            focus_content = self._block[start_index:end_index]
            self._block = self._block[end_index + 1 :]
        if focus_content == "":
            return None
        length = 0
        if matches := g16_log_patterns.FREQUENCIES.get_matches(focus_content):
            vib_dict["frequencies"] = (
                np.array(
                    list(
                        map(
                            float,
                            (freq for match in matches for freq in match[0].split()),
                        )
                    )
                )
                * atom_ureg.cm_1
            )
            length = len(vib_dict["frequencies"])
        if matches := g16_log_patterns.FREQUENCIES_REDUCED_MASS.get_matches(
            focus_content
        ):
            vib_dict["reduced_masses"] = (
                np.array(
                    list(
                        map(
                            float,
                            (freq for match in matches for freq in match[0].split()),
                        )
                    )
                )
                * atom_ureg.amu
            )
        if matches := g16_log_patterns.FREQUENCIES_FORCE_CONSTANTS.get_matches(
            focus_content
        ):
            vib_dict["force_constants"] = (
                np.array(
                    list(
                        map(
                            float,
                            (freq for match in matches for freq in match[0].split()),
                        )
                    )
                )
                * atom_ureg.mdyne
                / atom_ureg.angstrom
            )
        if matches := g16_log_patterns.FREQUENCIES_IR_INTENSITIES.get_matches(
            focus_content
        ):
            vib_dict["IR_intensities"] = (
                np.array(
                    list(
                        map(
                            float,
                            (freq for match in matches for freq in match[0].split()),
                        )
                    )
                )
                * atom_ureg.km
                / atom_ureg.mol
            )
        if matches := g16_log_patterns.FREQUENCIES_MODE.get_matches(focus_content):
            v = np.array(
                [float(freq) for match in matches for freq in match[0].split()]
            ).reshape(-1, 3)
            L = len(v) // length
            v1, v2, v3 = v[0::3], v[1::3], v[2::3]
            vib_dict["vibration_modes"] = [
                vn[i * L : i * L + L] for i in range(length // 3) for vn in [v1, v2, v3]
            ] * atom_ureg.angstrom
        if vib_dict:
            return Vibrations.model_validate(vib_dict)
        return None

    def _parse_thermal_infos(self) -> Optional[ThermalInformations]:
        thermal_dict: dict[str, Any] = {}
        focus_content, self._block = (
            g16_log_patterns.THERMOCHEMISTRY_PART.split_content(self._block)
        )
        if matches := g16_log_patterns.THERMOCHEMISTRY_CORRECTION.get_matches(
            focus_content
        ):
            mapping = {
                ("Zero-point", ""): "ZPVE",
                ("Thermal", " to Energy"): "TCE",
                ("Thermal", " to Enthalpy"): "TCH",
                ("Thermal", " to Gibbs Free Energy"): "TCG",
            }
            for match in matches:
                thermal_dict[mapping[(match[0], match[1])]] = float(
                    match[2]
                ) * atom_ureg.Unit("hartree/particle")
        if matches := g16_log_patterns.THERMOCHEMISTRY_SUM.get_matches(focus_content):
            mapping = {
                "zero-point Energies": "U_0",
                "thermal Energies": "U_T",
                "thermal Enthalpies": "H_T",
                "thermal Free Energies": "G_T",
            }
            for match in matches:
                thermal_dict[mapping[match[0]]] = float(match[1]) * atom_ureg.Unit(
                    "hartree/particle"
                )
        if matches := g16_log_patterns.THERMOCHEMISTRY_CV_S.get_matches(focus_content):
            thermal_dict["S"], thermal_dict["C_V"] = (
                float(matches[0][1]) * atom_ureg.Unit("cal/mol/K"),
                float(matches[0][2]) * atom_ureg.Unit("cal/mol/K"),
            )
        if thermal_dict:
            return ThermalInformations.model_validate(thermal_dict)
        return None

    def _parse_forces(self) -> Optional[NumpyQuantity]:
        focus_content, self._block = g16_log_patterns.FORCES_IN_CARTESIAN.split_content(
            self._block
        )
        if focus_content == "":
            return None
        if matches := g16_log_patterns.FORCES_IN_CARTESIAN.get_matches(focus_content):
            forces = (
                np.array([list(map(float, match)) for match in matches])
                * atom_ureg.hartree
                / atom_ureg.bohr
            )
            return forces
        return None

    def _parse_hessian(self) -> Optional[NumpyQuantity]:
        focus_content, self._block = (
            g16_log_patterns.HESSIAN_IN_CARTESIAN.split_content(self._block)
        )
        if matches := g16_log_patterns.HESSIAN_IN_CARTESIAN.get_matches(focus_content):
            hessian_dict: dict[int, list[float]] = {}
            for match in matches:
                row = int(match[0])
                elements = list(map(float, match[1].replace("D", "E").split()))
                if row not in hessian_dict:
                    hessian_dict[row] = []
                hessian_dict[row].extend(elements)
            hessian = (
                fill_symmetric_matrix(
                    np.array(
                        [element for row in hessian_dict.values() for element in row]
                    )
                )
                * atom_ureg.hartree
                / atom_ureg.bohr**2
            )
            return hessian
        return None

    def _parse_berny(self) -> Optional[GeometryOptimizationStatus]:
        berny_dict: dict[str, float] = {}
        start_index = self._block.find(
            "Item               Value     Threshold  Converged?"
        )
        end_index = self._block.find(
            "GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad",
            start_index,
        )
        if start_index == -1 or end_index == -1:
            focus_content, self._block = (
                g16_log_patterns.BERNY_STATE_MAJOR_PART.split_content(self._block)
            )
            if focus_content == "":
                focus_content, self._block = (
                    g16_log_patterns.BERNY_STATE_BACKUP_PART.split_content(self._block)
                )
        else:
            focus_content = self._block[start_index:end_index]
            self._block = self._block[end_index:]
        if focus_content == "":
            return None
        if matches := g16_log_patterns.BERNY_STATE.get_matches(focus_content):
            mapping = {
                "Maximum Force": "max_force",
                "RMS     Force": "rms_force",
                "Maximum Displacement": "max_displacement",
                "RMS     Displacement": "rms_displacement",
            }
            for match in matches:
                berny_dict[mapping[match[0]]] = float(match[1])
                berny_dict[f"{mapping[match[0]]}_threshold"] = float(match[2])
        if matches := g16_log_patterns.ENERGY_CHANGE.get_matches(focus_content):
            berny_dict["energy_change"] = float(matches[0][0].replace("D", "E"))
        if matches := g16_log_patterns.BERNY_CONCLUSION.get_matches(focus_content):
            berny_dict["geometry_optimized"] = True
        else:
            berny_dict["geometry_optimized"] = False
        if berny_dict:
            return GeometryOptimizationStatus.model_validate(berny_dict)
        return None

    def _parse_electric_dipole_moment_and_polarizability(
        self,
    ) -> Optional[Polarizability]:
        polarizability_dict: dict[str, Any] = {}
        focus_content, self._block = (
            g16_log_patterns.ELECTRIC_DIPOLE_PART.split_content(self._block)
        )
        if focus_content == "":
            return None
        if matches := g16_log_patterns.ELECTRIC_DIPOLE_MOMENT.get_matches(
            focus_content
        ):
            dipole = (
                np.array([float(match[1].replace("D", "E")) for match in matches[1:]])
                * atom_ureg.debye
            )
            polarizability_dict["dipole"] = dipole
        if matches := g16_log_patterns.DIPOLE_POLARIZABILITY.get_matches(focus_content):
            polarizability_dict["isotropic_polarizability"] = (
                float(matches[0][1].replace("D", "E")) * atom_ureg.bohr**3
            )
            polarizability_dict["anisotropic_polarizability"] = (
                float(matches[1][1].replace("D", "E")) * atom_ureg.bohr**3
            )
            polarizability_dict["polarizability_tensor"] = (
                np.array([float(match[1].replace("D", "E")) for match in matches[2:]])
                * atom_ureg.bohr**3
            )
        if polarizability_dict:
            return Polarizability.model_validate(polarizability_dict)
        return None

    def _parse_tail(self) -> dict[str, Any]:
        tail_dict: dict[str, Any] = {}
        focus_content, _ = g16_log_patterns.ARCHIVE_TAIL.split_content(self._block)
        focus_content = focus_content.replace("\n ", "")

        def parse_and_update(pattern: MolOPPattern, key: str):
            nonlocal focus_content
            try:
                sub_focus_content, focus_content = pattern.split_content(focus_content)
                moloplogger.debug(
                    f"{key} focus content: \n{sub_focus_content}\n{key} focus content end"
                )
                if matches := pattern.get_matches(sub_focus_content):
                    tail_dict[key] = matches[0][0]
            except Exception as e:
                moloplogger.error(f"Error in parsing {key}: {e}")

        parse_and_update(g16_log_patterns.JOB_TYPE_IN_ARCHIVE_TAIL, "job_type")
        parse_and_update(g16_log_patterns.FUNCTIONAL_IN_ARCHIVE_TAIL, "functional")
        parse_and_update(g16_log_patterns.BASIS_SET_IN_ARCHIVE_TAIL, "basis")
        parse_and_update(g16_log_patterns.KEYWORDS_IN_ARCHIVE_TAIL, "keywords")
        parse_and_update(g16_log_patterns.TITLE_IN_ARCHIVE_TAIL, "title_card")

        sub_focus_content, focus_content = (
            g16_log_patterns.CHARGE_SPIN_MULTIPLICITY_IN_ARCHIVE_TAIL.split_content(
                focus_content
            )
        )
        if (
            matches
            := g16_log_patterns.CHARGE_SPIN_MULTIPLICITY_IN_ARCHIVE_TAIL.get_matches(
                sub_focus_content
            )
        ):
            charge_multiplicity = matches[0]
            tail_dict["charge"] = int(charge_multiplicity[0])
            tail_dict["multiplicity"] = int(charge_multiplicity[1])
        sub_focus_content, focus_content = (
            g16_log_patterns.COORS_IN_ARCHIVE_TAIL.split_content(focus_content)
        )
        if matches := g16_log_patterns.COORS_IN_ARCHIVE_TAIL.get_matches(
            sub_focus_content
        ):
            tail_dict["atoms"] = [pt.GetAtomicNumber(match[0]) for match in matches]
            tail_dict["coords"] = (
                np.array([list(map(float, match[1:])) for match in matches])
                * atom_ureg.angstrom
            )

        parse_and_update(
            g16_log_patterns.VERSION_IN_ARCHIVE_TAIL, "qm_software_version"
        )
        self._block = focus_content
        return tail_dict

    def _parse_tail_energies(self) -> Optional[Energies]:
        energy_dict: dict[str, Any] = {}
        focus_content, self._block = (
            g16_log_patterns.ENERGIES_IN_ARCHIVE_TAIL.split_content(self._block)
        )
        if focus_content == "":
            return None
        if matches := g16_log_patterns.ENERGIES_IN_ARCHIVE_TAIL.get_matches(
            focus_content
        ):
            for match in matches:
                e = match[0]
                energies_value = float(match[1]) * atom_ureg.hartree
                if "HF" in e:
                    energy_dict["scf_energy"] = energies_value
                if "MP2" in e:
                    energy_dict["mp2_energy"] = energies_value
                if "MP3" in e:
                    energy_dict["mp3_energy"] = energies_value
                if "MP4" in e:
                    energy_dict["mp4_energy"] = energies_value
                if "CCSD" in e:
                    energy_dict["ccsd_energy"] = energies_value
        if energy_dict:
            return Energies.model_validate(energy_dict)
        return None

    def _parse_tail_thermal_infos(self) -> Optional[ThermalInformations]:
        thermal_dict: dict[str, Any] = {}
        thermal_mapping = {
            "ZeroPoint": "ZPVE",
            "Thermal": "TCE",
            "ETot": "U_T",
            "HTot": "H_T",
            "GTot": "G_T",
        }
        focus_content, self._block = (
            g16_log_patterns.THERMOCHEMISTRY_IN_ARCHIVE_TAIL.split_content(self._block)
        )
        if focus_content == "":
            return None
        if matches := g16_log_patterns.THERMOCHEMISTRY_IN_ARCHIVE_TAIL.get_matches(
            focus_content
        ):
            for match in matches:
                if match[0] in thermal_mapping:
                    thermal_dict[thermal_mapping[match[0]]] = float(
                        match[1]
                    ) * atom_ureg.Unit("hartree/particle")
        if thermal_dict:
            return ThermalInformations.model_validate(thermal_dict)
        return None

    def _parse_tail_polarizability(self) -> Optional[Polarizability]:
        polarizability_dict: dict[str, Any] = {}
        if matches := g16_log_patterns.DIPOLE_IN_ARCHIVE_TAIL.get_matches(self._block):
            polarizability_dict["dipole"] = (
                np.array([float(match[1]) for match in matches[1:]]) * atom_ureg.debye
            )
        if matches := g16_log_patterns.POLAR_IN_ARCHIVE_TAIL.get_matches(self._block):
            polarizability_dict["polarizability_tensor"] = (
                np.array([float(match[1]) for match in matches[1:]]) * atom_ureg.bohr**3
            )
        if matches := g16_log_patterns.QUADRUPOLE_IN_ARCHIVE_TAIL.get_matches(
            self._block
        ):
            polarizability_dict["quadrupole"] = (
                np.array([float(match[1]) for match in matches[1:]]) * atom_ureg.bohr**3
            )
        if polarizability_dict:
            return Polarizability.model_validate(polarizability_dict)
        return None

    def _parse_tail_hessian(self) -> Optional[NumpyQuantity]:
        focus_content, self._block = (
            g16_log_patterns.HESSIAN_IN_ARCHIVE_TAIL.split_content(self._block)
        )
        if matches := g16_log_patterns.HESSIAN_IN_ARCHIVE_TAIL.get_matches(
            focus_content
        ):
            return (
                fill_symmetric_matrix(np.array([float(match[0]) for match in matches]))
                * atom_ureg.hartree
                / atom_ureg.bohr**2
            )
        return None


class G16LogFileFrameParserMemory(
    G16LogFileFrameParserMixin, BaseFrameParser[G16LogFileFrameMemory]
):
    _file_frame_class_ = G16LogFileFrameMemory


class G16LogFileFrameParserDisk(
    G16LogFileFrameParserMixin, BaseFrameParser[G16LogFileFrameDisk]
):
    _file_frame_class_ = G16LogFileFrameDisk
