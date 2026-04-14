from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from pint.facets.numpy.quantity import NumpyQuantity
from pint.facets.plain import PlainQuantity, PlainUnit

from molop.config import moloplogger
from molop.io.base_models.DataClasses import GeometryOptimizationStatus
from molop.io.base_models.SearchPattern import MolOPPattern, MolOPPatternV2
from molop.io.logic.QM_frame_models.G16V3Components import (
    G16V3L202RotConstComponent,
    G16V3L9999ArchiveComponent,
)
from molop.io.logic.QM_frame_parsers._g16_v2_shared import (
    ARCHIVE_TAIL_V2,
    BERNY_STATE_BACKUP_PART_V2,
    BERNY_STATE_MAJOR_PART_V2,
    ELECTRIC_DIPOLE_PART_V2,
    ENERGIES_IN_ARCHIVE_TAIL_V2,
    FORCES_IN_CARTESIAN_V2,
    FREQUENCY_ANALYSIS_V2,
    HESSIAN_IN_ARCHIVE_TAIL_V2,
    HESSIAN_IN_CARTESIAN_V2,
    INPUT_COORDS_V2,
    ISOTROPIC_POLARIZABILITY_V2,
    POPULATION_ANALYSIS_V2,
    SCF_ENERGIES_V2,
    STANDARD_COORDS_V2,
    THERMOCHEMISTRY_IN_ARCHIVE_TAIL_V2,
    THERMOCHEMISTRY_PART_V2,
    _extract_molecular_orbital_payload_from_text,
    _summarize_parse_context,
    _trim_molecular_orbital_symmetries,
    extract_coords,
)
from molop.io.logic.QM_parsers._g16log_archive_tail import parse_archive_tail_metadata
from molop.io.patterns.G16Patterns import g16_log_patterns
from molop.unit import atom_ureg
from molop.utils.functions import fill_symmetric_matrix


@dataclass(slots=True)
class ParseState:
    content: str
    cursor: int = 0
    only_extract_structure: bool = False

    @property
    def remaining_content(self) -> str:
        return self.content[self.cursor :]

    def advance_to(self, next_cursor: int) -> None:
        self.cursor = max(self.cursor, next_cursor)


def _focus_from_state(pattern: MolOPPatternV2, state: ParseState) -> tuple[str, int]:
    return pattern.split_content_from(state.content, state.cursor)


def extract_input_coords_from_state(
    state: ParseState,
) -> tuple[list[int] | None, NumpyQuantity | None]:
    focus_content, next_cursor = _focus_from_state(INPUT_COORDS_V2, state)
    if focus_content == "":
        return None, None
    if coords_match := g16_log_patterns.INPUT_COORDS.get_matches(focus_content):
        state.advance_to(next_cursor)
        return extract_coords(coords_match)
    return None, None


def extract_standard_coords_from_state(
    state: ParseState,
) -> tuple[list[int] | None, NumpyQuantity | None]:
    focus_content, next_cursor = _focus_from_state(STANDARD_COORDS_V2, state)
    if focus_content == "":
        return None, None
    if coords_match := g16_log_patterns.STANDARD_COORDS.get_matches(focus_content):
        state.advance_to(next_cursor)
        return extract_coords(coords_match)
    state.advance_to(next_cursor)
    return None, None


def extract_energies_and_total_spin_from_state(
    state: ParseState,
) -> tuple[dict[str, Any] | None, dict[str, Any] | None]:
    scf_energies_dict: dict[str, PlainQuantity | None] = {}
    total_spin_dict: dict[str, float | None] = {}
    focus_content, next_cursor = _focus_from_state(SCF_ENERGIES_V2, state)
    if focus_content == "":
        return None, None
    state.advance_to(next_cursor)
    if matches := g16_log_patterns.SCF_ENERGY_AND_FUNCTIONAL.get_matches(focus_content):
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
        scf_energies_dict["mp5_energy"] = float(matches[0][1].replace("D", "E")) * atom_ureg.hartree
    if matches := g16_log_patterns.ENERGY_CCSD.get_matches(focus_content):
        scf_energies_dict["ccsd_energy"] = (
            float(matches[0][0].replace("D", "E")) * atom_ureg.hartree
        )
    if matches := g16_log_patterns.ENERGY_CCSD_T.get_matches(focus_content):
        scf_energies_dict["ccsd_energy"] = (
            float(matches[0][0].replace("D", "E")) * atom_ureg.hartree
        )
    return scf_energies_dict or None, total_spin_dict or None


def extract_polarizability_from_state(state: ParseState) -> dict[str, Any] | None:
    focus_content, next_cursor = _focus_from_state(ISOTROPIC_POLARIZABILITY_V2, state)
    if focus_content == "":
        return None
    state.advance_to(next_cursor)
    if matches := g16_log_patterns.ISOTROPIC_POLARIZABILITY.get_matches(focus_content):
        return {"isotropic_polarizability": float(matches[0][0]) * atom_ureg.bohr**3}
    return None


def extract_populations_from_state(state: ParseState) -> dict[str, Any]:
    infos: dict[str, Any] = {}
    mo: dict[str, Any] = {}
    pops: dict[str, Any] = {}
    polars: dict[str, Any] = {}
    focus_content, next_cursor = _focus_from_state(POPULATION_ANALYSIS_V2, state)
    if focus_content == "":
        return infos
    state.advance_to(next_cursor)
    remainder_content = state.remaining_content
    try:
        patterns_and_keys_1: list[tuple[MolOPPattern, str]] = [
            (g16_log_patterns.MOLECULAR_ORBITALS_SYMMETRY_ALPHA, "alpha_symmetries"),
            (g16_log_patterns.MOLECULAR_ORBITALS_SYMMETRY_BETA, "beta_symmetries"),
            (g16_log_patterns.MOLECULAR_ORBITALS_SYMMETRY, "alpha_symmetries"),
        ]
        for pattern, key in patterns_and_keys_1:
            sub_focus_content, focus_content = pattern.split_content(focus_content)
            if matches := pattern.get_matches(sub_focus_content):
                mo[key] = [sym for match in matches for sym in match[1].split()]
        mo.update(_extract_molecular_orbital_payload_from_text(focus_content))
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
                    float(match[0] if key != "mulliken_spins" else match[1]) for match in matches
                ]
        sub_focus_content, focus_content = g16_log_patterns.ELECTRONIC_SPATIAL_EXTENT.split_content(
            focus_content
        )
        if matches := g16_log_patterns.ELECTRONIC_SPATIAL_EXTENT.get_matches(sub_focus_content):
            polars["electronic_spatial_extent"] = float(matches[0][0]) * atom_ureg.bohr**2
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
            (g16_log_patterns.OCTAPOLE_MOMENT, "octapole", atom_ureg.debye * atom_ureg.angstrom**2),
            (
                g16_log_patterns.HEXADECAPOLE_MOMENT,
                "hexadecapole",
                atom_ureg.debye * atom_ureg.angstrom**3,
            ),
        ]
        for pattern, key, unit in patterns_and_keys_3:
            sub_focus_content, focus_content = pattern.split_content(focus_content)
            if matches := pattern.get_matches(sub_focus_content):
                polars[key] = np.array([float(match[0]) for match in matches]) * unit
        if matches := g16_log_patterns.EXACT_POLARIZABILITY.get_matches(remainder_content):
            polars["polarizability_tensor"] = (
                np.array([float(match) for match in matches[0]]) * atom_ureg.bohr**3
            )
        sub_focus_content, _ignored = g16_log_patterns.HIRSHFELD_POPULATION.split_content(
            remainder_content
        )
        if matches := g16_log_patterns.HIRSHFELD_POPULATION.get_matches(sub_focus_content):
            pops["hirshfeld_charges"] = [float(match[0]) for match in matches]
            pops["hirshfeld_spins"] = [float(match[1]) for match in matches]
            pops["hirshfeld_q_cm5"] = [float(match[5]) for match in matches]
        if matches := g16_log_patterns.DIPOLE_BEFORE_FORCE.match_content(remainder_content):
            polars["dipole"] = (
                np.array([float(match[0].replace("D", "E")) for match in matches])
                * atom_ureg.atomic_unit_of_current
                * atom_ureg.atomic_unit_of_time
                * atom_ureg.bohr
            )
        if matches := g16_log_patterns.POLARIZIABILITIES_BEFORE_FORCE.match_content(
            remainder_content
        ):
            polars["polarizability_tensor"] = (
                np.array([float(match[0].replace("D", "E")) for match in matches])
                * atom_ureg.bohr**3
            )
        if mo:
            mo = _trim_molecular_orbital_symmetries(mo)
            infos["molecular_orbitals"] = mo
        if pops:
            infos["charge_spin_populations"] = pops
        if polars:
            infos["polarizability"] = polars
    except (ValueError, IndexError) as exc:
        moloplogger.warning(
            "Error parsing populations: %s | context=%s",
            exc,
            _summarize_parse_context(remainder_content),
        )
    except Exception as exc:
        moloplogger.error(f"Unexpected error occurred while parsing populations: {exc}")
    return infos


def extract_vibrations_from_state(state: ParseState) -> dict[str, Any] | None:
    block = state.content
    start_index = block.find(
        "Harmonic frequencies (cm**-1), IR intensities (KM/Mole), Raman scattering",
        state.cursor,
    )
    end_index = block.find("-------------------", start_index)
    if start_index == -1 or end_index == -1:
        focus_content, next_cursor = _focus_from_state(FREQUENCY_ANALYSIS_V2, state)
    else:
        focus_content = block[start_index:end_index]
        next_cursor = end_index + 1
    if focus_content == "":
        return None
    state.advance_to(next_cursor)
    vib_dict: dict[str, Any] = {}
    length = 0
    if matches := g16_log_patterns.FREQUENCIES.get_matches(focus_content):
        vib_dict["frequencies"] = (
            np.array(list(map(float, (freq for match in matches for freq in match[0].split()))))
            * atom_ureg.cm_1
        )
        length = len(vib_dict["frequencies"])
    if matches := g16_log_patterns.FREQUENCIES_REDUCED_MASS.get_matches(focus_content):
        vib_dict["reduced_masses"] = (
            np.array(list(map(float, (freq for match in matches for freq in match[0].split()))))
            * atom_ureg.amu
        )
    if matches := g16_log_patterns.FREQUENCIES_FORCE_CONSTANTS.get_matches(focus_content):
        vib_dict["force_constants"] = (
            np.array(list(map(float, (freq for match in matches for freq in match[0].split()))))
            * atom_ureg.mdyne
            / atom_ureg.angstrom
        )
    if matches := g16_log_patterns.FREQUENCIES_IR_INTENSITIES.get_matches(focus_content):
        vib_dict["IR_intensities"] = (
            np.array(list(map(float, (freq for match in matches for freq in match[0].split()))))
            * atom_ureg.km
            / atom_ureg.mol
        )
    if matches := g16_log_patterns.FREQUENCIES_MODE.get_matches(focus_content):
        v = np.array([float(freq) for match in matches for freq in match[0].split()]).reshape(-1, 3)
        L = len(v) // length
        v1, v2, v3 = v[0::3], v[1::3], v[2::3]
        vib_dict["vibration_modes"] = [
            vn[i * L : i * L + L] for i in range(length // 3) for vn in [v1, v2, v3]
        ] * atom_ureg.angstrom
    return vib_dict or None


def extract_thermal_infos_from_state(state: ParseState) -> dict[str, Any] | None:
    focus_content, next_cursor = _focus_from_state(THERMOCHEMISTRY_PART_V2, state)
    if focus_content == "":
        return None
    state.advance_to(next_cursor)
    thermal_dict: dict[str, Any] = {}
    if matches := g16_log_patterns.THERMOCHEMISTRY_CORRECTION.get_matches(focus_content):
        mapping = {
            ("Zero-point", ""): "ZPVE",
            ("Thermal", " to Energy"): "TCE",
            ("Thermal", " to Enthalpy"): "TCH",
            ("Thermal", " to Gibbs Free Energy"): "TCG",
        }
        for match in matches:
            thermal_dict[mapping[(match[0], match[1])]] = float(match[2]) * atom_ureg.Unit(
                "hartree/particle"
            )
    if matches := g16_log_patterns.THERMOCHEMISTRY_SUM.get_matches(focus_content):
        mapping = {
            "zero-point Energies": "U_0",
            "thermal Energies": "U_T",
            "thermal Enthalpies": "H_T",
            "thermal Free Energies": "G_T",
        }
        for match in matches:
            thermal_dict[mapping[match[0]]] = float(match[1]) * atom_ureg.Unit("hartree/particle")
    if matches := g16_log_patterns.THERMOCHEMISTRY_CV_S.get_matches(focus_content):
        thermal_dict["S"], thermal_dict["C_V"] = (
            float(matches[0][1]) * atom_ureg.Unit("cal/mol/K"),
            float(matches[0][2]) * atom_ureg.Unit("cal/mol/K"),
        )
    return thermal_dict or None


def extract_forces_from_state(state: ParseState) -> NumpyQuantity | None:
    focus_content, next_cursor = _focus_from_state(FORCES_IN_CARTESIAN_V2, state)
    if focus_content == "":
        return None
    state.advance_to(next_cursor)
    if matches := g16_log_patterns.FORCES_IN_CARTESIAN.get_matches(focus_content):
        return (
            np.array([list(map(float, match)) for match in matches])
            * atom_ureg.hartree
            / atom_ureg.bohr
        )
    return None


def extract_hessian_from_state(state: ParseState) -> NumpyQuantity | None:
    focus_content, next_cursor = _focus_from_state(HESSIAN_IN_CARTESIAN_V2, state)
    if focus_content == "":
        return None
    state.advance_to(next_cursor)
    if matches := g16_log_patterns.HESSIAN_IN_CARTESIAN.get_matches(focus_content):
        hessian_dict: dict[int, list[float]] = {}
        for match in matches:
            row = int(match[0])
            elements = list(map(float, match[1].replace("D", "E").split()))
            if row not in hessian_dict:
                hessian_dict[row] = []
            hessian_dict[row].extend(elements)
        return (
            fill_symmetric_matrix(
                np.array([element for row in hessian_dict.values() for element in row])
            )
            * atom_ureg.hartree
            / atom_ureg.bohr**2
        )
    return None


def extract_berny_from_state(state: ParseState) -> GeometryOptimizationStatus | None:
    berny_dict: dict[str, float | bool] = {}
    start_index = state.content.find(
        "Item               Value     Threshold  Converged?", state.cursor
    )
    end_index = state.content.find(
        "GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad",
        start_index,
    )
    if start_index == -1 or end_index == -1:
        focus_content, next_cursor = _focus_from_state(BERNY_STATE_MAJOR_PART_V2, state)
        if focus_content == "":
            focus_content, next_cursor = _focus_from_state(BERNY_STATE_BACKUP_PART_V2, state)
    else:
        focus_content = state.content[start_index:end_index]
        next_cursor = end_index
    if focus_content == "":
        return None
    state.advance_to(next_cursor)
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
    berny_dict["geometry_optimized"] = bool(
        g16_log_patterns.BERNY_CONCLUSION.get_matches(focus_content)
    )
    return GeometryOptimizationStatus.model_validate(berny_dict) if berny_dict else None


def extract_electric_dipole_and_polarizability_from_state(
    state: ParseState,
) -> dict[str, Any] | None:
    focus_content, next_cursor = _focus_from_state(ELECTRIC_DIPOLE_PART_V2, state)
    if focus_content == "":
        return None
    state.advance_to(next_cursor)
    polarizability_dict: dict[str, Any] = {}
    if matches := g16_log_patterns.ELECTRIC_DIPOLE_MOMENT.get_matches(focus_content):
        polarizability_dict["dipole"] = (
            np.array([float(match[1].replace("D", "E")) for match in matches[1:]]) * atom_ureg.debye
        )
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
    return polarizability_dict or None


def extract_tail_metadata_from_state(state: ParseState) -> dict[str, Any]:
    focus_content, next_cursor = _focus_from_state(ARCHIVE_TAIL_V2, state)
    if focus_content == "":
        return {}
    state.advance_to(next_cursor)
    tail, _tail_remaining = parse_archive_tail_metadata(focus_content, include_coords=True)
    return tail


def extract_tail_energies_from_state(state: ParseState) -> dict[str, Any] | None:
    focus_content, next_cursor = _focus_from_state(ENERGIES_IN_ARCHIVE_TAIL_V2, state)
    if focus_content == "":
        return None
    state.advance_to(next_cursor)
    energy_dict: dict[str, Any] = {}
    if matches := g16_log_patterns.ENERGIES_IN_ARCHIVE_TAIL.get_matches(focus_content):
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
    return energy_dict or None


def extract_tail_thermal_infos_from_state(state: ParseState) -> dict[str, Any] | None:
    focus_content, next_cursor = _focus_from_state(THERMOCHEMISTRY_IN_ARCHIVE_TAIL_V2, state)
    if focus_content == "":
        return None
    state.advance_to(next_cursor)
    thermal_dict: dict[str, Any] = {}
    thermal_mapping = {
        "ZeroPoint": "ZPVE",
        "Thermal": "TCE",
        "ETot": "U_T",
        "HTot": "H_T",
        "GTot": "G_T",
    }
    if matches := g16_log_patterns.THERMOCHEMISTRY_IN_ARCHIVE_TAIL.get_matches(focus_content):
        for match in matches:
            if match[0] in thermal_mapping:
                thermal_dict[thermal_mapping[match[0]]] = float(match[1]) * atom_ureg.Unit(
                    "hartree/particle"
                )
    return thermal_dict or None


def extract_tail_polarizability_from_state(state: ParseState) -> dict[str, Any] | None:
    return G16V3L9999ArchiveComponent._extract_archive_polarizability(state.remaining_content)


def extract_tail_hessian_from_state(state: ParseState) -> NumpyQuantity | None:
    focus_content, next_cursor = _focus_from_state(HESSIAN_IN_ARCHIVE_TAIL_V2, state)
    if focus_content == "":
        return None
    state.advance_to(next_cursor)
    if matches := g16_log_patterns.HESSIAN_IN_ARCHIVE_TAIL.get_matches(focus_content):
        return (
            fill_symmetric_matrix(np.array([float(match[0]) for match in matches]))
            * atom_ureg.hartree
            / atom_ureg.bohr**2
        )
    return None


def extract_rotation_consts_from_state(state: ParseState) -> NumpyQuantity | None:
    return G16V3L202RotConstComponent._extract_rotation_constants(state.remaining_content)
