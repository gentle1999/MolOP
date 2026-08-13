from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np
from pint.facets.numpy.quantity import NumpyQuantity
from pint.facets.plain import PlainQuantity, PlainUnit

from molop.config import moloplogger
from molop.io.base_models.DataClasses import EnergyObservation
from molop.io.base_models.SearchPattern import MolOPPattern
from molop.io.logic.gaussian.log.frame_parsers._g16_shared import (
    ARCHIVE_TAIL,
    BERNY_STATE_BACKUP_PART,
    BERNY_STATE_MAJOR_PART,
    ELECTRIC_DIPOLE_PART,
    ENERGIES_IN_ARCHIVE_TAIL,
    FORCES_IN_CARTESIAN,
    FREQUENCY_ANALYSIS,
    HESSIAN_IN_ARCHIVE_TAIL,
    HESSIAN_IN_CARTESIAN,
    INPUT_COORDS,
    ISOTROPIC_POLARIZABILITY,
    POPULATION_ANALYSIS,
    SCF_ENERGIES,
    STANDARD_COORDS,
    THERMOCHEMISTRY_IN_ARCHIVE_TAIL,
    THERMOCHEMISTRY_PART,
    _extract_float_tokens,
    _extract_labeled_float_tokens,
    _extract_molecular_orbital_payload_from_text,
    _summarize_parse_context,
    _trim_molecular_orbital_symmetries,
    extract_coords,
    extract_rotation_constants,
    parse_rotational_constant,
)
from molop.io.logic.gaussian.log.parsers._g16_log_patterns import g16_log_patterns
from molop.io.logic.gaussian.log.parsers._g16log_archive_tail import (
    parse_archive_tail_energies,
    parse_archive_tail_hessian,
    parse_archive_tail_payload,
    parse_archive_tail_polarizability,
    parse_archive_tail_thermal_infos,
)
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


def _focus_from_state(pattern: MolOPPattern, state: ParseState) -> tuple[str, int]:
    return pattern.split_content_from(state.content, state.cursor)


def _last_pattern_matches(pattern: MolOPPattern, content: str) -> list[Any]:
    """Return matches from the last complete occurrence of a delimited table."""

    cursor = 0
    latest: list[Any] = []
    while located := pattern.locate_content_from(content, cursor):
        start_start, _start_end, _end_start, end_end = located
        matches = pattern.find_matches(content[start_start:end_end])
        if matches:
            latest = matches
        if end_end <= cursor:
            break
        cursor = end_end
    return latest


def extract_input_coords_from_state(
    state: ParseState,
) -> tuple[list[int] | None, NumpyQuantity | None]:
    focus_content, next_cursor = _focus_from_state(INPUT_COORDS, state)
    if focus_content == "":
        return None, None
    if coords_match := g16_log_patterns.INPUT_COORDS.find_matches(focus_content):
        state.advance_to(next_cursor)
        return extract_coords(coords_match)
    return None, None


def extract_standard_coords_from_state(
    state: ParseState,
) -> tuple[list[int] | None, NumpyQuantity | None]:
    focus_content, next_cursor = _focus_from_state(STANDARD_COORDS, state)
    if focus_content == "":
        return None, None
    if coords_match := g16_log_patterns.STANDARD_COORDS.find_matches(focus_content):
        state.advance_to(next_cursor)
        return extract_coords(coords_match)
    state.advance_to(next_cursor)
    return None, None


def extract_energies_and_total_spin_from_state(
    state: ParseState,
    *,
    capture_source_evidence: bool = False,
) -> tuple[dict[str, Any] | None, dict[str, Any] | None]:
    scf_energies_dict: dict[str, Any] = {}
    total_spin_dict: dict[str, float | None] = {}
    observations: dict[str, EnergyObservation] = {}

    def observe(method: str, value: PlainQuantity, source_label: str) -> None:
        if capture_source_evidence:
            observations[method] = EnergyObservation(
                method=method,
                quantity_semantics="total_energy",
                value=value,
                source_label=source_label,
            )

    focus_content, next_cursor = _focus_from_state(SCF_ENERGIES, state)
    external_matches = g16_log_patterns.EXTERNAL_ENERGY_RESULT.find_matches(state.remaining_content)
    if focus_content == "" and not external_matches:
        return None, None
    if focus_content:
        state.advance_to(next_cursor)
    else:
        focus_content = state.remaining_content
    if matches := g16_log_patterns.SCF_ENERGY_AND_FUNCTIONAL.find_matches(focus_content):
        reference_energy = float(matches[0].group("energy")) * atom_ureg.hartree
        scf_energies_dict["reference_energy"] = reference_energy
        observe("reference", reference_energy, "SCF Done")
    elif external_matches:
        reference_energy = (
            float(external_matches[0].group("energy").replace("D", "E").replace("d", "e"))
            * atom_ureg.hartree
        )
        scf_energies_dict["reference_energy"] = reference_energy
        observe("reference", reference_energy, "Gaussian External Energy")
    if matches := g16_log_patterns.SPIN_SPIN_SQUERE.find_matches(focus_content):
        total_spin_dict["spin_square"] = float(matches[0].group("spin_square"))
        total_spin_dict["spin_quantum_number"] = float(matches[0].group("spin_quantum_number"))
    if matches := g16_log_patterns.ENERGY_MP2_4.find_matches(focus_content):
        for matched in matches:
            method = matched.group("method").upper()
            energy = float(matched.group("energy").replace("D", "E")) * atom_ureg.hartree
            scf_energies_dict[f"{method.lower()}_energy"] = energy
            observe(method, energy, f"EU{method}")
    if matches := g16_log_patterns.ENERGY_MP5.find_matches(focus_content):
        mp5_energy = float(matches[0].group("energy").replace("D", "E")) * atom_ureg.hartree
        scf_energies_dict["mp5_energy"] = mp5_energy
        observe("MP5", mp5_energy, "MP5")
    if matches := g16_log_patterns.ENERGY_CCSD.find_matches(focus_content):
        ccsd_energy = float(matches[0].group("energy").replace("D", "E")) * atom_ureg.hartree
        scf_energies_dict["ccsd_energy"] = ccsd_energy
        observe("CCSD", ccsd_energy, "Wavefunction amplitudes converged. E(Corr)")
    if matches := g16_log_patterns.ENERGY_CCSD_T.find_matches(focus_content):
        ccsd_t_energy = float(matches[0].group("energy").replace("D", "E")) * atom_ureg.hartree
        scf_energies_dict["ccsd_t_energy"] = ccsd_t_energy
        observe("CCSD(T)", ccsd_t_energy, "CCSD(T)")
    if observations:
        scf_energies_dict["observations"] = list(observations.values())
    return scf_energies_dict or None, total_spin_dict or None


def extract_polarizability_from_state(state: ParseState) -> dict[str, Any] | None:
    focus_content, next_cursor = _focus_from_state(ISOTROPIC_POLARIZABILITY, state)
    if focus_content == "":
        return None
    state.advance_to(next_cursor)
    if matches := g16_log_patterns.ISOTROPIC_POLARIZABILITY.find_matches(focus_content):
        return {"isotropic_polarizability": float(matches[0].group("value")) * atom_ureg.bohr**3}
    return None


def extract_nmr_from_state(state: ParseState) -> dict[str, Any] | None:
    focus_content, next_cursor = _focus_from_state(g16_log_patterns.NMR_SHIELDING, state)
    if focus_content == "":
        return None
    state.advance_to(next_cursor)

    header_matches = g16_log_patterns.NMR_SHIELDING_HEADER.find_matches(focus_content)
    gauge = header_matches[0].group("gauge") if header_matches else None
    shielding_tensors: list[dict[str, Any]] = []
    component_names = (("xx", "xy", "xz"), ("yx", "yy", "yz"), ("zx", "zy", "zz"))
    for matched in g16_log_patterns.NMR_SHIELDING.find_matches(focus_content):
        tensor = np.asarray(
            [[float(matched.group(component)) for component in row] for row in component_names]
        )
        shielding_tensors.append(
            {
                "atom_index": int(matched.group("atom_index")) - 1,
                "atom_symbol": matched.group("atom_symbol"),
                "shielding_tensor": tensor * atom_ureg.ppm,
                "isotropic": float(matched.group("isotropic")) * atom_ureg.ppm,
                "anisotropy": float(matched.group("anisotropy")) * atom_ureg.ppm,
                "principal_values": np.asarray(
                    [
                        float(matched.group("eigenvalue_1")),
                        float(matched.group("eigenvalue_2")),
                        float(matched.group("eigenvalue_3")),
                    ]
                )
                * atom_ureg.ppm,
                "anisotropy_convention": "Gaussian",
                "orientation": "unknown",
            }
        )
    if not shielding_tensors:
        return None

    payload: dict[str, Any] = {"gauge": gauge, "shielding_tensors": shielding_tensors}
    coupling_size = max(item["atom_index"] for item in shielding_tensors) + 1
    component_patterns = {
        "spin_spin_coupling_k_components": (
            ("FC", g16_log_patterns.NMR_COUPLING_FC_K),
            ("SD", g16_log_patterns.NMR_COUPLING_SD_K),
            ("PSO", g16_log_patterns.NMR_COUPLING_PSO_K),
            ("DSO", g16_log_patterns.NMR_COUPLING_DSO_K),
        ),
        "spin_spin_coupling_j_components": (
            ("FC", g16_log_patterns.NMR_COUPLING_FC_J),
            ("SD", g16_log_patterns.NMR_COUPLING_SD_J),
            ("PSO", g16_log_patterns.NMR_COUPLING_PSO_J),
            ("DSO", g16_log_patterns.NMR_COUPLING_DSO_J),
        ),
    }
    for field_name, patterns in component_patterns.items():
        components: dict[str, Any] = {}
        for component_name, pattern in patterns:
            coupling_content, _ = pattern.split_content(focus_content)
            matrix = _extract_nmr_coupling_matrix(coupling_content, coupling_size)
            if matrix is not None:
                components[component_name] = matrix * atom_ureg.Hz
        if components:
            payload[field_name] = components
    for pattern, field_name in (
        (g16_log_patterns.NMR_TOTAL_COUPLING_K, "spin_spin_coupling_k"),
        (g16_log_patterns.NMR_TOTAL_COUPLING_J, "spin_spin_coupling_j"),
    ):
        coupling_content, _ = pattern.split_content(focus_content)
        matrix = _extract_nmr_coupling_matrix(coupling_content, coupling_size)
        if matrix is not None:
            payload[field_name] = matrix * atom_ureg.Hz
    if any(key.startswith("spin_spin_coupling_") for key in payload):
        payload["coupling_atom_indices"] = list(range(coupling_size))
    return payload


def _extract_nmr_coupling_matrix(content: str, size: int) -> np.ndarray | None:
    if not content or size <= 0:
        return None
    matrix = np.full((size, size), np.nan, dtype=float)
    columns: list[int] = []
    for line in content.splitlines():
        if column_match := g16_log_patterns.NMR_COUPLING_COLUMN_HEADER.match(line):
            columns = [int(token) - 1 for token in column_match.group("columns").split()]
            continue
        row_match = g16_log_patterns.NMR_COUPLING_ROW.match(line)
        if row_match is None or not columns:
            continue
        row_index = int(row_match.group("row")) - 1
        values = [
            float(token.replace("D", "E").replace("d", "e"))
            for token in row_match.group("values").split()
        ]
        for column_index, value in zip(columns, values, strict=False):
            if 0 <= row_index < size and 0 <= column_index < size:
                matrix[row_index, column_index] = value
                matrix[column_index, row_index] = value
    return matrix if np.isfinite(matrix).all() else None


def extract_populations_from_state(state: ParseState) -> dict[str, Any]:
    infos: dict[str, Any] = {}
    mo: dict[str, Any] = {}
    pops: dict[str, Any] = {}
    polars: dict[str, Any] = {}
    focus_content, next_cursor = _focus_from_state(POPULATION_ANALYSIS, state)
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
            if matches := pattern.find_matches(sub_focus_content):
                mo[key] = [sym for matched in matches for sym in matched.group("symbols").split()]
        mo.update(_extract_molecular_orbital_payload_from_text(focus_content))
        patterns_and_keys_2: list[tuple[MolOPPattern, str, str, str, str, str | None]] = [
            (
                g16_log_patterns.MULLIKEN_SPIN_DENSITY,
                "mulliken_spins",
                "spin",
                "mulliken",
                "spin_density",
                "total",
            ),
            (
                g16_log_patterns.MULLIKEN_POPULATION,
                "mulliken_charges",
                "charge",
                "mulliken",
                "charge",
                None,
            ),
            (
                g16_log_patterns.APT_POPULATION,
                "apt_charges",
                "charge",
                "apt",
                "charge",
                None,
            ),
            (
                g16_log_patterns.LOWDIN_POPULATION,
                "lowdin_charges",
                "charge",
                "lowdin",
                "charge",
                None,
            ),
        ]
        for pattern, key, group_name, scheme, quantity, spin_channel in patterns_and_keys_2:
            if matches := _last_pattern_matches(pattern, state.content):
                pops[key] = {
                    "scheme": scheme,
                    "quantity": quantity,
                    "values": [float(matched.group(group_name)) for matched in matches],
                    "spin_channel": spin_channel,
                    "source_label": pattern.description,
                }
        sub_focus_content, focus_content = g16_log_patterns.ELECTRONIC_SPATIAL_EXTENT.split_content(
            focus_content
        )
        if matches := g16_log_patterns.ELECTRONIC_SPATIAL_EXTENT.find_matches(sub_focus_content):
            polars["electronic_spatial_extent"] = (
                float(matches[0].group("value")) * atom_ureg.bohr**2
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
            (g16_log_patterns.OCTAPOLE_MOMENT, "octapole", atom_ureg.debye * atom_ureg.angstrom**2),
            (
                g16_log_patterns.HEXADECAPOLE_MOMENT,
                "hexadecapole",
                atom_ureg.debye * atom_ureg.angstrom**3,
            ),
        ]
        for pattern, key, unit in patterns_and_keys_3:
            sub_focus_content, focus_content = pattern.split_content(focus_content)
            if matches := pattern.find_matches(sub_focus_content):
                polars[key] = (
                    np.array([float(matched.group("value")) for matched in matches]) * unit
                )
        if exact_polarizability := _extract_labeled_float_tokens(
            remainder_content, "Exact polarizability:", expected_count=6, decimal_places=3
        ):
            polars["polarizability_tensor"] = np.array(exact_polarizability) * atom_ureg.bohr**3
        elif approx_polarizability := _extract_labeled_float_tokens(
            remainder_content, "Approx polarizability:", expected_count=6, decimal_places=3
        ):
            polars["polarizability_tensor"] = np.array(approx_polarizability) * atom_ureg.bohr**3
        if matches := _last_pattern_matches(g16_log_patterns.HIRSHFELD_POPULATION, state.content):
            hirshfeld_source = "Hirshfeld charges, spin densities, dipoles, and CM5 charges"
            pops["hirshfeld_charges"] = {
                "scheme": "hirshfeld",
                "quantity": "charge",
                "values": [float(matched.group("charge")) for matched in matches],
                "source_label": hirshfeld_source,
            }
            pops["hirshfeld_spins"] = {
                "scheme": "hirshfeld",
                "quantity": "spin_density",
                "values": [float(matched.group("spin")) for matched in matches],
                "spin_channel": "total",
                "source_label": hirshfeld_source,
            }
            pops["cm5_charges"] = {
                "scheme": "cm5",
                "quantity": "charge",
                "values": [float(matched.group("q_cm5")) for matched in matches],
                "source_label": hirshfeld_source,
            }
        if matches := _last_pattern_matches(g16_log_patterns.NPA_POPULATION, state.content):
            pops["npa_charges"] = {
                "scheme": "npa",
                "quantity": "charge",
                "values": [float(matched.group("charge")) for matched in matches],
                "source_label": "Summary of Natural Population Analysis",
            }
        if matches := _last_pattern_matches(g16_log_patterns.ESP_POPULATION, state.content):
            pops["esp_charges"] = {
                "scheme": "esp",
                "quantity": "charge",
                "values": [float(matched.group("charge")) for matched in matches],
                "source_label": "ESP charges",
            }
        if dipole_before_force := _extract_labeled_float_tokens(
            remainder_content, "Dipole        =", expected_count=3, decimal_places=8
        ):
            polars["dipole"] = (
                np.array(dipole_before_force)
                * atom_ureg.atomic_unit_of_current
                * atom_ureg.atomic_unit_of_time
                * atom_ureg.bohr
            )
        if polarizability_before_force := _extract_labeled_float_tokens(
            remainder_content, "Polarizability=", expected_count=6, decimal_places=8
        ):
            polars["polarizability_tensor"] = (
                np.array(polarizability_before_force) * atom_ureg.bohr**3
            )
        if mo:
            mo = _trim_molecular_orbital_symmetries(mo)
            infos["molecular_orbitals"] = mo
        if pops:
            infos["charge_spin_populations"] = {"populations": pops}
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
        focus_content, next_cursor = _focus_from_state(FREQUENCY_ANALYSIS, state)
    else:
        focus_content = block[start_index:end_index]
        next_cursor = end_index + 1
    if focus_content == "":
        return None
    state.advance_to(next_cursor)
    vib_dict: dict[str, Any] = {}
    length = 0
    if matches := g16_log_patterns.FREQUENCIES.find_matches(focus_content):
        vib_dict["frequencies"] = (
            np.array(
                [
                    value
                    for matched in matches
                    for value in _extract_float_tokens(matched.group("values"))
                ]
            )
            * atom_ureg.cm_1
        )
        length = len(vib_dict["frequencies"])
    if matches := g16_log_patterns.FREQUENCIES_REDUCED_MASS.find_matches(focus_content):
        vib_dict["reduced_masses"] = (
            np.array(
                [
                    value
                    for matched in matches
                    for value in _extract_float_tokens(matched.group("values"))
                ]
            )
            * atom_ureg.amu
        )
    if matches := g16_log_patterns.FREQUENCIES_FORCE_CONSTANTS.find_matches(focus_content):
        vib_dict["force_constants"] = (
            np.array(
                [
                    value
                    for matched in matches
                    for value in _extract_float_tokens(matched.group("values"))
                ]
            )
            * atom_ureg.mdyne
            / atom_ureg.angstrom
        )
    if matches := g16_log_patterns.FREQUENCIES_IR_INTENSITIES.find_matches(focus_content):
        vib_dict["IR_intensities"] = (
            np.array(
                [
                    value
                    for matched in matches
                    for value in _extract_float_tokens(matched.group("values"))
                ]
            )
            * atom_ureg.km
            / atom_ureg.mol
        )
    if matches := g16_log_patterns.FREQUENCIES_MODE.find_matches(focus_content):
        v = np.array(
            [
                value
                for matched in matches
                for value in _extract_float_tokens(matched.group("values"))
            ]
        ).reshape(-1, 3)
        L = len(v) // length
        v1, v2, v3 = v[0::3], v[1::3], v[2::3]
        vib_dict["vibration_modes"] = [
            vn[i * L : i * L + L] for i in range(length // 3) for vn in [v1, v2, v3]
        ] * atom_ureg.angstrom
        vib_dict.update(
            {
                "axis_order": ("mode", "atom", "cartesian"),
                "atom_order": "source",
                "normalization": "unknown",
                "mass_weighting": "unknown",
            }
        )
    return vib_dict or None


def extract_thermal_infos_from_state(state: ParseState) -> dict[str, Any] | None:
    source_content = state.remaining_content
    focus_content, next_cursor = _focus_from_state(THERMOCHEMISTRY_PART, state)
    if focus_content == "":
        return None
    state.advance_to(next_cursor)
    thermal_dict: dict[str, Any] = {}
    if matches := g16_log_patterns.MOLECULAR_MASS.find_matches(source_content):
        thermal_dict["molecular_mass"] = float(matches[0].group("mass")) * atom_ureg.amu
    if matches := g16_log_patterns.MOMENTS_OF_INERTIA.find_matches(source_content):
        values = _extract_float_tokens(matches[0].group("values"))
        if values:
            thermal_dict["moments_of_inertia"] = (
                np.array(values) * atom_ureg.amu * atom_ureg.bohr**2
            )
    if matches := g16_log_patterns.ROTATIONAL_SYMMETRY_NUMBER.find_matches(source_content):
        thermal_dict["rotational_symmetry_number"] = int(matches[0].group("number"))
    if matches := g16_log_patterns.ROTATIONAL_TEMPERATURE.find_matches(source_content):
        thermal_dict["rotational_temperatures"] = (
            np.array([float(matches[0].group(axis)) for axis in ("a", "b", "c")]) * atom_ureg.K
        )
    if matches := g16_log_patterns.ROTATIONAL_CONST_IN_FREQUENCY_ANALYSIS.find_matches(
        source_content
    ):
        thermal_dict["rotational_constants"] = (
            np.array(
                [parse_rotational_constant(matches[0].group(axis)) for axis in ("a", "b", "c")]
            )
            * atom_ureg.gigahertz
        )
    if matches := g16_log_patterns.VIBRATIONAL_TEMPERATURE.find_matches(source_content):
        numeric_tokens = [
            value for matched in matches for value in _extract_float_tokens(matched.group("values"))
        ]
        if numeric_tokens:
            thermal_dict["vibrational_temperatures"] = np.array(numeric_tokens) * atom_ureg.K
            frequencies = [
                value
                for matched in g16_log_patterns.FREQUENCIES.find_matches(source_content)
                for value in _extract_float_tokens(matched.group("values"))
            ]
            positive_mode_indices = [
                mode_index for mode_index, frequency in enumerate(frequencies) if frequency > 0
            ]
            if len(positive_mode_indices) == len(numeric_tokens):
                thermal_dict["vibrational_temperature_mode_indices"] = positive_mode_indices
            elif len(frequencies) == len(numeric_tokens):
                thermal_dict["vibrational_temperature_mode_indices"] = list(range(len(frequencies)))
    if matches := g16_log_patterns.THERMOCHEMISTRY_CORRECTION.find_matches(focus_content):
        correction_mapping: dict[tuple[str, str], str] = {
            ("Zero-point", ""): "ZPVE",
            ("Thermal", " to Energy"): "TCE",
            ("Thermal", " to Enthalpy"): "TCH",
            ("Thermal", " to Gibbs Free Energy"): "TCG",
        }
        for matched in matches:
            thermal_dict[correction_mapping[(matched.group("kind"), matched.group("suffix"))]] = (
                float(matched.group("value")) * atom_ureg.Unit("hartree/particle")
            )
    if matches := g16_log_patterns.THERMOCHEMISTRY_SUM.find_matches(focus_content):
        summary_mapping: dict[str, str] = {
            "zero-point Energies": "U_0",
            "thermal Energies": "U_T",
            "thermal Enthalpies": "H_T",
            "thermal Free Energies": "G_T",
        }
        for matched in matches:
            thermal_dict[summary_mapping[matched.group("term")]] = float(
                matched.group("value")
            ) * atom_ureg.Unit("hartree/particle")
    if matches := g16_log_patterns.THERMOCHEMISTRY_CV_S.find_matches(focus_content):
        thermal_dict["C_V"] = float(matches[0].group("cv")) * atom_ureg.Unit("cal/mol/K")
        thermal_dict["S"] = float(matches[0].group("entropy")) * atom_ureg.Unit("cal/mol/K")
    return thermal_dict or None


def extract_forces_from_state(state: ParseState) -> NumpyQuantity | None:
    focus_content, next_cursor = _focus_from_state(FORCES_IN_CARTESIAN, state)
    if focus_content == "":
        return None
    state.advance_to(next_cursor)
    if matches := g16_log_patterns.FORCES_IN_CARTESIAN.find_matches(focus_content):
        return (
            np.array(
                [
                    [
                        float(matched.group("x")),
                        float(matched.group("y")),
                        float(matched.group("z")),
                    ]
                    for matched in matches
                ]
            )
            * atom_ureg.hartree
            / atom_ureg.bohr
        )
    return None


def extract_hessian_from_state(state: ParseState) -> NumpyQuantity | None:
    focus_content, next_cursor = _focus_from_state(HESSIAN_IN_CARTESIAN, state)
    if focus_content == "":
        return None
    state.advance_to(next_cursor)
    if matches := g16_log_patterns.HESSIAN_IN_CARTESIAN.find_matches(focus_content):
        hessian_dict: dict[int, list[float]] = {}
        for matched in matches:
            row = int(matched.group("row"))
            elements = list(map(float, matched.group("values").replace("D", "E").split()))
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


def extract_berny_from_state(
    state: ParseState,
    *,
    capture_source_evidence: bool = False,
) -> dict[str, Any] | None:
    berny_dict: dict[str, Any] = {}
    source_converged: dict[str, bool | None] = {}
    source_labels: dict[str, str] = {}
    start_index = state.content.find(
        "Item               Value     Threshold  Converged?", state.cursor
    )
    end_index = state.content.find(
        "GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad",
        start_index,
    )
    if start_index == -1 or end_index == -1:
        focus_content, next_cursor = _focus_from_state(BERNY_STATE_MAJOR_PART, state)
        if focus_content == "":
            focus_content, next_cursor = _focus_from_state(BERNY_STATE_BACKUP_PART, state)
    else:
        focus_content = state.content[start_index:end_index]
        next_cursor = end_index
    if focus_content == "":
        return None
    state.advance_to(next_cursor)
    if matches := g16_log_patterns.BERNY_STATE.find_matches(focus_content):
        mapping = {
            "Maximum Force": "max_force",
            "RMS     Force": "rms_force",
            "Maximum Displacement": "max_displacement",
            "RMS     Displacement": "rms_displacement",
        }
        units = {
            "max_force": atom_ureg.hartree / atom_ureg.bohr,
            "rms_force": atom_ureg.hartree / atom_ureg.bohr,
            "max_displacement": atom_ureg.bohr,
            "rms_displacement": atom_ureg.bohr,
        }
        for matched in matches:
            key = mapping[matched.group("label")]
            berny_dict[key] = float(matched.group("value")) * units[key]
            berny_dict[f"{key}_threshold"] = float(matched.group("threshold")) * units[key]
            if capture_source_evidence:
                source_converged[key] = matched.group("converged") == "YES"
                source_labels[key] = " ".join(matched.group("label").split())
    if matches := g16_log_patterns.ENERGY_CHANGE.find_matches(focus_content):
        berny_dict["energy_change"] = (
            float(matches[0].group("value").replace("D", "E")) * atom_ureg.hartree
        )
        if capture_source_evidence:
            source_converged["energy_change"] = None
            source_labels["energy_change"] = "Predicted change in Energy"
    berny_dict["geometry_optimized"] = bool(
        g16_log_patterns.BERNY_CONCLUSION.find_matches(focus_content)
    )
    if source_converged:
        berny_dict["source_converged"] = source_converged
    if source_labels:
        berny_dict["source_labels"] = source_labels
    return berny_dict or None


def extract_electric_dipole_and_polarizability_from_state(
    state: ParseState,
) -> dict[str, Any] | None:
    focus_content, next_cursor = _focus_from_state(ELECTRIC_DIPOLE_PART, state)
    if focus_content == "":
        return None
    state.advance_to(next_cursor)
    polarizability_dict: dict[str, Any] = {}
    if matches := g16_log_patterns.ELECTRIC_DIPOLE_MOMENT.find_matches(focus_content):
        polarizability_dict["dipole"] = (
            np.array(
                [float(matched.group("eigenvalue").replace("D", "E")) for matched in matches[1:]]
            )
            * atom_ureg.debye
        )
    if matches := g16_log_patterns.DIPOLE_POLARIZABILITY.find_matches(focus_content):
        polarizability_dict["isotropic_polarizability"] = (
            float(matches[0].group("frequency_0").replace("D", "E")) * atom_ureg.bohr**3
        )
        polarizability_dict["anisotropic_polarizability"] = (
            float(matches[1].group("frequency_0").replace("D", "E")) * atom_ureg.bohr**3
        )
        polarizability_dict["polarizability_tensor"] = (
            np.array(
                [float(matched.group("frequency_0").replace("D", "E")) for matched in matches[2:]]
            )
            * atom_ureg.bohr**3
        )
    return polarizability_dict or None


def extract_tail_metadata_from_state(state: ParseState) -> dict[str, Any]:
    focus_content, next_cursor = _focus_from_state(ARCHIVE_TAIL, state)
    if focus_content == "":
        return {}
    state.advance_to(next_cursor)
    payload, _tail_remaining = parse_archive_tail_payload(focus_content, include_coords=True)
    return payload.get("metadata", {})


def _energy_observations(payload: Mapping[str, Any]) -> tuple[EnergyObservation, ...]:
    raw_observations = payload.get("observations")
    if raw_observations is None:
        return ()
    if not isinstance(raw_observations, Sequence) or isinstance(raw_observations, (str, bytes)):
        raise TypeError("energy observations must be a sequence")
    return tuple(
        observation
        if isinstance(observation, EnergyObservation)
        else EnergyObservation.model_validate(observation)
        for observation in raw_observations
    )


def _energy_observation_identity(
    observation: EnergyObservation,
) -> tuple[str, str, str, str]:
    """Return the complete observation identity in a canonical energy unit."""

    magnitude = float(observation.value.to("hartree").magnitude)
    return (
        observation.method,
        observation.quantity_semantics,
        observation.source_label,
        magnitude.hex(),
    )


def merge_g16_energy_payloads(
    primary: Mapping[str, Any] | None,
    fallback: Mapping[str, Any] | None,
) -> dict[str, Any]:
    """Merge G16 energies while retaining all distinct typed observations.

    Scalar fields keep the first non-None value, so live frame extraction stays
    authoritative and the archive tail only fills missing values. Observations
    preserve primary-then-fallback order and deduplicate by their full contract
    identity after normalizing values to hartree.
    """

    payloads = tuple(payload for payload in (primary, fallback) if payload is not None)
    merged: dict[str, Any] = {}
    for payload in payloads:
        for key, value in payload.items():
            if key == "observations" or value is None or key in merged:
                continue
            merged[key] = value

    observations: list[EnergyObservation] = []
    seen_observations: set[tuple[str, str, str, str]] = set()
    for payload in payloads:
        for observation in _energy_observations(payload):
            identity = _energy_observation_identity(observation)
            if identity in seen_observations:
                continue
            seen_observations.add(identity)
            observations.append(observation)
    if observations:
        merged["observations"] = observations
    return merged


def _archive_energy_observations(energies: dict[str, Any]) -> list[EnergyObservation]:
    method_by_field = {
        "electronic_energy": "electronic",
        "reference_energy": "reference",
        "mp2_energy": "MP2",
        "mp3_energy": "MP3",
        "mp4_energy": "MP4",
        "mp5_energy": "MP5",
        "ccsd_energy": "CCSD",
        "ccsd_t_energy": "CCSD(T)",
    }
    return [
        EnergyObservation(
            method=method,
            quantity_semantics="total_energy",
            value=energies[field_name],
            source_label=f"Gaussian archive {method}",
        )
        for field_name, method in method_by_field.items()
        if energies.get(field_name) is not None
    ]


def extract_archive_tail_payload_from_state(
    state: ParseState,
    *,
    capture_source_evidence: bool = False,
) -> dict[str, Any]:
    focus_content, next_cursor = _focus_from_state(ARCHIVE_TAIL, state)
    if focus_content == "":
        return {}

    # All archive sub-payloads live inside the same bounded archive block. Parse
    # them before advancing past the block, otherwise energy fallback data can be
    # skipped on terminal frequency frames without a live "SCF Done" line.
    state.advance_to(next_cursor)
    payload, _tail_remaining = parse_archive_tail_payload(focus_content, include_coords=True)
    energies = payload.get("energies")
    if capture_source_evidence and isinstance(energies, dict):
        observations = _archive_energy_observations(energies)
        if observations:
            energies["observations"] = observations
    return payload


def extract_tail_energies_from_state(state: ParseState) -> dict[str, Any] | None:
    focus_content, next_cursor = _focus_from_state(ENERGIES_IN_ARCHIVE_TAIL, state)
    if focus_content == "":
        return None
    state.advance_to(next_cursor)
    return parse_archive_tail_energies(focus_content)


def extract_tail_thermal_infos_from_state(state: ParseState) -> dict[str, Any] | None:
    focus_content, next_cursor = _focus_from_state(THERMOCHEMISTRY_IN_ARCHIVE_TAIL, state)
    if focus_content == "":
        return None
    state.advance_to(next_cursor)
    return parse_archive_tail_thermal_infos(focus_content)


def extract_tail_polarizability_from_state(state: ParseState) -> dict[str, Any] | None:
    return parse_archive_tail_polarizability(state.remaining_content)


def extract_tail_hessian_from_state(state: ParseState) -> NumpyQuantity | None:
    focus_content, next_cursor = _focus_from_state(HESSIAN_IN_ARCHIVE_TAIL, state)
    if focus_content == "":
        return None
    state.advance_to(next_cursor)
    return parse_archive_tail_hessian(focus_content)


def extract_rotation_consts_from_state(state: ParseState) -> NumpyQuantity | None:
    return extract_rotation_constants(state.remaining_content)
