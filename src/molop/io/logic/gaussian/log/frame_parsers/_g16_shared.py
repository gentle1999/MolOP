from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import numpy as np
from pint.facets.numpy.quantity import NumpyQuantity
from pint.facets.plain import PlainQuantity

from molop.io.base_models.SearchPattern import MolOPMatch
from molop.io.logic.gaussian.log.parsers._g16_log_patterns import g16_log_patterns
from molop.unit import atom_ureg


INPUT_COORDS = g16_log_patterns.INPUT_COORDS
STANDARD_COORDS = g16_log_patterns.STANDARD_COORDS
SCF_ENERGIES = g16_log_patterns.SCF_ENERGIES
ISOTROPIC_POLARIZABILITY = g16_log_patterns.ISOTROPIC_POLARIZABILITY
POPULATION_ANALYSIS = g16_log_patterns.POPULATION_ANALYSIS
FREQUENCY_ANALYSIS = g16_log_patterns.FREQUENCY_ANALYSIS
THERMOCHEMISTRY_PART = g16_log_patterns.THERMOCHEMISTRY_PART
FORCES_IN_CARTESIAN = g16_log_patterns.FORCES_IN_CARTESIAN
HESSIAN_IN_CARTESIAN = g16_log_patterns.HESSIAN_IN_CARTESIAN
BERNY_STATE_MAJOR_PART = g16_log_patterns.BERNY_STATE_MAJOR_PART
BERNY_STATE_BACKUP_PART = g16_log_patterns.BERNY_STATE_BACKUP_PART
ELECTRIC_DIPOLE_PART = g16_log_patterns.ELECTRIC_DIPOLE_PART
ENERGIES_IN_ARCHIVE_TAIL = g16_log_patterns.ENERGIES_IN_ARCHIVE_TAIL
THERMOCHEMISTRY_IN_ARCHIVE_TAIL = g16_log_patterns.THERMOCHEMISTRY_IN_ARCHIVE_TAIL
HESSIAN_IN_ARCHIVE_TAIL = g16_log_patterns.HESSIAN_IN_ARCHIVE_TAIL
ARCHIVE_TAIL = g16_log_patterns.ARCHIVE_TAIL


def extract_coords(
    coords_match: Sequence[MolOPMatch],
) -> tuple[list[int], NumpyQuantity] | tuple[None, None]:
    if len(coords_match) == 0:
        return None, None
    atoms: list[int] = []
    coords: list[list[float]] = []
    for matched in coords_match:
        atoms.append(int(matched.group("atomic_number")))
        coords.append(
            [
                float(matched.group("x")),
                float(matched.group("y")),
                float(matched.group("z")),
            ]
        )
    return atoms, np.array(coords) * atom_ureg.angstrom


def extract_rotation_constants(block: str) -> NumpyQuantity | None:
    if matched := g16_log_patterns.ROTATIONAL_CONST.search(block):
        return (
            np.array(
                [
                    float(matched.group("a")),
                    float(matched.group("b")),
                    float(matched.group("c")),
                ],
            )
            * atom_ureg.gigahertz
        )
    return None


def _extract_float_tokens(text: str, *, decimal_places: int | None = None) -> list[float]:
    pattern = g16_log_patterns.float_token(decimal_places)
    return [
        float(matched.group("value").replace("D", "E").replace("d", "E"))
        for matched in pattern.find_matches(text)
    ]


def _parse_orbital_line_values(energies: str) -> list[float]:
    precise_values = [
        float(matched.group("value").replace("D", "E").replace("d", "E"))
        for matched in g16_log_patterns.FLOAT_TOKEN_5DP.find_matches(energies)
    ]
    if precise_values:
        return precise_values
    return [
        float(matched.group("value").replace("D", "E").replace("d", "E"))
        for matched in g16_log_patterns.FLOAT_TOKEN.find_matches(energies)
    ]


def _parse_frequency_line_values(line: str) -> list[float]:
    if "--" in line:
        line = line.split("--", 1)[1]
    return _parse_orbital_line_values(line)


def _extract_labeled_float_tokens(
    text: str,
    label: str,
    *,
    expected_count: int | None = None,
    decimal_places: int | None = None,
) -> list[float]:
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped.startswith(label):
            continue
        body = stripped.removeprefix(label)
        values = _extract_float_tokens(body)
        if expected_count is None or len(values) >= expected_count:
            return values if expected_count is None else values[:expected_count]
        if decimal_places is not None:
            precise_values = _extract_float_tokens(body, decimal_places=decimal_places)
            if expected_count is None or len(precise_values) >= expected_count:
                return precise_values if expected_count is None else precise_values[:expected_count]
    return []


def _extract_molecular_orbital_payload_from_text(focus_content: str) -> dict[str, Any]:
    mo: dict[str, Any] = {}
    temp_alpha_orbitals: list[float] = []
    temp_alpha_occupancies: list[bool] = []
    temp_beta_orbitals: list[float] = []
    temp_beta_occupancies: list[bool] = []

    for line in focus_content.splitlines():
        stripped = line.strip()
        if stripped.startswith("The electronic state is "):
            mo["electronic_state"] = stripped.removeprefix("The electronic state is ").removesuffix(
                "."
            )
            continue
        if "eigenvalues --" not in stripped:
            continue

        if stripped.startswith("Alpha"):
            orbital_type = "Alpha"
            remainder = stripped.removeprefix("Alpha").lstrip()
        elif stripped.startswith("Beta"):
            orbital_type = "Beta"
            remainder = stripped.removeprefix("Beta").lstrip()
        else:
            continue

        if remainder.startswith("occ. eigenvalues --"):
            occ_stat = "occ."
            energies = remainder.removeprefix("occ. eigenvalues --")
        elif remainder.startswith("virt. eigenvalues --"):
            occ_stat = "virt."
            energies = remainder.removeprefix("virt. eigenvalues --")
        else:
            continue

        values = _extract_float_tokens(energies, decimal_places=5)
        if orbital_type == "Alpha":
            temp_alpha_orbitals.extend(values)
            temp_alpha_occupancies.extend([occ_stat == "occ."] * len(values))
        else:
            temp_beta_orbitals.extend(values)
            temp_beta_occupancies.extend([occ_stat == "occ."] * len(values))

    if temp_alpha_orbitals:
        mo["alpha_energies"] = np.array(temp_alpha_orbitals) * atom_ureg.hartree
        mo["alpha_occupancies"] = temp_alpha_occupancies
    if temp_beta_orbitals:
        mo["beta_energies"] = np.array(temp_beta_orbitals) * atom_ureg.hartree
        mo["beta_occupancies"] = temp_beta_occupancies

    return mo


def _trim_molecular_orbital_symmetries(mo: dict[str, Any]) -> dict[str, Any]:
    alpha_energies = mo.get("alpha_energies")
    beta_energies = mo.get("beta_energies")
    alpha_symmetries = mo.get("alpha_symmetries")
    beta_symmetries = mo.get("beta_symmetries")

    if alpha_energies is not None and alpha_symmetries is not None:
        count = len(alpha_energies)
        if len(alpha_symmetries) > count:
            mo["alpha_symmetries"] = alpha_symmetries[-count:]
    if beta_energies is not None and beta_symmetries is not None:
        count = len(beta_energies)
        if len(beta_symmetries) > count:
            mo["beta_symmetries"] = beta_symmetries[-count:]
    return mo


def _parse_running_time(block: str) -> PlainQuantity | None:
    if matches := g16_log_patterns.PROCEDURE_TIME.find_matches(block):
        total_seconds = 0.0
        for matched in matches:
            total_seconds += float(matched.group("cpu")) + float(matched.group("elapsed"))
        return total_seconds * atom_ureg.second
    if matches := g16_log_patterns.JOB_TIME.find_matches(block):
        total_seconds = 0.0
        for matched in matches:
            total_seconds += (
                float(matched.group("days")) * 24 * 3600
                + float(matched.group("hours")) * 3600
                + float(matched.group("minutes")) * 60
                + float(matched.group("seconds"))
            )
        return total_seconds * atom_ureg.second
    return None


def _temperature_and_pressure_from_block(block: str) -> dict[str, Any]:
    if matches := g16_log_patterns.TEMPEREATURE_PRESSURE.find_content_matches(block):
        matched = matches[0]
        return {
            "temperature": float(matched.group("temperature")) * atom_ureg.K,
            "pressure": float(matched.group("pressure")) * atom_ureg.atm,
        }
    return {}


def _summarize_parse_context(text: str, *, limit: int = 240) -> str:
    normalized = " ".join(text.split())
    return normalized[:limit] + ("..." if len(normalized) > limit else "")
