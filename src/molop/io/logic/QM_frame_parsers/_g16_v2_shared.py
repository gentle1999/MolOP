from __future__ import annotations

from typing import Any

import numpy as np
from pint.facets.numpy.quantity import NumpyQuantity

from molop.io.base_models.SearchPattern import MolOPPatternV2
from molop.io.patterns.G16Patterns import g16_log_patterns
from molop.unit import atom_ureg


INPUT_COORDS_V2 = MolOPPatternV2.from_pattern(g16_log_patterns.INPUT_COORDS)
STANDARD_COORDS_V2 = MolOPPatternV2.from_pattern(g16_log_patterns.STANDARD_COORDS)
SCF_ENERGIES_V2 = MolOPPatternV2.from_pattern(g16_log_patterns.SCF_ENERGIES)
ISOTROPIC_POLARIZABILITY_V2 = MolOPPatternV2.from_pattern(g16_log_patterns.ISOTROPIC_POLARIZABILITY)
POPULATION_ANALYSIS_V2 = MolOPPatternV2.from_pattern(g16_log_patterns.POPULATION_ANALYSIS)
FREQUENCY_ANALYSIS_V2 = MolOPPatternV2.from_pattern(g16_log_patterns.FREQUENCY_ANALYSIS)
THERMOCHEMISTRY_PART_V2 = MolOPPatternV2.from_pattern(g16_log_patterns.THERMOCHEMISTRY_PART)
FORCES_IN_CARTESIAN_V2 = MolOPPatternV2.from_pattern(g16_log_patterns.FORCES_IN_CARTESIAN)
HESSIAN_IN_CARTESIAN_V2 = MolOPPatternV2.from_pattern(g16_log_patterns.HESSIAN_IN_CARTESIAN)
BERNY_STATE_MAJOR_PART_V2 = MolOPPatternV2.from_pattern(g16_log_patterns.BERNY_STATE_MAJOR_PART)
BERNY_STATE_BACKUP_PART_V2 = MolOPPatternV2.from_pattern(g16_log_patterns.BERNY_STATE_BACKUP_PART)
ELECTRIC_DIPOLE_PART_V2 = MolOPPatternV2.from_pattern(g16_log_patterns.ELECTRIC_DIPOLE_PART)
ENERGIES_IN_ARCHIVE_TAIL_V2 = MolOPPatternV2.from_pattern(g16_log_patterns.ENERGIES_IN_ARCHIVE_TAIL)
THERMOCHEMISTRY_IN_ARCHIVE_TAIL_V2 = MolOPPatternV2.from_pattern(
    g16_log_patterns.THERMOCHEMISTRY_IN_ARCHIVE_TAIL
)
HESSIAN_IN_ARCHIVE_TAIL_V2 = MolOPPatternV2.from_pattern(g16_log_patterns.HESSIAN_IN_ARCHIVE_TAIL)
ARCHIVE_TAIL_V2 = MolOPPatternV2.from_pattern(g16_log_patterns.ARCHIVE_TAIL)


def extract_coords(
    coords_match: list[tuple[str | Any, ...]],
) -> tuple[list[int], NumpyQuantity] | tuple[None, None]:
    if len(coords_match) == 0:
        return None, None
    atoms: list[int] = []
    coords: list[list[float]] = []
    for match in coords_match:
        atoms.append(int(match[0]))
        coords.append([float(match[1]), float(match[2]), float(match[3])])
    return atoms, np.array(coords) * atom_ureg.angstrom


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

        if stripped.startswith("Alpha "):
            orbital_type = "Alpha"
            remainder = stripped.removeprefix("Alpha ")
        elif stripped.startswith("Beta "):
            orbital_type = "Beta"
            remainder = stripped.removeprefix("Beta ")
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

        values = [float(value.replace("D", "E").replace("d", "E")) for value in energies.split()]
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


def _summarize_parse_context(text: str, *, limit: int = 240) -> str:
    normalized = " ".join(text.split())
    return normalized[:limit] + ("..." if len(normalized) > limit else "")
