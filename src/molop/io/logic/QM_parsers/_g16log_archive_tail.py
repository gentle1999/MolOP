from __future__ import annotations

from typing import Any

import numpy as np
from pint.facets.numpy.quantity import NumpyQuantity
from rdkit import Chem

from molop.config import moloplogger
from molop.io.base_models.SearchPattern import MolOPPattern
from molop.io.patterns.G16Patterns import g16_log_patterns
from molop.unit import atom_ureg
from molop.utils.functions import fill_symmetric_matrix


pt = Chem.GetPeriodicTable()


def extract_archive_tail_block(content: str) -> tuple[str | None, str]:
    located = g16_log_patterns.ARCHIVE_TAIL.locate_content(content)
    if located is None:
        return None, content
    start_start, _start_end, end_start, end_end = located
    return content[start_start:end_start], content[end_end:]


def parse_archive_tail(content: str, *, include_coords: bool = False) -> tuple[dict[str, Any], str]:
    raw_tail_text, remaining_content = extract_archive_tail_block(content)
    if raw_tail_text is None:
        return {}, content
    payload, _remaining_tail = parse_archive_tail_payload(
        raw_tail_text,
        include_coords=include_coords,
        structure_only=True,
    )
    return payload.get("metadata", {}), remaining_content


def extract_archive_tail_payload(
    content: str,
    *,
    include_coords: bool = False,
    structure_only: bool = False,
) -> tuple[dict[str, Any], str]:
    raw_tail_text, remaining_content = extract_archive_tail_block(content)
    if raw_tail_text is None:
        return {}, content
    payload, _remaining_tail = parse_archive_tail_payload(
        raw_tail_text,
        include_coords=include_coords,
        structure_only=structure_only,
    )
    return payload, remaining_content


def parse_archive_tail_for_frame(
    content: str, *, include_coords: bool
) -> tuple[dict[str, Any], str]:
    raw_tail_text, _remaining_content = extract_archive_tail_block(content)
    if raw_tail_text is None:
        return {}, content
    return parse_archive_tail_metadata(raw_tail_text, include_coords=include_coords)


def parse_archive_tail_energies(content: str) -> dict[str, Any] | None:
    energy_dict: dict[str, Any] = {}
    normalized_content = content.replace("\n ", "")
    if matches := g16_log_patterns.ENERGIES_IN_ARCHIVE_TAIL.get_matches(normalized_content):
        for match in matches:
            energy_name = match[0]
            energies_value = float(match[1]) * atom_ureg.hartree
            if "HF" in energy_name:
                energy_dict["reference_energy"] = energies_value
            if "MP2" in energy_name:
                energy_dict["mp2_energy"] = energies_value
            if "MP3" in energy_name:
                energy_dict["mp3_energy"] = energies_value
            if "MP4" in energy_name:
                energy_dict["mp4_energy"] = energies_value
            if "CCSD" in energy_name:
                energy_dict["ccsd_energy"] = energies_value
    return energy_dict or None


def parse_archive_tail_thermal_infos(content: str) -> dict[str, Any] | None:
    thermal_dict: dict[str, Any] = {}
    thermal_mapping = {
        "ZeroPoint": "ZPVE",
        "Thermal": "TCE",
        "ETot": "U_T",
        "HTot": "H_T",
        "GTot": "G_T",
    }
    normalized_content = content.replace("\n ", "")
    if matches := g16_log_patterns.THERMOCHEMISTRY_IN_ARCHIVE_TAIL.get_matches(normalized_content):
        for match in matches:
            if match[0] in thermal_mapping:
                thermal_dict[thermal_mapping[match[0]]] = float(match[1]) * atom_ureg.Unit(
                    "hartree/particle"
                )
    return thermal_dict or None


def parse_archive_tail_polarizability(content: str) -> dict[str, Any] | None:
    polarizability_dict: dict[str, Any] = {}
    normalized_content = content.replace("\n ", "")
    if matches := g16_log_patterns.DIPOLE_IN_ARCHIVE_TAIL.get_matches(normalized_content):
        polarizability_dict["dipole"] = (
            np.array([float(value) for value in matches[0]]) * atom_ureg.debye
        )
    if matches := g16_log_patterns.POLAR_IN_ARCHIVE_TAIL.get_matches(normalized_content):
        polarizability_dict["polarizability_tensor"] = (
            np.array([float(value) for value in matches[0]]) * atom_ureg.bohr**3
        )
    if matches := g16_log_patterns.QUADRUPOLE_IN_ARCHIVE_TAIL.get_matches(normalized_content):
        polarizability_dict["quadrupole"] = (
            np.array([float(value) for value in matches[0]]) * atom_ureg.debye * atom_ureg.angstrom
        )
    return polarizability_dict or None


def parse_archive_tail_hessian(content: str) -> NumpyQuantity | None:
    normalized_content = content.replace("\n ", "")
    focus_content, _remaining_content = g16_log_patterns.HESSIAN_IN_ARCHIVE_TAIL.split_content(
        normalized_content
    )
    if matches := g16_log_patterns.HESSIAN_IN_ARCHIVE_TAIL.get_matches(focus_content):
        try:
            return (
                fill_symmetric_matrix(np.array([float(match[0]) for match in matches]))
                * atom_ureg.hartree
                / atom_ureg.bohr**2
            )
        except (AssertionError, ValueError) as exc:
            moloplogger.warning(
                "Skipping invalid Gaussian archive hessian: %s | values=%d",
                exc,
                len(matches),
            )
    return None


def parse_archive_tail_payload(
    raw_tail_text: str,
    *,
    include_coords: bool = False,
    structure_only: bool = False,
) -> tuple[dict[str, Any], str]:
    metadata, remaining_tail = parse_archive_tail_metadata(
        raw_tail_text,
        include_coords=include_coords,
    )
    payload: dict[str, Any] = {"metadata": metadata}
    if structure_only:
        return payload, remaining_tail

    if energies := parse_archive_tail_energies(raw_tail_text):
        payload["energies"] = energies
    if thermal_infos := parse_archive_tail_thermal_infos(raw_tail_text):
        payload["thermal_informations"] = thermal_infos
    if polarizability := parse_archive_tail_polarizability(raw_tail_text):
        payload["polarizability"] = polarizability
    if (hessian := parse_archive_tail_hessian(raw_tail_text)) is not None:
        payload["hessian"] = hessian
    return payload, remaining_tail


def _parse_and_update_scalar(
    focus_content: str,
    tail_dict: dict[str, Any],
    *,
    pattern: MolOPPattern,
    key: str,
) -> str:
    try:
        sub_focus_content, sub_continued_content = pattern.split_content(focus_content)
        if matches := pattern.get_matches(sub_focus_content):
            tail_dict[key] = matches[0][0]
            return sub_continued_content
    except Exception as exc:
        moloplogger.error(f"Error in parsing {key}: {exc}")
    return focus_content


def parse_archive_tail_metadata(
    raw_tail_text: str,
    *,
    include_coords: bool = False,
) -> tuple[dict[str, Any], str]:
    tail_dict: dict[str, Any] = {}
    focus_content = raw_tail_text.replace("\n ", "")

    focus_content = _parse_and_update_scalar(
        focus_content,
        tail_dict,
        pattern=g16_log_patterns.JOB_TYPE_IN_ARCHIVE_TAIL,
        key="job_type",
    )
    focus_content = _parse_and_update_scalar(
        focus_content,
        tail_dict,
        pattern=g16_log_patterns.FUNCTIONAL_IN_ARCHIVE_TAIL,
        key="functional",
    )
    focus_content = _parse_and_update_scalar(
        focus_content,
        tail_dict,
        pattern=g16_log_patterns.BASIS_SET_IN_ARCHIVE_TAIL,
        key="basis_set",
    )
    focus_content = _parse_and_update_scalar(
        focus_content,
        tail_dict,
        pattern=g16_log_patterns.KEYWORDS_IN_ARCHIVE_TAIL,
        key="keywords",
    )
    focus_content = _parse_and_update_scalar(
        focus_content,
        tail_dict,
        pattern=g16_log_patterns.TITLE_IN_ARCHIVE_TAIL,
        key="title_card",
    )

    sub_focus_content, sub_continued_content = (
        g16_log_patterns.CHARGE_SPIN_MULTIPLICITY_IN_ARCHIVE_TAIL.split_content(focus_content)
    )
    if matches := g16_log_patterns.CHARGE_SPIN_MULTIPLICITY_IN_ARCHIVE_TAIL.get_matches(
        sub_focus_content
    ):
        charge_multiplicity = matches[0]
        tail_dict["charge"] = int(charge_multiplicity[0])
        tail_dict["multiplicity"] = int(charge_multiplicity[1])
        focus_content = sub_continued_content

    if include_coords:
        sub_focus_content, sub_continued_content = (
            g16_log_patterns.COORS_IN_ARCHIVE_TAIL.split_content(focus_content)
        )
        if matches := g16_log_patterns.COORS_IN_ARCHIVE_TAIL.get_matches(sub_focus_content):
            tail_dict["atoms"] = [pt.GetAtomicNumber(match[0]) for match in matches]
            tail_dict["coords"] = (
                np.array([list(map(float, match[1:])) for match in matches]) * atom_ureg.angstrom
            )
            focus_content = sub_continued_content

    focus_content = _parse_and_update_scalar(
        focus_content,
        tail_dict,
        pattern=g16_log_patterns.VERSION_IN_ARCHIVE_TAIL,
        key="qm_software_version",
    )

    moloplogger.debug(f"parsed tail_dict: {tail_dict}")
    return tail_dict, focus_content
