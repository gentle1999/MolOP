from __future__ import annotations

from typing import Any

import numpy as np
from rdkit import Chem

from molop.config import moloplogger
from molop.io.base_models.SearchPattern import MolOPPattern
from molop.io.patterns.G16Patterns import g16_log_patterns
from molop.unit import atom_ureg


pt = Chem.GetPeriodicTable()


def extract_archive_tail_block(content: str) -> tuple[str | None, str]:
    located = g16_log_patterns.ARCHIVE_TAIL.locate_content(content)
    if located is None:
        return None, content
    start_start, _start_end, end_start, end_end = located
    return content[start_start:end_start], content[end_end:]


def parse_archive_tail(content: str, *, include_coords: bool) -> tuple[dict[str, Any], str]:
    raw_tail_text, remaining_content = extract_archive_tail_block(content)
    if raw_tail_text is None:
        return {}, content
    tail_dict, _remaining_tail = parse_archive_tail_metadata(
        raw_tail_text, include_coords=include_coords
    )
    return tail_dict, remaining_content


def parse_archive_tail_for_frame(
    content: str, *, include_coords: bool
) -> tuple[dict[str, Any], str]:
    raw_tail_text, _remaining_content = extract_archive_tail_block(content)
    if raw_tail_text is None:
        return {}, content
    return parse_archive_tail_metadata(raw_tail_text, include_coords=include_coords)


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
    include_coords: bool,
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
