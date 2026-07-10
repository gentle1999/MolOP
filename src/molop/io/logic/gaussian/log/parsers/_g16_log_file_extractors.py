from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import numpy as np
from pint.facets.plain import PlainQuantity

from molop.io.base_models.DataClasses import ImplicitSolvation, Status
from molop.io.base_models.ParseContainers import TextParseContext
from molop.io.base_models.SearchPattern import MolOPPattern
from molop.io.codec_exceptions import FormatMismatchError
from molop.io.logic.gaussian.log.parsers._g16_log_patterns import g16_log_patterns
from molop.unit import atom_ureg
from molop.utils.functions import find_rigid_transform


SPLIT_PATTERN = "Input orientation:"
SPLIT_PATTERN_2 = "Standard orientation:"
_GAUSSIAN_PROBE_BYTES = 20000
_GAUSSIAN_FINGERPRINTS = (
    "Entering Gaussian System",
    "Gaussian 16:",
    "This is part of the Gaussian(R)",
    "Gaussian, Inc.",
)


def ensure_g16_output_content(file_content: str) -> None:
    prefix = file_content[:_GAUSSIAN_PROBE_BYTES]
    if not any(fingerprint in prefix for fingerprint in _GAUSSIAN_FINGERPRINTS):
        raise FormatMismatchError("Not a Gaussian output file.")


def split_g16_sections(file_content: str) -> list[str]:
    matches = g16_log_patterns.LINK1_SECTION.find_matches(file_content)
    if not matches:
        return [file_content]
    sections: list[str] = []
    for idx, matched in enumerate(matches):
        next_start = matches[idx + 1].start() if idx + 1 < len(matches) else len(file_content)
        sections.append(file_content[matched.start() : next_start])
    return [section for section in sections if section.strip()]


def split_g16_section_frames(section_content: str) -> list[str]:
    split_ = SPLIT_PATTERN
    if SPLIT_PATTERN in section_content:
        split_ = SPLIT_PATTERN
    elif SPLIT_PATTERN_2 in section_content:
        split_ = SPLIT_PATTERN_2
    else:
        return []
    fragments = section_content.split(split_)
    header = fragments[0]
    frame_contents: list[str] = []
    for idx, fragment in enumerate(fragments[1:]):
        frame_block = f"{split_}\n{fragment}"
        if idx == 0 and header.strip():
            frame_block = f"{header}{frame_block}"
        frame_contents.append(frame_block)
    return frame_contents


def first_frame_value(frames: Sequence[Any], field: str) -> Any:
    for frame in frames:
        value = getattr(frame, field, None)
        if value is not None:
            return value
    return None


def last_frame_value(frames: Sequence[Any], field: str) -> Any:
    for frame in reversed(frames):
        value = getattr(frame, field, None)
        if value is not None:
            return value
    return None


def extract_g16_version(context: TextParseContext) -> str | None:
    focus_content = context.split(g16_log_patterns.VERSION)
    if matches := g16_log_patterns.VERSION.find_matches(focus_content):
        return matches[0].group("version")
    return None


def extract_g16_options(context: TextParseContext) -> str | None:
    focus_content, continued_content = g16_log_patterns.OPTIONS.split_content(context.content)
    if matches := g16_log_patterns.OPTIONS.find_matches(focus_content):
        options = "\n".join(
            [f"{matched.group('key')}={matched.group('value')}" for matched in matches]
        )
        context.content = continued_content
        return options
    return None


def extract_g16_keywords(context: TextParseContext) -> str | None:
    focus_content, continued_content = g16_log_patterns.KEYWORDS.split_content(context.content)
    if len(keyword_lines := focus_content.splitlines()) >= 3:
        context.content = continued_content
        return "\n".join(keyword_lines[1:-1]).replace("\n ", "")
    return None


def extract_g16_title(context: TextParseContext) -> str | None:
    focus_content, continued_content = g16_log_patterns.TITLE.split_content(context.content)
    if len(title_lines := focus_content.splitlines()) >= 3:
        context.content = continued_content
        return "\n".join(title_lines[1:-1]).replace("\n ", "")
    return None


def extract_g16_charge_multiplicity(
    context: TextParseContext,
) -> tuple[int, int] | tuple[None, None]:
    if matched := g16_log_patterns.CHARGE_MULTIPLICITY.search(context.content):
        charge = int(matched.group("charge"))
        multiplicity = int(matched.group("multiplicity"))
        return charge, multiplicity
    return None, None


def extract_g16_coordinates(
    context: TextParseContext, pattern: MolOPPattern
) -> np.ndarray[Any, Any] | None:
    focus_content = context.split(pattern)
    if matches := pattern.find_matches(focus_content):
        return np.array(
            [
                [
                    float(matched.group("x")),
                    float(matched.group("y")),
                    float(matched.group("z")),
                ]
                for matched in matches
            ],
        )
    return None


def extract_g16_standard_orientation_transformation_matrix(
    context: TextParseContext,
) -> np.ndarray[Any, Any] | None:
    coords = extract_g16_coordinates(context, g16_log_patterns.INITIAL_INPUT_COORDS)
    if coords is None:
        return None
    standard_coords = extract_g16_coordinates(context, g16_log_patterns.STANDARD_COORDS)
    if standard_coords is None:
        return None
    return find_rigid_transform(coords, standard_coords)


def extract_g16_solvent(context: TextParseContext) -> ImplicitSolvation | None:
    solvent_dict: dict[str, Any] = {}
    focus_content = context.split(g16_log_patterns.SOLVENT_PARAMETERS)
    if matches := g16_log_patterns.SOLVENT_MODEL.find_matches(focus_content):
        solvent_dict["solvent_model"] = matches[0].group("model")
    if matches := g16_log_patterns.SOLVENT_ATOM_RADII.find_matches(focus_content):
        solvent_dict["atomic_radii"] = matches[0].group("radii")
    if matches := g16_log_patterns.SOLVENT_TYPE.find_matches(focus_content):
        solvent_dict["solvent"] = matches[0].group("solvent")
    if matches := g16_log_patterns.SOLVENT_EPS.find_matches(focus_content):
        solvent_dict["solvent_epsilon"] = float(matches[0].group("value"))
    if matches := g16_log_patterns.SOLVENT_EPS_INF.find_matches(focus_content):
        solvent_dict["solvent_epsilon_infinite"] = float(matches[0].group("value"))
    if solvent_dict:
        return ImplicitSolvation.model_validate(solvent_dict)
    return None


def extract_g16_temperature_and_pressure(
    content: str,
) -> tuple[PlainQuantity, PlainQuantity] | tuple[None, None]:
    index = content.find(" - Thermochemistry -")
    if index == -1:
        return None, None
    focus_content = content[index:]
    if matches := g16_log_patterns.TEMPEREATURE_PRESSURE.find_content_matches(focus_content):
        matched = matches[0]
        temperature, pressure = (
            float(matched.group("temperature")) * atom_ureg.K,
            float(matched.group("pressure")) * atom_ureg.atm,
        )
        return temperature, pressure
    return None, None


def extract_g16_running_time(context: TextParseContext) -> PlainQuantity | None:
    if matches := g16_log_patterns.JOB_TIME.find_matches(context.content):
        total_seconds = 0.0
        for matched in matches:
            days = float(matched.group("days"))
            hours = float(matched.group("hours"))
            minutes = float(matched.group("minutes"))
            seconds = float(matched.group("seconds"))
            total_seconds += (days * 24 + hours) * 3600 + minutes * 60 + seconds
        return total_seconds * atom_ureg.second
    return None


def extract_g16_termination_status(context: TextParseContext) -> Status | None:
    status_dict: dict[str, Any] = {}
    if matches := g16_log_patterns.TERMINATION_STATUS.find_matches(context.content):
        status = matches[-1].group("status")
        if status == "Normal":
            status_dict["normal_terminated"] = True
            status_dict["scf_converged"] = True
        else:
            status_dict["normal_terminated"] = False
    if status_dict:
        return Status.model_validate(status_dict)
    return None
