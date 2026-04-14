from __future__ import annotations

import re
from collections.abc import Iterable, Iterator, Mapping
from dataclasses import dataclass, field
from functools import cache
from typing import TYPE_CHECKING, Any, ClassVar, Literal, Protocol

import numpy as np
from pint.facets.numpy.quantity import NumpyQuantity
from pint.facets.plain import PlainQuantity, PlainUnit
from pydantic import Field
from rdkit import Chem

from molop.config import moloplogger
from molop.io.base_models.Bases import BaseDataClassWithUnit
from molop.io.base_models.DataClasses import (
    ChargeSpinPopulations,
    Energies,
    GeometryOptimizationStatus,
    MolecularOrbitals,
    Polarizability,
    Status,
    ThermalInformations,
    TotalSpin,
    Vibrations,
)
from molop.io.base_models.SearchPattern import MolOPPattern, MolOPPatternV2
from molop.io.logic.QM_frame_parsers._g16_v2_shared import (
    ARCHIVE_TAIL_V2,
    FORCES_IN_CARTESIAN_V2,
    HESSIAN_IN_CARTESIAN_V2,
    INPUT_COORDS_V2,
    POPULATION_ANALYSIS_V2,
    STANDARD_COORDS_V2,
    _summarize_parse_context,
    extract_coords,
)
from molop.io.logic.QM_parsers._g16log_archive_tail import (
    parse_archive_tail,
)
from molop.io.patterns.G16Patterns import g16_log_patterns
from molop.unit import atom_ureg
from molop.utils.functions import fill_symmetric_matrix


RAW_MODEL_KEY = "__molop_model__"
RAW_DATA_KEY = "data"
RAW_QUANTITY_KEY = "__molop_quantity__"

pt = Chem.GetPeriodicTable()

_RAW_MODEL_REGISTRY: dict[str, type[BaseDataClassWithUnit]] = {}


def _register_raw_model(model_cls: type[BaseDataClassWithUnit]) -> type[BaseDataClassWithUnit]:
    _RAW_MODEL_REGISTRY[model_cls.__name__] = model_cls
    return model_cls


def _format_float(value: Any, digits: int = 6) -> str:
    if hasattr(value, "m"):
        value = value.m
    if isinstance(value, np.ndarray):
        value = float(value.reshape(-1)[0]) if value.size == 1 else float(value.flat[0])
    return f"{float(value):.{digits}f}"


def _format_job_time_seconds(seconds: float) -> str:
    days = int(seconds // 86400)
    seconds -= days * 86400
    hours = int(seconds // 3600)
    seconds -= hours * 3600
    minutes = int(seconds // 60)
    seconds -= minutes * 60
    return f"Job cpu time: {days} days {hours} hours {minutes} minutes {_format_float(seconds, 2)} seconds."


def _render_orientation_block(title: str, payload: Mapping[str, Any]) -> str:
    coords = payload.get("coords")
    if coords is None:
        coords = payload.get("standard_coords")
    atoms = payload.get("atoms", [])
    if coords is None or not atoms:
        return ""
    coords_array = coords.m if hasattr(coords, "m") else np.asarray(coords)
    lines = [
        f"{title}:",
        " ---------------------------------------------------------------------",
        " Center     Atomic      Atomic             Coordinates (Angstroms)",
        " Number     Number       Type             X           Y           Z",
        " ---------------------------------------------------------------------",
    ]
    for idx, (atom, xyz) in enumerate(zip(atoms, coords_array, strict=False), start=1):
        lines.append(
            f" {idx:5d} {int(atom):11d} {0:11d}"
            f" {float(xyz[0]):13.6f} {float(xyz[1]):11.6f} {float(xyz[2]):11.6f}"
        )
    lines.append(" ---------------------------------------------------------------------")
    return "\n".join(lines)


def _render_frequency_payload(payload: Mapping[str, Any]) -> str:
    vibrations = payload.get("vibrations")
    if vibrations is None:
        return ""
    vib_map = _payload_mapping(_restore_raw_payload(vibrations))
    lines: list[str] = [
        " Harmonic frequencies (cm**-1), IR intensities (KM/Mole), Raman scattering",
        " activities (A**4/AMU), depolarization ratios for plane and unpolarized",
        " incident light, reduced masses (AMU), force constants (mDyne/A),",
        " and normal coordinates:",
    ]
    frequencies = vib_map.get("frequencies")
    if frequencies is None or not len(frequencies):
        return ""
    atoms = payload.get("atoms") or []
    frequency_values = np.asarray(
        frequencies.m if hasattr(frequencies, "m") else frequencies
    ).reshape(-1)
    reduced_masses = vib_map.get("reduced_masses")
    reduced_mass_values = (
        np.asarray(reduced_masses.m if hasattr(reduced_masses, "m") else reduced_masses).reshape(-1)
        if reduced_masses is not None and len(reduced_masses)
        else None
    )
    force_constants = vib_map.get("force_constants")
    force_constant_values = (
        np.asarray(force_constants.m if hasattr(force_constants, "m") else force_constants).reshape(
            -1
        )
        if force_constants is not None and len(force_constants)
        else None
    )
    ir_intensities = vib_map.get("IR_intensities")
    ir_intensity_values = (
        np.asarray(ir_intensities.m if hasattr(ir_intensities, "m") else ir_intensities).reshape(-1)
        if ir_intensities is not None and len(ir_intensities)
        else None
    )
    vibration_modes = vib_map.get("vibration_modes")
    mode_arrays: list[np.ndarray[Any, Any]] = []
    if vibration_modes is not None and len(vibration_modes):
        restored_modes = vibration_modes.m if hasattr(vibration_modes, "m") else vibration_modes
        mode_arrays = [
            np.asarray(mode.m if hasattr(mode, "m") else mode) for mode in restored_modes
        ]

    for start in range(0, len(frequency_values), 3):
        stop = min(start + 3, len(frequency_values))
        chunk = frequency_values[start:stop]
        lines.append("")
        lines.append("".join(f"{mode_number:23d}" for mode_number in range(start + 1, stop + 1)))
        lines.append("".join(f"{'A':>23}" for _ in range(len(chunk))))
        lines.append(" Frequencies --" + "".join(f"{float(value):12.4f}" for value in chunk))
        if reduced_mass_values is not None:
            lines.append(
                " Red. masses --"
                + "".join(f"{float(value):12.4f}" for value in reduced_mass_values[start:stop])
            )
        if force_constant_values is not None:
            lines.append(
                " Frc consts  --"
                + "".join(f"{float(value):12.4f}" for value in force_constant_values[start:stop])
            )
        if ir_intensity_values is not None:
            lines.append(
                " IR Inten    --"
                + "".join(f"{float(value):12.4f}" for value in ir_intensity_values[start:stop])
            )
        lines.append("  Atom  AN" + "".join("      X      Y      Z" for _ in range(len(chunk))))
        if atoms and mode_arrays:
            for atom_index, atomic_number in enumerate(atoms, start=1):
                row = f"{atom_index:6d}{int(atomic_number):4d}"
                for mode_idx in range(start, stop):
                    if mode_idx < len(mode_arrays):
                        mode = np.asarray(mode_arrays[mode_idx])
                        vector = (
                            mode[atom_index - 1]
                            if mode.ndim > 1 and atom_index - 1 < len(mode)
                            else [0.0, 0.0, 0.0]
                        )
                        row += "".join(
                            f"{float(value):7.2f}" for value in np.asarray(vector).reshape(-1)[:3]
                        )
                    else:
                        row += f"{0.0:7.2f}{0.0:7.2f}{0.0:7.2f}"
                lines.append(row)
    lines.append(" -------------------")
    return "\n".join(lines)


for _model_cls in (
    Energies,
    TotalSpin,
    MolecularOrbitals,
    ChargeSpinPopulations,
    Polarizability,
    Vibrations,
    ThermalInformations,
    GeometryOptimizationStatus,
    Status,
):
    _register_raw_model(_model_cls)


@dataclass(frozen=True, slots=True)
class G16BoundarySpec:
    component_name: str
    parent_block: str = "g16.frame"
    start_patterns: tuple[str, ...] = ()
    literal_start_markers: tuple[str, ...] = ()
    contains_markers: tuple[str, ...] = ()
    end_patterns: tuple[str, ...] = ()
    end_offset: int = 0
    repeatable: bool = True
    kind: Literal["major", "child"] = "major"
    overlap_policy: Literal["forbid", "allow_shared_end", "parent_child"] = "forbid"
    skip_records: int = 0
    boundary_source: Literal["iochem", "molop", "temporary"] = "iochem"
    match_mode: Literal["line", "section", "molop_pattern"] = "line"
    include_end: bool = False
    molop_patterns: tuple[MolOPPatternV2, ...] = ()


@dataclass(slots=True)
class G16ScanContext:
    full_text: str
    _marker_presence_cache: dict[tuple[str, ...], bool] = field(default_factory=dict, repr=False)
    _marker_line_starts_cache: dict[tuple[str, ...], tuple[int, ...]] = field(
        default_factory=dict, repr=False
    )
    _single_marker_presence_cache: dict[str, bool] = field(default_factory=dict, repr=False)
    _single_marker_line_starts_cache: dict[str, tuple[int, ...]] = field(
        default_factory=dict, repr=False
    )

    def has_marker(self, marker: str) -> bool:
        if marker not in self._single_marker_presence_cache:
            self._single_marker_presence_cache[marker] = marker in self.full_text
        return self._single_marker_presence_cache[marker]

    def has_any_markers(self, markers: tuple[str, ...]) -> bool:
        if markers not in self._marker_presence_cache:
            self._marker_presence_cache[markers] = any(
                self.has_marker(marker) for marker in markers
            )
        return self._marker_presence_cache[markers]

    def marker_line_starts_for_marker(self, marker: str) -> tuple[int, ...]:
        if marker not in self._single_marker_line_starts_cache:
            self._single_marker_line_starts_cache[marker] = tuple(
                _iter_single_marker_line_starts(self.full_text, marker)
            )
        return self._single_marker_line_starts_cache[marker]

    def marker_line_starts(self, markers: tuple[str, ...]) -> tuple[int, ...]:
        if markers not in self._marker_line_starts_cache:
            seen: set[int] = set()
            line_starts: list[int] = []
            for marker in markers:
                for line_start in self.marker_line_starts_for_marker(marker):
                    if line_start not in seen:
                        seen.add(line_start)
                        line_starts.append(line_start)
            self._marker_line_starts_cache[markers] = tuple(line_starts)
        return self._marker_line_starts_cache[markers]


@dataclass(frozen=True, slots=True)
class G16SyntheticChildSpec:
    component_cls: type[G16V3BaseComponent]
    payload: Mapping[str, Any]
    node_name: str | None = None
    include_in_aggregation: bool = False
    include_in_render: bool = False


def _iter_pattern_spans(pattern: MolOPPatternV2, full_text: str) -> Iterator[tuple[int, int]]:
    cursor = 0
    while located := pattern.locate_content_from(full_text, cursor):
        start_start, _start_end, _end_start, end_end = located
        yield start_start, end_end
        cursor = max(end_end, start_start + 1)


def _line_end(full_text: str, start_pos: int) -> int:
    line_end = full_text.find("\n", start_pos)
    if line_end == -1:
        return len(full_text)
    return line_end + 1


def _line_start(full_text: str, position: int) -> int:
    line_start = full_text.rfind("\n", 0, position)
    return 0 if line_start == -1 else line_start + 1


def _advance_lines(full_text: str, position: int, count: int) -> int:
    end = position
    for _ in range(max(count, 0)):
        if end >= len(full_text):
            return len(full_text)
        end = _line_end(full_text, end)
    return end


@cache
def _compile_multiline_regex(pattern: str) -> re.Pattern[str]:
    return re.compile(pattern, re.MULTILINE)


def _normalize_line_match_start(full_text: str, position: int) -> int:
    normalized = position
    while normalized < len(full_text) and full_text[normalized] in "\r\n":
        normalized += 1
    return normalized


def _iter_single_marker_line_starts(full_text: str, marker: str) -> Iterator[int]:
    cursor = 0
    seen: set[int] = set()
    while True:
        position = full_text.find(marker, cursor)
        if position == -1:
            return
        line_start = _normalize_line_match_start(full_text, _line_start(full_text, position))
        if line_start not in seen:
            seen.add(line_start)
            yield line_start
        cursor = position + max(1, len(marker))


def _iter_marker_line_starts(full_text: str, markers: tuple[str, ...]) -> Iterator[int]:
    seen: set[int] = set()
    for marker in markers:
        for line_start in _iter_single_marker_line_starts(full_text, marker):
            if line_start not in seen:
                seen.add(line_start)
                yield line_start


def _match_literal_line_start(
    full_text: str, line_start: int, literal_markers: tuple[str, ...]
) -> tuple[int, int] | None:
    stripped_start = line_start
    while stripped_start < len(full_text) and full_text[stripped_start] in " \t":
        stripped_start += 1
    for marker in literal_markers:
        if full_text.startswith(marker, stripped_start):
            return line_start, stripped_start + len(marker)
    return None


def _iter_line_pattern_spans(
    pattern: str,
    full_text: str,
    candidate_starts: tuple[int, ...] | None = None,
    literal_start_markers: tuple[str, ...] = (),
) -> Iterator[tuple[int, int]]:
    if candidate_starts is not None:
        for line_start in candidate_starts:
            if literal_start_markers:
                matched = _match_literal_line_start(full_text, line_start, literal_start_markers)
                if matched is not None:
                    yield matched[0], _line_end(full_text, matched[0])
                    continue
            regex = _compile_multiline_regex(pattern)
            if regex.match(full_text, pos=line_start):
                yield line_start, _line_end(full_text, line_start)
        return
    regex = _compile_multiline_regex(pattern)
    for match in regex.finditer(full_text):
        line_start = _normalize_line_match_start(full_text, match.start())
        yield line_start, _line_end(full_text, line_start)


def _iter_section_spans_from_spec(
    spec: G16BoundarySpec, full_text: str, context: G16ScanContext | None = None
) -> Iterator[tuple[int, int]]:
    start_patterns = spec.start_patterns
    end_patterns = spec.end_patterns
    if not start_patterns:
        return

    start_regexes = tuple(_compile_multiline_regex(pattern) for pattern in start_patterns)
    end_regexes = tuple(_compile_multiline_regex(pattern) for pattern in end_patterns)

    starts: list[tuple[int, int]] = []
    candidate_starts = (
        context.marker_line_starts(spec.contains_markers)
        if context is not None and spec.contains_markers
        else tuple(_iter_marker_line_starts(full_text, spec.contains_markers))
    )
    if candidate_starts:
        for line_start in candidate_starts:
            if spec.literal_start_markers:
                matched = _match_literal_line_start(
                    full_text, line_start, spec.literal_start_markers
                )
                if matched is not None:
                    starts.append(matched)
                    continue
            for regex in start_regexes:
                match = regex.match(full_text, pos=line_start)
                if match is not None:
                    starts.append((line_start, match.end()))
                    break
    if not starts:
        for regex in start_regexes:
            starts.extend(
                (_normalize_line_match_start(full_text, match.start()), match.end())
                for match in regex.finditer(full_text)
            )
        starts.sort()

    for start, start_match_end in starts:
        if not end_patterns:
            yield start, len(full_text)
            continue

        matched_end: re.Match[str] | None = None
        for end_regex in end_regexes:
            current_match = end_regex.search(full_text, pos=start_match_end)
            if current_match is None:
                continue
            if matched_end is None or current_match.start() < matched_end.start():
                matched_end = current_match

        if matched_end is None:
            yield start, len(full_text)
            continue

        matched_end_start = _normalize_line_match_start(full_text, matched_end.start())
        if spec.include_end:
            end = _line_end(full_text, matched_end_start)
            end = _advance_lines(full_text, end, spec.end_offset)
        else:
            end = matched_end_start
        if end > start:
            yield start, end


def _iter_spans_from_spec(
    spec: G16BoundarySpec, full_text: str, context: G16ScanContext | None = None
) -> Iterator[tuple[int, int]]:
    has_markers = (
        context.has_any_markers(spec.contains_markers)
        if context is not None and spec.contains_markers
        else any(marker in full_text for marker in spec.contains_markers)
    )
    if spec.contains_markers and not has_markers:
        return
    if spec.match_mode == "molop_pattern":
        for pattern in spec.molop_patterns:
            yield from _iter_pattern_spans(pattern, full_text)
        return
    if spec.match_mode == "section":
        yield from _iter_section_spans_from_spec(spec, full_text, context=context)
        return
    candidate_starts = (
        (
            context.marker_line_starts(spec.contains_markers)
            if context is not None
            else tuple(_iter_marker_line_starts(full_text, spec.contains_markers))
        )
        if spec.contains_markers
        else None
    )
    for pattern in spec.start_patterns:
        yielded = False
        for span in _iter_line_pattern_spans(
            pattern,
            full_text,
            candidate_starts=candidate_starts,
            literal_start_markers=spec.literal_start_markers,
        ):
            yielded = True
            yield span
        if candidate_starts and not yielded:
            for span in _iter_line_pattern_spans(pattern, full_text, candidate_starts=None):
                yield span


def _temperature_and_pressure_from_block(block: str) -> dict[str, Any]:
    if matches := g16_log_patterns.TEMPEREATURE_PRESSURE.match_content(block):
        return {
            "temperature": float(matches[0][0]) * atom_ureg.K,
            "pressure": float(matches[0][1]) * atom_ureg.atm,
        }
    return {}


def _extract_input_coords_impl(
    block: str,
) -> tuple[list[int] | None, NumpyQuantity | None, str]:
    focus_content, continued_content = g16_log_patterns.INPUT_COORDS.split_content(block)
    if focus_content == "":
        return None, None, block
    if coords_match := g16_log_patterns.INPUT_COORDS.get_matches(focus_content):
        atoms, coords = extract_coords(coords_match)
        return atoms, coords, continued_content
    return None, None, block


def _extract_standard_coords_impl(
    block: str,
) -> tuple[list[int] | None, NumpyQuantity | None, str]:
    focus_content, continued_content = g16_log_patterns.STANDARD_COORDS.split_content(block)
    if focus_content == "":
        return None, None, block
    if coords_match := g16_log_patterns.STANDARD_COORDS.get_matches(focus_content):
        atoms, coords = extract_coords(coords_match)
        return atoms, coords, continued_content
    return None, None, continued_content


def _extract_rotation_consts_impl(block: str) -> NumpyQuantity | None:
    if matches := g16_log_patterns.ROTATIONAL_CONST.get_matches(block):
        return np.array(list(map(float, matches[0]))) * atom_ureg.gigahertz
    return None


def _parse_grouped_float_matches(matches: list[tuple[str, ...]]) -> np.ndarray:
    joined = " ".join(match[0] for match in matches)
    return np.fromstring(joined, sep=" ", dtype=float)


def _parse_orbital_line_values(energies: str) -> list[float]:
    return [
        float(value.replace("D", "E").replace("d", "E"))
        for value in re.findall(r"[-+]?\d+\.\d+(?:[DEde][-+]?\d+)?", energies)
    ]


def _parse_frequency_line_values(line: str) -> list[float]:
    if "--" in line:
        line = line.split("--", 1)[1]
    return [float(value.replace("D", "E").replace("d", "E")) for value in line.split()]


def _parse_version_payload_impl(block: str) -> str | None:
    focus_content, _continued_content = g16_log_patterns.VERSION.split_content(block)
    if matches := g16_log_patterns.VERSION.get_matches(focus_content):
        return matches[0][0]
    return None


def extract_input_coords_from_state(
    state: ParseState,
) -> tuple[list[int] | None, NumpyQuantity | None]:
    return _v2_extractors_module().extract_input_coords_from_state(state)


def extract_standard_coords_from_state(
    state: ParseState,
) -> tuple[list[int] | None, NumpyQuantity | None]:
    return _v2_extractors_module().extract_standard_coords_from_state(state)


def extract_energies_and_total_spin_from_state(
    state: ParseState,
) -> tuple[dict[str, Any] | None, dict[str, Any] | None]:
    return _v2_extractors_module().extract_energies_and_total_spin_from_state(state)


def extract_polarizability_from_state(state: ParseState) -> dict[str, Any] | None:
    return _v2_extractors_module().extract_polarizability_from_state(state)


def extract_populations_from_state(state: ParseState) -> dict[str, Any]:
    from molop.io.logic.QM_frame_parsers._g16_v2_extractors import (
        extract_populations_from_state as _impl,
    )

    return _impl(state)


def extract_vibrations_from_state(state: ParseState) -> dict[str, Any] | None:
    from molop.io.logic.QM_frame_parsers._g16_v2_extractors import (
        extract_vibrations_from_state as _impl,
    )

    return _impl(state)


def extract_thermal_infos_from_state(state: ParseState) -> dict[str, Any] | None:
    from molop.io.logic.QM_frame_parsers._g16_v2_extractors import (
        extract_thermal_infos_from_state as _impl,
    )

    return _impl(state)


def extract_forces_from_state(state: ParseState) -> NumpyQuantity | None:
    from molop.io.logic.QM_frame_parsers._g16_v2_extractors import (
        extract_forces_from_state as _impl,
    )

    return _impl(state)


def extract_hessian_from_state(state: ParseState) -> NumpyQuantity | None:
    from molop.io.logic.QM_frame_parsers._g16_v2_extractors import (
        extract_hessian_from_state as _impl,
    )

    return _impl(state)


def extract_berny_from_state(state: ParseState) -> GeometryOptimizationStatus | None:
    from molop.io.logic.QM_frame_parsers._g16_v2_extractors import extract_berny_from_state as _impl

    return _impl(state)


def extract_electric_dipole_and_polarizability_from_state(
    state: ParseState,
) -> dict[str, Any] | None:
    return _v2_extractors_module().extract_electric_dipole_and_polarizability_from_state(state)


def extract_tail_metadata_from_state(state: ParseState) -> dict[str, Any]:
    from molop.io.logic.QM_frame_parsers._g16_v2_extractors import (
        extract_tail_metadata_from_state as _impl,
    )

    return _impl(state)


def extract_tail_energies_from_state(state: ParseState) -> dict[str, Any] | None:
    from molop.io.logic.QM_frame_parsers._g16_v2_extractors import (
        extract_tail_energies_from_state as _impl,
    )

    return _impl(state)


def extract_tail_thermal_infos_from_state(state: ParseState) -> dict[str, Any] | None:
    from molop.io.logic.QM_frame_parsers._g16_v2_extractors import (
        extract_tail_thermal_infos_from_state as _impl,
    )

    return _impl(state)


def extract_tail_polarizability_from_state(state: ParseState) -> dict[str, Any] | None:
    from molop.io.logic.QM_frame_parsers._g16_v2_extractors import (
        extract_tail_polarizability_from_state as _impl,
    )

    return _impl(state)


def extract_tail_hessian_from_state(state: ParseState) -> NumpyQuantity | None:
    from molop.io.logic.QM_frame_parsers._g16_v2_extractors import (
        extract_tail_hessian_from_state as _impl,
    )

    return _impl(state)


def extract_rotation_consts_from_state(state: ParseState) -> NumpyQuantity | None:
    from molop.io.logic.QM_frame_parsers._g16_v2_extractors import (
        extract_rotation_consts_from_state as _impl,
    )

    return _impl(state)


def _parse_running_time(block: str) -> PlainQuantity | None:
    return G16V3JobCPUComponent._extract_running_time(block)


def _parse_running_time_impl(block: str) -> PlainQuantity | None:
    if matches := g16_log_patterns.PROCEDURE_TIME.match_content(block):
        total_seconds = 0.0
        for match in matches:
            total_seconds += float(match[1]) + float(match[2])
        return total_seconds * atom_ureg.second
    if matches := g16_log_patterns.JOB_TIME.match_content(block):
        total_seconds = 0.0
        for match in matches:
            total_seconds += (
                float(match[1]) * 24 * 3600
                + float(match[2]) * 3600
                + float(match[3]) * 60
                + float(match[4])
            )
        return total_seconds * atom_ureg.second
    return None


def _parse_jobcpu_payload(block: str) -> dict[str, Any]:
    payload: dict[str, Any] = {}
    if matches := g16_log_patterns.JOB_TIME.match_content(block):
        for match in matches:
            total_seconds = (
                float(match[1]) * 24 * 3600
                + float(match[2]) * 3600
                + float(match[3]) * 60
                + float(match[4])
            )
            label = match[0].strip().lower()
            if label == "job cpu time":
                payload["job_cpu_time"] = total_seconds * atom_ureg.second
            elif label == "elapsed time":
                payload["elapsed_time"] = total_seconds * atom_ureg.second
    if payload.get("job_cpu_time") is not None:
        payload["running_time"] = payload["job_cpu_time"]
    elif payload.get("elapsed_time") is not None:
        payload["running_time"] = payload["elapsed_time"]
    return payload


def _parse_termination_status_payload(block: str) -> dict[str, Any]:
    return G16V3L9999FinalComponent._extract_termination_status(block)


def _parse_termination_status_payload_impl(block: str) -> dict[str, Any]:
    if matches := g16_log_patterns.TERMINATION_STATUS.match_content(block):
        termination_kind = matches[-1][0]
        termination_line = block.strip().splitlines()[-1] if block.strip() else ""
        payload: dict[str, Any] = {
            "termination_kind": termination_kind,
            "termination_line": termination_line,
        }
        if termination_kind == "Normal":
            payload["status"] = Status(normal_terminated=True, scf_converged=True)
        else:
            payload["status"] = Status(normal_terminated=False, scf_converged=False)
        return payload
    return {}


def _payload_mapping(value: Any) -> Mapping[str, Any]:
    if isinstance(value, Mapping):
        return value
    if hasattr(value, "model_dump"):
        return value.model_dump()
    return {}


def _contains_any_marker(block: str, markers: tuple[str, ...]) -> bool:
    return any(marker in block for marker in markers)


if TYPE_CHECKING:
    from molop.io.logic.QM_frame_parsers._g16_v2_extractors import ParseState


def _v2_extractors_module():
    from molop.io.logic.QM_frame_parsers import _g16_v2_extractors

    return _g16_v2_extractors


def _rawify_payload(value: Any) -> Any:
    if isinstance(value, BaseDataClassWithUnit):
        model_fields = type(value).model_fields
        return {
            RAW_MODEL_KEY: value.__class__.__name__,
            RAW_DATA_KEY: {
                field_name: _rawify_payload(getattr(value, field_name))
                for field_name in model_fields
                if getattr(value, field_name) is not None
            },
        }
    if hasattr(value, "m") and hasattr(value, "units"):
        magnitude = value.m.tolist() if isinstance(value.m, np.ndarray) else value.m
        return {
            RAW_QUANTITY_KEY: True,
            "magnitude": _rawify_payload(magnitude),
            "unit": str(value.units),
        }
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, list):
        return [_rawify_payload(item) for item in value]
    if isinstance(value, tuple):
        return [_rawify_payload(item) for item in value]
    if isinstance(value, Mapping):
        return {key: _rawify_payload(item) for key, item in value.items()}
    return value


def _restore_raw_payload(value: Any) -> Any:
    if isinstance(value, list):
        return [_restore_raw_payload(item) for item in value]
    if isinstance(value, Mapping):
        if value.get(RAW_QUANTITY_KEY):
            magnitude = _restore_raw_payload(value["magnitude"])
            unit = atom_ureg.Unit(value["unit"])
            if isinstance(magnitude, list):
                return np.array(magnitude) * unit
            return magnitude * unit
        if RAW_MODEL_KEY in value:
            model_cls = _RAW_MODEL_REGISTRY[value[RAW_MODEL_KEY]]
            restored_data = {
                key: _restore_raw_payload(item) for key, item in value[RAW_DATA_KEY].items()
            }
            return model_cls.model_validate(restored_data)
        return {key: _restore_raw_payload(item) for key, item in value.items()}
    return value


def _merge_if_present(target: dict[str, Any], key: str, value: Any) -> None:
    if value is None:
        return
    if key in target and target[key] is not None:
        from molop.utils.functions import merge_models

        target[key] = merge_models(target[key], value, force_update=True)
    else:
        target[key] = value


def _merge_if_present_no_force(target: dict[str, Any], key: str, value: Any) -> None:
    if value is None:
        return
    if key in target and target[key] is not None:
        from molop.utils.functions import merge_models

        target[key] = merge_models(target[key], value)
    else:
        target[key] = value


def _has_meaningful_value(value: Any) -> bool:
    if value is None:
        return False
    if isinstance(value, str):
        return bool(value.strip())
    try:
        return len(value) > 0  # type: ignore[arg-type]
    except Exception:
        return True


class G16V3ComponentTreeProtocol(Protocol):
    only_extract_structure: bool


@dataclass(slots=True)
class G16V3BaseComponent:
    component_name: ClassVar[str] = "g16.v3.component"
    gaussian_block: ClassVar[str] = ""
    parent_block: ClassVar[str] = "g16.frame"
    repeatable: ClassVar[bool] = True
    boundaries: ClassVar[tuple[G16BoundarySpec, ...]] = ()
    allowed_child_component_names: ClassVar[tuple[str, ...]] = ()
    required_child_component_names: ClassVar[tuple[str, ...]] = ()
    repeatable_child_component_names: ClassVar[tuple[str, ...]] = ()
    required_frame_fields: ClassVar[tuple[str, ...]] = ()
    optional_frame_fields: ClassVar[tuple[str, ...]] = ()

    span_start: int = 0
    span_end: int = 0
    raw_text: str = ""
    payload: dict[str, Any] = field(default_factory=dict)

    @classmethod
    def boundary_specs(cls) -> tuple[G16BoundarySpec, ...]:
        return cls.boundaries

    @classmethod
    def primary_boundary_spec(cls) -> G16BoundarySpec:
        if cls.boundaries:
            return cls.boundaries[0]
        return G16BoundarySpec(
            component_name=cls.component_name,
            parent_block=cls.parent_block,
            repeatable=cls.repeatable,
            boundary_source="temporary",
        )

    @classmethod
    def boundary_markers(cls) -> tuple[str, ...]:
        markers: list[str] = []
        for spec in cls.boundary_specs():
            markers.extend(spec.contains_markers)
        return tuple(dict.fromkeys(markers))

    @classmethod
    def iter_spans(
        cls, full_text: str, *, context: G16ScanContext | None = None
    ) -> Iterable[tuple[int, int]]:
        spans: list[tuple[int, int]] = []
        seen: set[tuple[int, int]] = set()
        for spec in cls.boundary_specs():
            for span in _iter_spans_from_spec(spec, full_text, context=context):
                if span not in seen:
                    seen.add(span)
                    spans.append(span)
        return spans

    def parse_payload(
        self, full_text: str, tree_context: G16V3ComponentTreeProtocol
    ) -> dict[str, Any]:
        self.payload = {}
        return self.payload

    def build_synthetic_children(self) -> tuple[G16SyntheticChildSpec, ...]:
        return ()

    @classmethod
    def build_payload_from_frame(cls, frame: Any) -> dict[str, Any]:
        return {}

    @classmethod
    def can_build_from_frame(cls, frame: Any) -> bool:
        return all(
            getattr(frame, field_name, None) is not None for field_name in cls.required_frame_fields
        )

    def aggregate_into(self, infos: dict[str, Any], payload: Mapping[str, Any]) -> None:
        return None

    def _render_fakeg(self, **kwargs) -> str:
        return self.raw_text

    def render_fakeg(self, **kwargs) -> str:
        return self._render_fakeg(**kwargs)


class G16V3L1HeaderComponent(G16V3BaseComponent):
    component_name = "l1.header"
    gaussian_block = "l1.header"
    allowed_child_component_names = (
        "l1.options",
        "l1.keywords",
        "l101.title",
        "l101.charge_multiplicity",
    )
    repeatable = False
    optional_frame_fields = (
        "qm_software_version",
        "options",
        "keywords",
        "title_card",
        "charge",
        "multiplicity",
    )
    boundaries = (
        G16BoundarySpec(
            component_name="l1.header",
            start_patterns=(
                r"^\s*(?:Entering Link 1 = .*|Link1:\s+Proceeding to internal job step number\s+\d+\.)\s*$",
            ),
            literal_start_markers=(
                "Entering Link 1 =",
                "Link1:  Proceeding to internal job step number",
            ),
            contains_markers=(
                "Entering Link 1 =",
                "Link1:  Proceeding to internal job step number",
            ),
            end_patterns=(
                r"^\s*Input orientation:",
                r"^\s*Standard orientation:",
                r"^\s*1[\\|]1[\\|]GINC",
                r"^\s*(Job cpu time|Elapsed time):",
                r"^\s*Normal termination of Gaussian",
                r"^\s*Leave Link",
            ),
            repeatable=False,
            boundary_source="iochem",
            match_mode="section",
            include_end=False,
        ),
    )

    @staticmethod
    def _extract_options(block: str) -> tuple[str | None, str]:
        focus_content, continued_content = g16_log_patterns.OPTIONS.split_content(block)
        if matches := g16_log_patterns.OPTIONS.get_matches(focus_content):
            return "\n".join(f"{match[0]}={match[1]}" for match in matches), continued_content
        return None, block

    @staticmethod
    def _extract_keywords(block: str) -> tuple[str | None, str]:
        focus_content, continued_content = g16_log_patterns.KEYWORDS.split_content(block)
        keyword_lines = focus_content.splitlines()
        if len(keyword_lines) >= 3:
            return "\n".join(keyword_lines[1:-1]).replace("\n ", ""), continued_content
        return None, block

    @staticmethod
    def _extract_title(block: str) -> tuple[str | None, str]:
        focus_content, continued_content = g16_log_patterns.TITLE.split_content(block)
        title_lines = focus_content.splitlines()
        if len(title_lines) >= 3:
            return "\n".join(title_lines[1:-1]).replace("\n ", ""), continued_content
        return None, block

    @staticmethod
    def _extract_charge_multiplicity(block: str) -> tuple[tuple[int, int] | tuple[None, None], str]:
        if match := g16_log_patterns.CHARGE_MULTIPLICITY.match_content(block):
            return (int(match[0][0]), int(match[0][1])), block
        return (None, None), block

    @classmethod
    def _parse_header_payload(cls, block: str) -> dict[str, Any]:
        payload: dict[str, Any] = {}
        working_block = block
        if version := _parse_version_payload_impl(working_block):
            payload["qm_software_version"] = version
        options, working_block = cls._extract_options(working_block)
        if options is not None:
            payload["options"] = options
        keywords, working_block = cls._extract_keywords(working_block)
        if keywords is not None:
            payload["keywords"] = keywords
        title, working_block = cls._extract_title(working_block)
        if title is not None:
            payload["title_card"] = title
        charge_multiplicity, _working_block = cls._extract_charge_multiplicity(working_block)
        charge, multiplicity = charge_multiplicity
        if charge is not None:
            payload["charge"] = charge
        if multiplicity is not None:
            payload["multiplicity"] = multiplicity
        return payload

    def parse_payload(
        self, full_text: str, tree_context: G16V3ComponentTreeProtocol
    ) -> dict[str, Any]:
        payload: dict[str, Any] = self._parse_header_payload(self.raw_text)
        if running_time := _parse_running_time(self.raw_text):
            payload["running_time"] = running_time
        self.payload = payload
        return self.payload

    def build_synthetic_children(self) -> tuple[G16SyntheticChildSpec, ...]:
        specs: list[G16SyntheticChildSpec] = []
        if self.payload.get("options") is not None:
            specs.append(
                G16SyntheticChildSpec(
                    component_cls=G16V3L1OptionsComponent,
                    payload={"options": self.payload["options"]},
                )
            )
        if self.payload.get("keywords") is not None:
            specs.append(
                G16SyntheticChildSpec(
                    component_cls=G16V3L1KeywordsComponent,
                    payload={"keywords": self.payload["keywords"]},
                )
            )
        if self.payload.get("title_card") is not None:
            specs.append(
                G16SyntheticChildSpec(
                    component_cls=G16V3L101TitleComponent,
                    payload={"title_card": self.payload["title_card"]},
                )
            )
        charge_mult_payload = {
            key: self.payload[key]
            for key in ("charge", "multiplicity")
            if self.payload.get(key) is not None
        }
        if charge_mult_payload:
            specs.append(
                G16SyntheticChildSpec(
                    component_cls=G16V3L101ChargeMultiplicityComponent,
                    payload=charge_mult_payload,
                )
            )
        return tuple(specs)

    def aggregate_into(self, infos: dict[str, Any], payload: Mapping[str, Any]) -> None:
        for key in (
            "qm_software_version",
            "options",
            "keywords",
            "title_card",
            "charge",
            "multiplicity",
        ):
            if key not in infos and payload.get(key) is not None:
                infos[key] = payload[key]
        if payload.get("running_time") is not None and "running_time" not in infos:
            infos["running_time"] = payload["running_time"]

    @classmethod
    def build_payload_from_frame(cls, frame: Any) -> dict[str, Any]:
        payload: dict[str, Any] = {}
        for key in (
            "qm_software_version",
            "options",
            "keywords",
            "title_card",
            "charge",
            "multiplicity",
            "running_time",
        ):
            if (value := getattr(frame, key, None)) is not None:
                payload[key] = value
        return payload

    def _render_fakeg(self, **kwargs) -> str:
        running_time = self.payload.get("running_time")
        if running_time is None:
            return self.raw_text
        seconds = float(
            running_time.to("second").m if hasattr(running_time, "to") else running_time
        )
        return _format_job_time_seconds(seconds)


class G16V3L101TitleComponent(G16V3BaseComponent):
    component_name = "l101.title"
    gaussian_block = "l101.title"
    repeatable = False


class G16V3L1OptionsComponent(G16V3BaseComponent):
    component_name = "l1.options"
    gaussian_block = "l1.options"
    repeatable = False


class G16V3L1KeywordsComponent(G16V3BaseComponent):
    component_name = "l1.keywords"
    gaussian_block = "l1.keywords"
    repeatable = False


class G16V3L101ChargeMultiplicityComponent(G16V3BaseComponent):
    component_name = "l101.charge_multiplicity"
    gaussian_block = "l101.charge_multiplicity"
    repeatable = False

    def _render_fakeg(self, **kwargs) -> str:
        charge = self.payload.get("charge")
        multiplicity = self.payload.get("multiplicity")
        if charge is None or multiplicity is None:
            return self.raw_text
        return f" Charge = {int(charge):4d} Multiplicity = {int(multiplicity)}"


class G16V3L202OrientComponent(G16V3BaseComponent):
    component_name = "l202.orient"
    gaussian_block = "l202.orient"
    allowed_child_component_names = (
        "l202.orient.input",
        "l202.orient.standard",
        "l202.distmat",
        "l202.stoich",
    )
    required_frame_fields = ("atoms",)
    optional_frame_fields = ("coords", "standard_coords")
    boundaries = (
        G16BoundarySpec(
            component_name="l202.orient",
            start_patterns=(r"^\s*Input orientation:", r"^\s*Standard orientation:"),
            literal_start_markers=("Input orientation:", "Standard orientation:"),
            contains_markers=("Input orientation:", "Standard orientation:"),
            end_patterns=(r"^\s*-+\s*$",),
            end_offset=2,
            repeatable=True,
            boundary_source="iochem",
            match_mode="molop_pattern",
            molop_patterns=(INPUT_COORDS_V2, STANDARD_COORDS_V2),
        ),
    )

    @staticmethod
    def _extract_input_orientation(
        block: str,
    ) -> tuple[list[int] | None, NumpyQuantity | None, str]:
        return _extract_input_coords_impl(block)

    @staticmethod
    def _extract_standard_orientation(
        block: str,
    ) -> tuple[list[int] | None, NumpyQuantity | None, str]:
        return _extract_standard_coords_impl(block)

    def parse_payload(
        self, full_text: str, tree_context: G16V3ComponentTreeProtocol
    ) -> dict[str, Any]:
        payload: dict[str, Any] = {}
        atoms, coords, _remaining = self._extract_input_orientation(self.raw_text)
        if atoms and coords is not None:
            payload["atoms"] = atoms
            payload["coords"] = coords
        atoms, standard_coords, _remaining = self._extract_standard_orientation(self.raw_text)
        if atoms and standard_coords is not None:
            payload["atoms"] = atoms
            payload["standard_coords"] = standard_coords
        self.payload = payload
        return self.payload

    def _render_fakeg(self, **kwargs) -> str:
        if "standard_coords" in self.payload:
            return _render_orientation_block("Standard orientation", self.payload)
        if "coords" in self.payload:
            return _render_orientation_block("Input orientation", self.payload)
        return self.raw_text

    def aggregate_into(self, infos: dict[str, Any], payload: Mapping[str, Any]) -> None:
        if payload.get("atoms"):
            infos["atoms"] = payload["atoms"]
        if payload.get("coords") is not None:
            infos["coords"] = payload["coords"]
        if payload.get("standard_coords") is not None:
            infos["standard_coords"] = payload["standard_coords"]
            if payload.get("atoms"):
                infos["atoms"] = payload["atoms"]

    @classmethod
    def build_payload_from_frame(cls, frame: Any) -> dict[str, Any]:
        payload: dict[str, Any] = {}
        atoms = getattr(frame, "atoms", None)
        coords = getattr(frame, "coords", None)
        standard_coords = getattr(frame, "standard_coords", None)
        if atoms and coords is not None:
            payload["atoms"] = atoms
            payload["coords"] = coords
        if atoms and standard_coords is not None:
            payload["atoms"] = atoms
            payload["standard_coords"] = standard_coords
        return payload


class G16V3L202RotConstComponent(G16V3BaseComponent):
    component_name = "l202.rotconst"
    gaussian_block = "l202.rotconst"
    boundaries = (
        G16BoundarySpec(
            component_name="l202.rotconst",
            start_patterns=(r"^\s*Rotational constants \(GHZ\):",),
            literal_start_markers=("Rotational constants (GHZ):",),
            contains_markers=("Rotational constants (GHZ):",),
            repeatable=True,
            boundary_source="molop",
            match_mode="line",
        ),
    )
    required_frame_fields = ("rotation_constants",)

    @staticmethod
    def _extract_rotation_constants(block: str) -> NumpyQuantity | None:
        return _extract_rotation_consts_impl(block)

    def parse_payload(
        self, full_text: str, tree_context: G16V3ComponentTreeProtocol
    ) -> dict[str, Any]:
        payload: dict[str, Any] = {}
        rotation_constants = self._extract_rotation_constants(self.raw_text)
        if rotation_constants is not None:
            payload["rotation_constants"] = rotation_constants
        self.payload = payload
        return self.payload

    def aggregate_into(self, infos: dict[str, Any], payload: Mapping[str, Any]) -> None:
        if payload.get("rotation_constants") is not None:
            infos["rotation_constants"] = payload["rotation_constants"]

    @classmethod
    def build_payload_from_frame(cls, frame: Any) -> dict[str, Any]:
        if (rotation_constants := getattr(frame, "rotation_constants", None)) is not None:
            return {"rotation_constants": rotation_constants}
        return {}


class G16V3L502CycleComponent(G16V3BaseComponent):
    component_name = "l502.cycle"
    gaussian_block = "l502.cycle"
    boundaries = (
        G16BoundarySpec(
            component_name="l502.cycle",
            start_patterns=(r"^\s*Cycle\s+.*",),
            literal_start_markers=("Cycle",),
            contains_markers=("Cycle", "SCF Done:"),
            end_patterns=(r"^\s*SCF Done:.*",),
            repeatable=True,
            boundary_source="iochem",
            match_mode="section",
            include_end=True,
        ),
        G16BoundarySpec(
            component_name="l502.cycle",
            start_patterns=(r"^\s*SCF Done:",),
            literal_start_markers=("SCF Done:",),
            contains_markers=("SCF Done:",),
            repeatable=True,
            boundary_source="temporary",
            match_mode="line",
        ),
    )
    required_frame_fields = ("energies",)
    optional_frame_fields = ("total_spin",)

    @classmethod
    def iter_spans(
        cls, full_text: str, *, context: G16ScanContext | None = None
    ) -> Iterable[tuple[int, int]]:
        for spec in cls.boundary_specs():
            spans = list(dict.fromkeys(_iter_spans_from_spec(spec, full_text, context=context)))
            if spans:
                return spans
        return []

    @staticmethod
    def _extract_energies_and_total_spin(
        block: str,
    ) -> tuple[Energies | None, TotalSpin | None, str]:
        scf_energies_dict: dict[str, PlainQuantity | None] = {}
        total_spin_dict: dict[str, float | None] = {}
        focus_content, continued_content = g16_log_patterns.SCF_ENERGIES.split_content(block)
        if focus_content == "":
            return (
                Energies.model_validate(scf_energies_dict) if scf_energies_dict else None,
                TotalSpin.model_validate(total_spin_dict) if total_spin_dict else None,
                block,
            )
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
            Energies.model_validate(scf_energies_dict) if scf_energies_dict else None,
            TotalSpin.model_validate(total_spin_dict) if total_spin_dict else None,
            continued_content,
        )

    @classmethod
    def _parse_cycle_payload(cls, block: str) -> dict[str, Any]:
        payload: dict[str, Any] = {}
        energies, total_spin, _remaining = cls._extract_energies_and_total_spin(block)
        if energies is not None:
            payload["energies"] = energies
        if total_spin is not None:
            payload["total_spin"] = total_spin
        return payload

    def parse_payload(
        self, full_text: str, tree_context: G16V3ComponentTreeProtocol
    ) -> dict[str, Any]:
        if tree_context.only_extract_structure:
            self.payload = {}
            return self.payload
        self.payload = self._parse_cycle_payload(self.raw_text)
        return self.payload

    def _render_fakeg(self, **kwargs) -> str:
        lines: list[str] = []
        energies = self.payload.get("energies")
        if energies is not None and getattr(energies, "scf_energy", None) is not None:
            lines.append(
                f" SCF Done:  E(SCF) =  {_format_float(energies.scf_energy, 9)} A.U. after   1 cycles"
            )
        total_spin = self.payload.get("total_spin")
        if total_spin is not None and getattr(total_spin, "spin_square", None) is not None:
            line = f" S**2 = {_format_float(total_spin.spin_square, 6)}"
            if getattr(total_spin, "spin_quantum_number", None) is not None:
                line += f"                 S = {_format_float(total_spin.spin_quantum_number, 6)}"
            lines.append(line)
        return "\n".join(lines) if lines else self.raw_text

    def aggregate_into(self, infos: dict[str, Any], payload: Mapping[str, Any]) -> None:
        _merge_if_present(infos, "energies", payload.get("energies"))
        _merge_if_present(infos, "total_spin", payload.get("total_spin"))

    @classmethod
    def build_payload_from_frame(cls, frame: Any) -> dict[str, Any]:
        payload: dict[str, Any] = {}
        if (energies := getattr(frame, "energies", None)) is not None:
            payload["energies"] = energies
        if (total_spin := getattr(frame, "total_spin", None)) is not None:
            payload["total_spin"] = total_spin
        return payload


class G16V3L601PopAnalComponent(G16V3BaseComponent):
    component_name = "l601.popanal"
    gaussian_block = "l601.popanal"
    allowed_child_component_names = (
        "l601.state",
        "l601.molecular_orbitals",
        "l601.charge_spin_populations",
        "l601.polarizability",
    )
    optional_frame_fields = (
        "molecular_orbitals",
        "charge_spin_populations",
        "polarizability",
    )
    boundaries = (
        G16BoundarySpec(
            component_name="l601.popanal",
            start_patterns=(r"^\s*Population analysis using the (SCF|CC) [Dd]ensity\.",),
            literal_start_markers=("Population analysis using the",),
            contains_markers=("Population analysis using the",),
            end_patterns=(r"^\s*N-N=.*",),
            repeatable=True,
            boundary_source="iochem",
            match_mode="molop_pattern",
            molop_patterns=(POPULATION_ANALYSIS_V2,),
        ),
    )

    _MO_SYMMETRY_MARKERS = ("Orbital symmetries:",)
    _ELECTRONIC_STATE_MARKERS = ("The electronic state is",)
    _MO_ENERGY_MARKERS = (
        "Alpha  occ. eigenvalues --",
        "Beta  occ. eigenvalues --",
        "Alpha virt. eigenvalues --",
        "Beta virt. eigenvalues --",
    )
    _POPULATION_SECTION_MARKERS = {
        "mulliken_spins": ("Mulliken charges and spin densities:",),
        "mulliken_charges": (
            "Mulliken charges:",
            "Mulliken atomic charges",
            "Mulliken charges and spin densities:",
        ),
        "apt_charges": ("APT charges:",),
        "lowdin_charges": ("Lowdin charges",),
    }
    _POPULATION_SECTION_END_MARKERS = {
        "apt_charges": ("Sum of APT charges",),
        "lowdin_charges": ("Sum of Lowdin charges",),
    }
    _POLAR_SECTION_MARKERS = {
        "electronic_spatial_extent": ("Electronic spatial extent (au):",),
        "dipole": ("Dipole moment (field-independent basis",),
        "quadrupole": ("Quadrupole moment (field-independent basis",),
        "traceless_quadrupole": ("Traceless Quadrupole moment (field-independent basis",),
        "octapole": ("Octapole moment (field-independent basis",),
        "hexadecapole": ("Hexadecapole moment (field-independent basis",),
    }
    _REMAINING_BLOCK_MARKERS = {
        "exact_polarizability": ("Exact polarizability:",),
        "hirshfeld": ("Hirshfeld charges, spin densities, dipoles, and CM5 charges",),
        "dipole_before_force": ("Dipole        =",),
        "polarizability_before_force": ("Polarizability=",),
    }

    @staticmethod
    def _extract_population_payload(block: str) -> tuple[dict[str, Any], str]:
        infos: dict[str, Any] = {}
        mo: dict[str, Any] = {}
        pops: dict[str, Any] = {}
        polars: dict[str, Any] = {}
        remaining_block = block
        try:
            focus_content, remaining_block = g16_log_patterns.POPULATION_ANALYSIS.split_content(
                block
            )
            if focus_content == "":
                return infos, block
            patterns_and_keys_1: list[tuple[MolOPPattern, str]] = [
                (g16_log_patterns.MOLECULAR_ORBITALS_SYMMETRY_ALPHA, "alpha_symmetries"),
                (g16_log_patterns.MOLECULAR_ORBITALS_SYMMETRY_BETA, "beta_symmetries"),
                (g16_log_patterns.MOLECULAR_ORBITALS_SYMMETRY, "alpha_symmetries"),
            ]
            if _contains_any_marker(focus_content, G16V3L601PopAnalComponent._MO_SYMMETRY_MARKERS):
                for pattern, key in patterns_and_keys_1:
                    sub_focus_content, focus_content = pattern.split_content(focus_content)
                    if matches := pattern.get_matches(sub_focus_content):
                        mo[key] = [sym for match in matches for sym in match[1].split()]
            if _contains_any_marker(
                focus_content, G16V3L601PopAnalComponent._ELECTRONIC_STATE_MARKERS
            ) and (matches := g16_log_patterns.ELECTRONIC_STATE.get_matches(focus_content)):
                mo["electronic_state"] = matches[0][0]
            if _contains_any_marker(
                focus_content, G16V3L601PopAnalComponent._MO_ENERGY_MARKERS
            ) and (matches := g16_log_patterns.MOLECULAR_ORBITALS.get_matches(focus_content)):
                temp_alpha_orbitals = []
                temp_alpha_occupancies = []
                temp_beta_orbitals = []
                temp_beta_occupancies = []
                for orbital_type, occ_stat, energies in matches:
                    parsed_energies = _parse_orbital_line_values(energies)
                    if occ_stat not in ["occ.", "virt."]:
                        raise ValueError(f"Invalid {orbital_type} occupancy status: {occ_stat}")
                    if orbital_type == "Alpha":
                        temp_alpha_orbitals.extend(parsed_energies)
                        temp_alpha_occupancies.extend([occ_stat == "occ."] * len(parsed_energies))
                    elif orbital_type == "Beta":
                        temp_beta_orbitals.extend(parsed_energies)
                        temp_beta_occupancies.extend([occ_stat == "occ."] * len(parsed_energies))
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
                if not _contains_any_marker(
                    focus_content, G16V3L601PopAnalComponent._POPULATION_SECTION_MARKERS[key]
                ):
                    continue
                sub_focus_content, focus_content = pattern.split_content(focus_content)
                if matches := pattern.get_matches(sub_focus_content):
                    pops[key] = [
                        float(match[0] if key != "mulliken_spins" else match[1])
                        for match in matches
                    ]
            if _contains_any_marker(
                focus_content,
                G16V3L601PopAnalComponent._POLAR_SECTION_MARKERS["electronic_spatial_extent"],
            ):
                sub_focus_content, focus_content = (
                    g16_log_patterns.ELECTRONIC_SPATIAL_EXTENT.split_content(focus_content)
                )
            else:
                sub_focus_content = ""
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
                if not _contains_any_marker(
                    focus_content, G16V3L601PopAnalComponent._POLAR_SECTION_MARKERS[key]
                ):
                    continue
                sub_focus_content, focus_content = pattern.split_content(focus_content)
                if matches := pattern.get_matches(sub_focus_content):
                    polars[key] = np.array([float(match[0]) for match in matches]) * unit
            if _contains_any_marker(
                remaining_block,
                G16V3L601PopAnalComponent._REMAINING_BLOCK_MARKERS["exact_polarizability"],
            ) and (matches := g16_log_patterns.EXACT_POLARIZABILITY.get_matches(remaining_block)):
                polars["polarizability_tensor"] = (
                    np.array([float(match) for match in matches[0]]) * atom_ureg.bohr**3
                )
            if _contains_any_marker(
                remaining_block, G16V3L601PopAnalComponent._REMAINING_BLOCK_MARKERS["hirshfeld"]
            ):
                sub_focus_content, remaining_block = (
                    g16_log_patterns.HIRSHFELD_POPULATION.split_content(remaining_block)
                )
            else:
                sub_focus_content = ""
            if matches := g16_log_patterns.HIRSHFELD_POPULATION.get_matches(sub_focus_content):
                pops["hirshfeld_charges"] = [float(match[0]) for match in matches]
                pops["hirshfeld_spins"] = [float(match[1]) for match in matches]
                pops["hirshfeld_q_cm5"] = [float(match[5]) for match in matches]
            if _contains_any_marker(
                remaining_block,
                G16V3L601PopAnalComponent._REMAINING_BLOCK_MARKERS["dipole_before_force"],
            ) and (matches := g16_log_patterns.DIPOLE_BEFORE_FORCE.match_content(remaining_block)):
                polars["dipole"] = (
                    np.array([float(match[0].replace("D", "E")) for match in matches])
                    * atom_ureg.atomic_unit_of_current
                    * atom_ureg.atomic_unit_of_time
                    * atom_ureg.bohr
                )
            if _contains_any_marker(
                remaining_block,
                G16V3L601PopAnalComponent._REMAINING_BLOCK_MARKERS["polarizability_before_force"],
            ) and (
                matches := g16_log_patterns.POLARIZIABILITIES_BEFORE_FORCE.match_content(
                    remaining_block
                )
            ):
                polars["polarizability_tensor"] = (
                    np.array([float(match[0].replace("D", "E")) for match in matches])
                    * atom_ureg.bohr**3
                )
            if mo:
                infos["molecular_orbitals"] = MolecularOrbitals.model_validate(mo)
            if pops:
                infos["charge_spin_populations"] = ChargeSpinPopulations.model_validate(pops)
            if polars:
                infos["polarizability"] = Polarizability.model_validate(polars)
        except (ValueError, IndexError) as exc:
            moloplogger.warning(
                "Error parsing populations: %s | context=%s",
                exc,
                _summarize_parse_context(remaining_block),
            )
        except Exception as exc:
            moloplogger.error(f"Unexpected error occurred while parsing populations: {exc}")
        return infos, remaining_block

    @classmethod
    def _parse_population_payload(cls, block: str) -> dict[str, Any]:
        payload, _remaining = cls._extract_population_payload(block)
        return payload

    def parse_payload(
        self, full_text: str, tree_context: G16V3ComponentTreeProtocol
    ) -> dict[str, Any]:
        if tree_context.only_extract_structure:
            self.payload = {}
            return self.payload
        self.payload = self._parse_population_payload(self.raw_text)
        return self.payload

    def aggregate_into(self, infos: dict[str, Any], payload: Mapping[str, Any]) -> None:
        _merge_if_present(infos, "molecular_orbitals", payload.get("molecular_orbitals"))
        _merge_if_present(infos, "charge_spin_populations", payload.get("charge_spin_populations"))
        _merge_if_present(infos, "polarizability", payload.get("polarizability"))

    @classmethod
    def build_payload_from_frame(cls, frame: Any) -> dict[str, Any]:
        payload: dict[str, Any] = {}
        if (molecular_orbitals := getattr(frame, "molecular_orbitals", None)) is not None:
            payload["molecular_orbitals"] = molecular_orbitals
        if (charge_spin_populations := getattr(frame, "charge_spin_populations", None)) is not None:
            payload["charge_spin_populations"] = charge_spin_populations
        if (polarizability := getattr(frame, "polarizability", None)) is not None:
            payload["polarizability"] = polarizability
        return payload


class G16V3L716FreqComponent(G16V3BaseComponent):
    component_name = "l716.freq"
    gaussian_block = "l716.freq"
    allowed_child_component_names = (
        "l716.forceconstants",
        "l716.diagvib",
        "l716.irspectrum",
        "l716.vibration.mode",
    )

    repeatable_child_component_names = ("l716.vibration.mode",)
    required_frame_fields = ("vibrations",)
    boundaries = (
        G16BoundarySpec(
            component_name="l716.freq",
            start_patterns=(r"Harmonic frequencies \(cm\*\*-1\)",),
            literal_start_markers=("Harmonic frequencies (cm**-1)",),
            contains_markers=("Harmonic frequencies (cm**-1)", "Frequencies --"),
            end_patterns=(
                r"^\s*- Thermochemistry -",
                r"^\s*Center\s+Atomic\s+Forces \(Hartrees/Bohr\)",
                r"^\s*Center\s+Atomic\s+\s+\s+Forces \(Hartrees/Bohr\)",
                r"^\s*Force constants in Cartesian coordinates",
                r"^\s*Item\s+Value\s+Threshold\s+Converged\?",
                r"^\s*1[\\|]1[\\|]GINC",
                r"^\s*(Job cpu time|Elapsed time):",
                r"^\s*Normal termination of Gaussian",
                r"^\s*Leave Link",
            ),
            repeatable=True,
            boundary_source="iochem",
            match_mode="section",
            include_end=False,
        ),
    )

    @staticmethod
    def _extract_vibration_payload(block: str) -> tuple[Vibrations | None, str]:
        vib_dict: dict[str, Any] = {}
        start_index = block.find(
            "Harmonic frequencies (cm**-1), IR intensities (KM/Mole), Raman scattering"
        )
        end_index = block.find("-------------------", start_index)
        if start_index == -1 or end_index == -1:
            focus_content, remaining_block = g16_log_patterns.FREQUENCY_ANALYSIS.split_content(
                block
            )
        else:
            focus_content = block[start_index:end_index]
            remaining_block = block[end_index + 1 :]
        if focus_content == "":
            return None, block
        frequencies: list[float] = []
        reduced_masses: list[float] = []
        force_constants: list[float] = []
        ir_intensities: list[float] = []
        vibration_modes: list[list[list[float]]] = []
        lines = focus_content.splitlines()
        idx = 0
        total_modes = 0
        while idx < len(lines):
            line = lines[idx]
            stripped = line.strip()
            if not stripped.startswith("Frequencies --"):
                idx += 1
                continue

            chunk_freqs = _parse_frequency_line_values(stripped)
            chunk_size = len(chunk_freqs)
            if chunk_size == 0:
                idx += 1
                continue
            frequencies.extend(chunk_freqs)
            total_modes += chunk_size

            if idx + 1 < len(lines) and lines[idx + 1].strip().startswith("Red. masses --"):
                reduced_masses.extend(_parse_frequency_line_values(lines[idx + 1].strip()))
            if idx + 2 < len(lines) and lines[idx + 2].strip().startswith("Frc consts  --"):
                force_constants.extend(_parse_frequency_line_values(lines[idx + 2].strip()))
            if idx + 3 < len(lines) and lines[idx + 3].strip().startswith("IR Inten    --"):
                ir_intensities.extend(_parse_frequency_line_values(lines[idx + 3].strip()))

            mode_rows: list[list[list[float]]] = [[] for _ in range(chunk_size)]
            scan_idx = idx + 4
            if scan_idx < len(lines) and lines[scan_idx].strip().startswith("Atom  AN"):
                scan_idx += 1
            while scan_idx < len(lines):
                row = lines[scan_idx]
                row_stripped = row.strip()
                if not row_stripped:
                    break
                if row_stripped.startswith("Frequencies --") or row_stripped.startswith(
                    "- Thermochemistry"
                ):
                    break
                parts = row.split()
                if len(parts) >= 2 + 3 * chunk_size and parts[0].isdigit() and parts[1].isdigit():
                    values = [
                        float(value.replace("D", "E").replace("d", "E"))
                        for value in parts[2 : 2 + 3 * chunk_size]
                    ]
                    for mode_idx in range(chunk_size):
                        start = mode_idx * 3
                        mode_rows[mode_idx].append(values[start : start + 3])
                    scan_idx += 1
                    continue
                if parts and all(token.isdigit() for token in parts):
                    scan_idx += 1
                    continue
                if parts and all(re.fullmatch(r"[A-Za-z][\w-]*", token) for token in parts):
                    scan_idx += 1
                    continue
                break
            vibration_modes.extend(mode_rows)
            idx = scan_idx

        if frequencies:
            vib_dict["frequencies"] = np.array(frequencies) * atom_ureg.cm_1
        if reduced_masses:
            vib_dict["reduced_masses"] = np.array(reduced_masses) * atom_ureg.amu
        if force_constants:
            vib_dict["force_constants"] = (
                np.array(force_constants) * atom_ureg.mdyne / atom_ureg.angstrom
            )
        if ir_intensities:
            vib_dict["IR_intensities"] = np.array(ir_intensities) * atom_ureg.km / atom_ureg.mol
        if vibration_modes and total_modes:
            vib_dict["vibration_modes"] = [
                np.array(mode_rows) for mode_rows in vibration_modes
            ] * atom_ureg.angstrom
        return (Vibrations.model_validate(vib_dict) if vib_dict else None, remaining_block)

    @classmethod
    def _parse_frequency_payload(cls, block: str) -> dict[str, Any]:
        payload: dict[str, Any] = {}
        vibrations, _remaining = cls._extract_vibration_payload(block)
        if vibrations is not None:
            payload["vibrations"] = vibrations
        return payload

    def parse_payload(
        self, full_text: str, tree_context: G16V3ComponentTreeProtocol
    ) -> dict[str, Any]:
        if tree_context.only_extract_structure:
            self.payload = {}
            return self.payload
        self.payload = self._parse_frequency_payload(self.raw_text)
        return self.payload

    def _render_fakeg(self, **kwargs) -> str:
        rendered = _render_frequency_payload(self.payload)
        return rendered or self.raw_text

    def build_synthetic_children(self) -> tuple[G16SyntheticChildSpec, ...]:
        vibrations = self.payload.get("vibrations")
        if vibrations is None:
            return ()
        specs: list[G16SyntheticChildSpec] = []
        force_constants = getattr(vibrations, "force_constants", None)
        if force_constants is not None:
            specs.append(
                G16SyntheticChildSpec(
                    component_cls=G16V3L716ForceConstantsComponent,
                    payload={"force_constants": force_constants},
                )
            )
        vibration_modes = getattr(vibrations, "vibration_modes", None)
        if vibration_modes is not None:
            specs.append(
                G16SyntheticChildSpec(
                    component_cls=G16V3L716DiagVibComponent,
                    payload={"vibration_modes": vibration_modes},
                )
            )
        ir_intensities = getattr(vibrations, "IR_intensities", None)
        if ir_intensities is not None:
            specs.append(
                G16SyntheticChildSpec(
                    component_cls=G16V3L716IRSpectrumComponent,
                    payload={"IR_intensities": ir_intensities},
                )
            )
        for idx, vibration in enumerate(vibrations):
            specs.append(
                G16SyntheticChildSpec(
                    component_cls=G16V3L716VibrationModeComponent,
                    payload={
                        "mode_index": idx,
                        "is_imaginary": getattr(vibration, "is_imaginary", False),
                        **{
                            attr: value
                            for attr in (
                                "frequency",
                                "reduced_mass",
                                "force_constant",
                                "IR_intensity",
                                "vibration_mode",
                            )
                            if (value := getattr(vibration, attr, None)) is not None
                        },
                    },
                    node_name=f"l716.vibration.mode[{idx}]",
                )
            )
        return tuple(specs)

    def aggregate_into(self, infos: dict[str, Any], payload: Mapping[str, Any]) -> None:
        _merge_if_present(infos, "vibrations", payload.get("vibrations"))

    @classmethod
    def build_payload_from_frame(cls, frame: Any) -> dict[str, Any]:
        payload: dict[str, Any] = {}
        if (vibrations := getattr(frame, "vibrations", None)) is not None:
            payload["vibrations"] = vibrations
        if (atoms := getattr(frame, "atoms", None)) is not None:
            payload["atoms"] = atoms
        return payload


class G16V3L716FreqChunkComponent(G16V3BaseComponent):
    component_name = "l716.freq.chunkx"
    gaussian_block = "l716.freq.chunkx"


class G16V3L716ForceConstantsComponent(G16V3BaseComponent):
    component_name = "l716.forceconstants"
    gaussian_block = "l716.forceconstants"

    def _render_fakeg(self, **kwargs) -> str:
        values = _restore_raw_payload(self.payload.get("force_constants"))
        if values is None or not len(values):
            return self.raw_text
        arr = values.m if hasattr(values, "m") else values
        joined = " ".join(f"{float(v):10.4f}" for v in np.asarray(arr).reshape(-1))
        return f" Frc consts  -- {joined}"


class G16V3L716DiagVibComponent(G16V3BaseComponent):
    component_name = "l716.diagvib"
    gaussian_block = "l716.diagvib"

    def _render_fakeg(self, **kwargs) -> str:
        modes = _restore_raw_payload(self.payload.get("vibration_modes"))
        if modes is None or len(modes) == 0:
            return self.raw_text
        lines = [" Atom  AN      X      Y      Z"]
        for idx, mode in enumerate(modes[:3], start=1):
            arr = mode.m if hasattr(mode, "m") else mode
            vec = np.asarray(arr)
            first = vec[0] if vec.ndim > 1 else vec[:3]
            first = np.asarray(first).reshape(-1)
            padded = list(first[:3]) + [0.0] * max(0, 3 - len(first[:3]))
            lines.append(
                f" {idx:4d}   0 {float(padded[0]):7.2f} {float(padded[1]):7.2f} {float(padded[2]):7.2f}"
            )
        return "\n".join(lines)


class G16V3L716IRSpectrumComponent(G16V3BaseComponent):
    component_name = "l716.irspectrum"
    gaussian_block = "l716.irspectrum"

    def _render_fakeg(self, **kwargs) -> str:
        values = _restore_raw_payload(self.payload.get("IR_intensities"))
        if values is None or not len(values):
            return self.raw_text
        arr = values.m if hasattr(values, "m") else values
        joined = " ".join(f"{float(v):10.4f}" for v in np.asarray(arr).reshape(-1))
        return f" IR Inten    -- {joined}"


class G16V3L716VibrationModeComponent(G16V3BaseComponent):
    component_name = "l716.vibration.mode"
    gaussian_block = "l716.vibration.mode"

    def _render_fakeg(self, **kwargs) -> str:
        lines: list[str] = []
        if (mode_index := self.payload.get("mode_index")) is not None:
            lines.append(f"{int(mode_index) + 1:23d}")
            lines.append(f"{'A':>23}")
        if (value := _restore_raw_payload(self.payload.get("frequency"))) is not None:
            lines.append(f" Frequencies -- {_format_float(value, 4)}")
        if (value := _restore_raw_payload(self.payload.get("reduced_mass"))) is not None:
            lines.append(f" Red. masses -- {_format_float(value, 4)}")
        if (value := _restore_raw_payload(self.payload.get("force_constant"))) is not None:
            lines.append(f" Frc consts  -- {_format_float(value, 4)}")
        if (value := _restore_raw_payload(self.payload.get("IR_intensity"))) is not None:
            lines.append(f" IR Inten    -- {_format_float(value, 4)}")
        mode = _restore_raw_payload(self.payload.get("vibration_mode"))
        if mode is not None:
            arr = mode.m if hasattr(mode, "m") else mode
            vec = np.asarray(arr)
            first = vec[0] if vec.ndim > 1 else vec[:3]
            first = np.asarray(first).reshape(-1)
            padded = list(first[:3]) + [0.0] * max(0, 3 - len(first[:3]))
            lines.append(" Atom  AN      X      Y      Z")
            lines.append(
                f" {1:4d}   0 {float(padded[0]):7.2f} {float(padded[1]):7.2f} {float(padded[2]):7.2f}"
            )
        return "\n".join(lines) if lines else self.raw_text


class G16V3L716ThermochemistryComponent(G16V3BaseComponent):
    component_name = "l716.thermochemistry"
    gaussian_block = "l716.thermochemistry"
    allowed_child_component_names = (
        "l716.thermochemistry.temperature",
        "l716.thermochemistry.mass",
        "l716.thermochemistry.moi",
        "l716.thermochemistry.rotsymnum",
        "l716.thermochemistry.rottemp",
        "l716.thermochemistry.rotconsts",
        "l716.thermochemistry.vibtemp",
        "l716.thermochemistry.zpe",
        "l716.thermoprops",
        "l716.thermochemistry.energy",
        "l716.thermochemistry.enthalpy",
        "l716.thermochemistry.gibbs",
        "l716.thermochemistry.entropy",
        "l716.thermochemistry.heatcapacity",
    )

    required_frame_fields = ("thermal_informations",)
    optional_frame_fields = ("temperature", "pressure")
    boundaries = (
        G16BoundarySpec(
            component_name="l716.thermochemistry",
            start_patterns=(r"^\s*- Thermochemistry -",),
            literal_start_markers=("- Thermochemistry -",),
            contains_markers=("- Thermochemistry -",),
            end_patterns=(
                r"^\s*Center\s+Atomic\s+Forces \(Hartrees/Bohr\)",
                r"^\s*Center\s+Atomic\s+\s+\s+Forces \(Hartrees/Bohr\)",
                r"^\s*Force constants in Cartesian coordinates",
                r"^\s*Item\s+Value\s+Threshold\s+Converged\?",
                r"^\s*1[\\|]1[\\|]GINC",
                r"^\s*(Job cpu time|Elapsed time):",
                r"^\s*Normal termination of Gaussian",
                r"^\s*Leave Link",
            ),
            repeatable=True,
            boundary_source="iochem",
            match_mode="section",
            include_end=False,
        ),
    )

    @staticmethod
    def _extract_thermal_payload(block: str) -> tuple[ThermalInformations | None, str]:
        thermal_dict: dict[str, Any] = {}
        focus_content, continued_content = g16_log_patterns.THERMOCHEMISTRY_PART.split_content(
            block
        )
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
                thermal_dict[mapping[match[0]]] = float(match[1]) * atom_ureg.Unit(
                    "hartree/particle"
                )
        if matches := g16_log_patterns.THERMOCHEMISTRY_CV_S.get_matches(focus_content):
            thermal_dict["S"], thermal_dict["C_V"] = (
                float(matches[0][1]) * atom_ureg.Unit("cal/mol/K"),
                float(matches[0][2]) * atom_ureg.Unit("cal/mol/K"),
            )
        return (
            ThermalInformations.model_validate(thermal_dict) if thermal_dict else None,
            continued_content,
        )

    @classmethod
    def _parse_thermochemistry_payload(cls, raw_text: str) -> dict[str, Any]:
        payload = _temperature_and_pressure_from_block(raw_text)
        thermal, _remaining = cls._extract_thermal_payload(raw_text)
        if thermal is not None:
            payload["thermal_informations"] = thermal
        return payload

    def parse_payload(
        self, full_text: str, tree_context: G16V3ComponentTreeProtocol
    ) -> dict[str, Any]:
        if tree_context.only_extract_structure:
            self.payload = {}
            return self.payload
        self.payload = self._parse_thermochemistry_payload(self.raw_text)
        return self.payload

    def _render_fakeg(self, **kwargs) -> str:
        thermal = self.payload.get("thermal_informations")
        temperature = self.payload.get("temperature")
        pressure = self.payload.get("pressure")
        if thermal is None and temperature is None and pressure is None:
            return self.raw_text
        lines: list[str] = [" - Thermochemistry -"]
        if temperature is not None or pressure is not None:
            temperature_text = (
                f" Temperature   {_format_float(temperature, 3)} Kelvin."
                if temperature is not None
                else ""
            )
            pressure_text = (
                f"  Pressure   {_format_float(pressure, 5)} Atm." if pressure is not None else ""
            )
            lines.append(f"{temperature_text}{pressure_text}".rstrip())
        if thermal is None:
            return "\n".join(lines) if lines else self.raw_text
        for key, label in (
            ("ZPVE", "Zero-point correction="),
            ("TCE", "Thermal correction to Energy="),
            ("TCH", "Thermal correction to Enthalpy="),
            ("TCG", "Thermal correction to Gibbs Free Energy="),
        ):
            value = getattr(thermal, key, None)
            if value is not None:
                lines.append(f"{label} {_format_float(value, 6)}")
        for key, label in (
            ("U_0", "Sum of electronic and zero-point Energies="),
            ("U_T", "Sum of electronic and thermal Energies="),
            ("H_T", "Sum of electronic and thermal Enthalpies="),
            ("G_T", "Sum of electronic and thermal Free Energies="),
        ):
            value = getattr(thermal, key, None)
            if value is not None:
                lines.append(f"{label} {_format_float(value, 6)}")
        entropy = getattr(thermal, "S", None)
        heat_capacity = getattr(thermal, "C_V", None)
        if entropy is not None or heat_capacity is not None:
            thermal_energy = getattr(thermal, "TCE", None)
            lines.append(" E (Thermal)             CV                S")
            lines.append(
                " Total"
                f"{_format_float(thermal_energy, 3) if thermal_energy is not None else '0.000':>13}"
                f"{_format_float(heat_capacity, 3) if heat_capacity is not None else '0.000':>18}"
                f"{_format_float(entropy, 3) if entropy is not None else '0.000':>17}"
            )
        lines.append(" Rotational           0.000000D+00    0.0000    0.0000")
        return "\n".join(lines) if lines else self.raw_text

    def build_synthetic_children(self) -> tuple[G16SyntheticChildSpec, ...]:
        specs: list[G16SyntheticChildSpec] = []
        if self.payload.get("temperature") is not None or self.payload.get("pressure") is not None:
            specs.append(
                G16SyntheticChildSpec(
                    component_cls=G16V3L716ThermochemistryTemperatureComponent,
                    payload={
                        key: self.payload[key]
                        for key in ("temperature", "pressure")
                        if self.payload.get(key) is not None
                    },
                )
            )
        thermal_info = self.payload.get("thermal_informations")
        if thermal_info is not None:
            if getattr(thermal_info, "ZPVE", None) is not None:
                specs.append(
                    G16SyntheticChildSpec(
                        component_cls=G16V3L716ThermochemistryZPEComponent,
                        payload={"ZPVE": thermal_info.ZPVE},
                    )
                )
            specs.append(
                G16SyntheticChildSpec(
                    component_cls=G16V3L716ThermoPropsComponent,
                    payload={
                        attr: value
                        for attr in (
                            "U_0",
                            "U_T",
                            "TCE",
                            "H_T",
                            "TCH",
                            "G_T",
                            "TCG",
                            "S",
                            "C_V",
                        )
                        if (value := getattr(thermal_info, attr, None)) is not None
                    },
                )
            )
            for component_cls, attrs in (
                (G16V3L716ThermochemistryEnergyComponent, ("U_0", "U_T", "TCE")),
                (G16V3L716ThermochemistryEnthalpyComponent, ("H_T", "TCH")),
                (G16V3L716ThermochemistryGibbsComponent, ("G_T", "TCG")),
                (G16V3L716ThermochemistryEntropyComponent, ("S",)),
                (G16V3L716ThermochemistryHeatCapacityComponent, ("C_V",)),
            ):
                payload = {
                    attr: value
                    for attr in attrs
                    if (value := getattr(thermal_info, attr, None)) is not None
                }
                if payload:
                    specs.append(
                        G16SyntheticChildSpec(component_cls=component_cls, payload=payload)
                    )
        if matches := g16_log_patterns.MOLECULAR_MASS.get_matches(self.raw_text):
            specs.append(
                G16SyntheticChildSpec(
                    component_cls=G16V3L716ThermochemistryMassComponent,
                    payload={"molecular_mass": float(matches[0][0]) * atom_ureg.amu},
                )
            )
        if matches := g16_log_patterns.MOMENTS_OF_INERTIA.get_matches(self.raw_text):
            values = [float(v) for v in re.findall(r"-?\d+\.\d+", matches[0][0])]
            if values:
                specs.append(
                    G16SyntheticChildSpec(
                        component_cls=G16V3L716ThermochemistryMOIComponent,
                        payload={
                            "moments_of_inertia": np.array(values)
                            * atom_ureg.amu
                            * atom_ureg.bohr**2
                        },
                    )
                )
        if matches := g16_log_patterns.ROTATIONAL_SYMMETRY_NUMBER.get_matches(self.raw_text):
            specs.append(
                G16SyntheticChildSpec(
                    component_cls=G16V3L716ThermochemistryRotSymNumComponent,
                    payload={"rotational_symmetry_number": int(matches[0][0])},
                )
            )
        if matches := g16_log_patterns.ROTATIONAL_TEMPERATURE.get_matches(self.raw_text):
            specs.append(
                G16SyntheticChildSpec(
                    component_cls=G16V3L716ThermochemistryRotTempComponent,
                    payload={
                        "rotational_temperatures": np.array(list(map(float, matches[0])))
                        * atom_ureg.K
                    },
                )
            )
        if matches := g16_log_patterns.ROTATIONAL_CONST_IN_FREQUENCY_ANALYSIS.get_matches(
            self.raw_text
        ):
            specs.append(
                G16SyntheticChildSpec(
                    component_cls=G16V3L716ThermochemistryRotConstsComponent,
                    payload={
                        "rotational_constants": np.array(list(map(float, matches[0])))
                        * atom_ureg.gigahertz
                    },
                )
            )
        if matches := g16_log_patterns.VIBRATIONAL_TEMPERATURE.get_matches(self.raw_text):
            numeric_tokens = [
                float(token) for token in matches[0][1:] if token is not None and str(token).strip()
            ]
            if numeric_tokens:
                specs.append(
                    G16SyntheticChildSpec(
                        component_cls=G16V3L716ThermochemistryVibTempComponent,
                        payload={
                            "vibrational_temperatures": np.array(numeric_tokens) * atom_ureg.K
                        },
                    )
                )
        return tuple(specs)

    def aggregate_into(self, infos: dict[str, Any], payload: Mapping[str, Any]) -> None:
        _merge_if_present(infos, "thermal_informations", payload.get("thermal_informations"))

    @classmethod
    def build_payload_from_frame(cls, frame: Any) -> dict[str, Any]:
        payload: dict[str, Any] = {}
        if (thermal := getattr(frame, "thermal_informations", None)) is not None:
            payload["thermal_informations"] = thermal
        if (temperature := getattr(frame, "temperature", None)) is not None:
            payload["temperature"] = temperature
        if (pressure := getattr(frame, "pressure", None)) is not None:
            payload["pressure"] = pressure
        return payload


class G16V3L716ThermochemistryTemperatureComponent(G16V3BaseComponent):
    component_name = "l716.thermochemistry.temperature"
    gaussian_block = "l716.thermochemistry.temperature"

    def _render_fakeg(self, **kwargs) -> str:
        temperature = _restore_raw_payload(self.payload.get("temperature"))
        pressure = _restore_raw_payload(self.payload.get("pressure"))
        if temperature is None and pressure is None:
            return self.raw_text
        if temperature is not None and pressure is not None:
            return (
                f" Temperature   {_format_float(temperature, 3)} Kelvin.  "
                f"Pressure   {_format_float(pressure, 5)} Atm."
            )
        if temperature is not None:
            return f" Temperature   {_format_float(temperature, 3)} Kelvin."
        return f" Pressure   {_format_float(pressure, 5)} Atm."


class G16V3L716ThermochemistryMassComponent(G16V3BaseComponent):
    component_name = "l716.thermochemistry.mass"
    gaussian_block = "l716.thermochemistry.mass"

    def _render_fakeg(self, **kwargs) -> str:
        value = _restore_raw_payload(self.payload.get("molecular_mass"))
        if value is not None:
            return f" Molecular mass: {_format_float(value, 4)} amu."
        return self.raw_text


class G16V3L716ThermochemistryMOIComponent(G16V3BaseComponent):
    component_name = "l716.thermochemistry.moi"
    gaussian_block = "l716.thermochemistry.moi"

    def _render_fakeg(self, **kwargs) -> str:
        values = _restore_raw_payload(self.payload.get("moments_of_inertia"))
        if values is None or not len(values):
            return self.raw_text
        arr = values.m if hasattr(values, "m") else values
        joined = " ".join(f"{float(v):10.4f}" for v in np.asarray(arr).reshape(-1))
        return f" Eigenvalues -- {joined}"


class G16V3L716ThermochemistryRotSymNumComponent(G16V3BaseComponent):
    component_name = "l716.thermochemistry.rotsymnum"
    gaussian_block = "l716.thermochemistry.rotsymnum"

    def _render_fakeg(self, **kwargs) -> str:
        value = _restore_raw_payload(self.payload.get("rotational_symmetry_number"))
        if value is not None:
            return f" Rotational symmetry number {int(value)}."
        return self.raw_text


class G16V3L716ThermochemistryRotTempComponent(G16V3BaseComponent):
    component_name = "l716.thermochemistry.rottemp"
    gaussian_block = "l716.thermochemistry.rottemp"

    def _render_fakeg(self, **kwargs) -> str:
        rotational_temperatures = _restore_raw_payload(self.payload.get("rotational_temperatures"))
        if rotational_temperatures is None:
            return self.raw_text
        arr = (
            rotational_temperatures.m
            if hasattr(rotational_temperatures, "m")
            else rotational_temperatures
        )
        joined = "".join(f" {float(v):10.4f}" for v in np.asarray(arr).reshape(-1))
        return f" Rotational temperatures (Kelvin){joined}"


class G16V3L716ThermochemistryRotConstsComponent(G16V3BaseComponent):
    component_name = "l716.thermochemistry.rotconsts"
    gaussian_block = "l716.thermochemistry.rotconsts"

    def _render_fakeg(self, **kwargs) -> str:
        rotational_constants = _restore_raw_payload(self.payload.get("rotational_constants"))
        if rotational_constants is None:
            return self.raw_text
        arr = rotational_constants.m if hasattr(rotational_constants, "m") else rotational_constants
        joined = "".join(f" {float(v):10.4f}" for v in np.asarray(arr).reshape(-1))
        return f" Rotational constants (GHZ):{joined}"


class G16V3L716ThermochemistryVibTempComponent(G16V3BaseComponent):
    component_name = "l716.thermochemistry.vibtemp"
    gaussian_block = "l716.thermochemistry.vibtemp"

    def _render_fakeg(self, **kwargs) -> str:
        vibrational_temperatures = _restore_raw_payload(
            self.payload.get("vibrational_temperatures")
        )
        if vibrational_temperatures is None:
            return self.raw_text
        arr = (
            vibrational_temperatures.m
            if hasattr(vibrational_temperatures, "m")
            else vibrational_temperatures
        )
        joined = "".join(f" {float(v):10.4f}" for v in np.asarray(arr).reshape(-1))
        return f" Vibrational temperatures:{joined}"


class G16V3L716ThermochemistryZPEComponent(G16V3BaseComponent):
    component_name = "l716.thermochemistry.zpe"
    gaussian_block = "l716.thermochemistry.zpe"

    def _render_fakeg(self, **kwargs) -> str:
        zpe = _restore_raw_payload(self.payload.get("ZPVE"))
        if zpe is not None:
            return f"Zero-point correction= {_format_float(zpe, 6)} (Hartree/Particle)"
        return self.raw_text


class G16V3L716ThermoPropsComponent(G16V3BaseComponent):
    component_name = "l716.thermoprops"
    gaussian_block = "l716.thermoprops"


class G16V3L716ThermochemistryEnergyComponent(G16V3BaseComponent):
    component_name = "l716.thermochemistry.energy"
    gaussian_block = "l716.thermochemistry.energy"

    def _render_fakeg(self, **kwargs) -> str:
        lines = []
        value = _restore_raw_payload(self.payload.get("TCE"))
        if value is not None:
            lines.append(f"Thermal correction to Energy= {_format_float(value, 6)}")
        value = _restore_raw_payload(self.payload.get("U_0"))
        if value is not None:
            lines.append(f"Sum of electronic and zero-point Energies= {_format_float(value, 6)}")
        value = _restore_raw_payload(self.payload.get("U_T"))
        if value is not None:
            lines.append(f"Sum of electronic and thermal Energies= {_format_float(value, 6)}")
        return "\n".join(lines) if lines else self.raw_text


class G16V3L716ThermochemistryEnthalpyComponent(G16V3BaseComponent):
    component_name = "l716.thermochemistry.enthalpy"
    gaussian_block = "l716.thermochemistry.enthalpy"

    def _render_fakeg(self, **kwargs) -> str:
        lines = []
        value = _restore_raw_payload(self.payload.get("TCH"))
        if value is not None:
            lines.append(f"Thermal correction to Enthalpy= {_format_float(value, 6)}")
        value = _restore_raw_payload(self.payload.get("H_T"))
        if value is not None:
            lines.append(f"Sum of electronic and thermal Enthalpies= {_format_float(value, 6)}")
        return "\n".join(lines) if lines else self.raw_text


class G16V3L716ThermochemistryGibbsComponent(G16V3BaseComponent):
    component_name = "l716.thermochemistry.gibbs"
    gaussian_block = "l716.thermochemistry.gibbs"

    def _render_fakeg(self, **kwargs) -> str:
        lines = []
        value = _restore_raw_payload(self.payload.get("TCG"))
        if value is not None:
            lines.append(f"Thermal correction to Gibbs Free Energy= {_format_float(value, 6)}")
        value = _restore_raw_payload(self.payload.get("G_T"))
        if value is not None:
            lines.append(f"Sum of electronic and thermal Free Energies= {_format_float(value, 6)}")
        return "\n".join(lines) if lines else self.raw_text


class G16V3L716ThermochemistryEntropyComponent(G16V3BaseComponent):
    component_name = "l716.thermochemistry.entropy"
    gaussian_block = "l716.thermochemistry.entropy"

    def _render_fakeg(self, **kwargs) -> str:
        value = _restore_raw_payload(self.payload.get("S"))
        if value is not None:
            return f"Total Entropy= {_format_float(value, 4)} cal/mol-K"
        return self.raw_text


class G16V3L716ThermochemistryHeatCapacityComponent(G16V3BaseComponent):
    component_name = "l716.thermochemistry.heatcapacity"
    gaussian_block = "l716.thermochemistry.heatcapacity"

    def _render_fakeg(self, **kwargs) -> str:
        value = _restore_raw_payload(self.payload.get("C_V"))
        if value is not None:
            return f"Total Heat Capacity= {_format_float(value, 4)} cal/mol-K"
        return self.raw_text


class G16V3L716PolarizabilityComponent(G16V3BaseComponent):
    component_name = "l716.polarizability"
    gaussian_block = "l716.polarizability"
    allowed_child_component_names = (
        "l716.dipole",
        "l716.polarizability.detail",
    )
    required_frame_fields = ("polarizability",)
    boundaries = (
        G16BoundarySpec(
            component_name="l716.polarizability",
            start_patterns=(r"^\s*Isotropic polarizability for W=",),
            literal_start_markers=("Isotropic polarizability for W=",),
            contains_markers=("Isotropic polarizability for W=",),
            repeatable=True,
            boundary_source="molop",
            match_mode="line",
        ),
        G16BoundarySpec(
            component_name="l716.polarizability",
            start_patterns=(r"^\s*Electric dipole moment \(input orientation\):",),
            literal_start_markers=("Electric dipole moment (input orientation):",),
            contains_markers=("Electric dipole moment (input orientation):",),
            end_patterns=(
                r"^\s*Center\s+Atomic\s+Forces \(Hartrees/Bohr\)",
                r"^\s*Center\s+Atomic\s+\s+\s+Forces \(Hartrees/Bohr\)",
                r"^\s*Force constants in Cartesian coordinates",
                r"Harmonic frequencies \(cm\*\*-1\)",
                r"^\s*- Thermochemistry -",
                r"^\s*1[\\|]1[\\|]GINC",
                r"^\s*(Job cpu time|Elapsed time):",
                r"^\s*Normal termination of Gaussian",
                r"^\s*Leave Link",
            ),
            repeatable=True,
            boundary_source="molop",
            match_mode="section",
            include_end=False,
        ),
    )

    @staticmethod
    def _extract_polarizability_payload(block: str) -> tuple[Polarizability | None, str]:
        focus_content, continued_content = g16_log_patterns.ISOTROPIC_POLARIZABILITY.split_content(
            block
        )
        if focus_content == "":
            return None, block
        if matches := g16_log_patterns.ISOTROPIC_POLARIZABILITY.get_matches(focus_content):
            return (
                Polarizability(isotropic_polarizability=float(matches[0][0]) * atom_ureg.bohr**3),
                continued_content,
            )
        return None, continued_content

    @classmethod
    def _parse_polarizability_payload(cls, raw_text: str) -> dict[str, Any]:
        payload: dict[str, Any] = {}
        polarizability, _remaining = cls._extract_polarizability_payload(raw_text)
        if polarizability is not None:
            payload.update(_payload_mapping(polarizability))
        return payload

    def parse_payload(
        self, full_text: str, tree_context: G16V3ComponentTreeProtocol
    ) -> dict[str, Any]:
        if tree_context.only_extract_structure:
            self.payload = {}
            return self.payload
        self.payload = self._parse_polarizability_payload(self.raw_text)
        return self.payload

    def build_synthetic_children(self) -> tuple[G16SyntheticChildSpec, ...]:
        specs: list[G16SyntheticChildSpec] = []
        if self.payload.get("dipole") is not None:
            specs.append(
                G16SyntheticChildSpec(
                    component_cls=G16V3L716DipoleComponent,
                    payload={"dipole": self.payload.get("dipole")},
                )
            )
        detail_payload = {
            key: value
            for key, value in self.payload.items()
            if key != "dipole" and value is not None
        }
        if detail_payload:
            specs.append(
                G16SyntheticChildSpec(
                    component_cls=G16V3L716PolarizabilityDetailComponent,
                    payload=detail_payload,
                )
            )
        return tuple(specs)

    def aggregate_into(self, infos: dict[str, Any], payload: Mapping[str, Any]) -> None:
        if payload:
            _merge_if_present(infos, "polarizability", payload)

    @classmethod
    def build_payload_from_frame(cls, frame: Any) -> dict[str, Any]:
        if (polarizability := getattr(frame, "polarizability", None)) is not None:
            return dict(_payload_mapping(polarizability))
        return {}


class G16V3L716DipoleComponent(G16V3BaseComponent):
    component_name = "l716.dipole"
    gaussian_block = "l716.dipole"


class G16V3L716PolarizabilityDetailComponent(G16V3BaseComponent):
    component_name = "l716.polarizability.detail"
    gaussian_block = "l716.polarizability.detail"


class G16V3L716ForcesComponent(G16V3BaseComponent):
    component_name = "l716.forces"
    gaussian_block = "l716.forces"
    boundaries = (
        G16BoundarySpec(
            component_name="l716.forces",
            start_patterns=(r"^\s*Center\s+Atomic\s+Forces \(Hartrees/Bohr\)",),
            contains_markers=("Forces (Hartrees/Bohr)",),
            repeatable=True,
            boundary_source="molop",
            match_mode="molop_pattern",
            molop_patterns=(FORCES_IN_CARTESIAN_V2,),
        ),
    )
    required_frame_fields = ("forces",)

    @staticmethod
    def _extract_forces_payload(block: str) -> tuple[NumpyQuantity | None, str]:
        focus_content, continued_content = g16_log_patterns.FORCES_IN_CARTESIAN.split_content(block)
        if focus_content == "":
            return None, block
        if matches := g16_log_patterns.FORCES_IN_CARTESIAN.get_matches(focus_content):
            return (
                np.array([list(map(float, match)) for match in matches])
                * atom_ureg.hartree
                / atom_ureg.bohr,
                continued_content,
            )
        return None, continued_content

    @classmethod
    def _parse_forces_payload(cls, raw_text: str) -> dict[str, Any]:
        payload: dict[str, Any] = {}
        forces, _remaining = cls._extract_forces_payload(raw_text)
        if forces is not None:
            payload["forces"] = forces
        return payload

    def parse_payload(
        self, full_text: str, tree_context: G16V3ComponentTreeProtocol
    ) -> dict[str, Any]:
        if tree_context.only_extract_structure:
            self.payload = {}
            return self.payload
        self.payload = self._parse_forces_payload(self.raw_text)
        return self.payload

    def aggregate_into(self, infos: dict[str, Any], payload: Mapping[str, Any]) -> None:
        if payload.get("forces") is not None:
            infos["forces"] = payload["forces"]

    @classmethod
    def build_payload_from_frame(cls, frame: Any) -> dict[str, Any]:
        if (forces := getattr(frame, "forces", None)) is not None:
            return {"forces": forces}
        return {}


class G16V3L716SecondDerivComponent(G16V3BaseComponent):
    component_name = "l716.secondderiv"
    gaussian_block = "l716.secondderiv"
    boundaries = (
        G16BoundarySpec(
            component_name="l716.secondderiv",
            start_patterns=(r"^\s*Force constants in Cartesian coordinates",),
            literal_start_markers=("Force constants in Cartesian coordinates",),
            contains_markers=("Force constants in Cartesian coordinates",),
            repeatable=True,
            boundary_source="molop",
            match_mode="molop_pattern",
            molop_patterns=(HESSIAN_IN_CARTESIAN_V2,),
        ),
    )
    required_frame_fields = ("hessian",)

    @staticmethod
    def _extract_hessian_payload(block: str) -> tuple[NumpyQuantity | None, str]:
        focus_content, continued_content = g16_log_patterns.HESSIAN_IN_CARTESIAN.split_content(
            block
        )
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
                / atom_ureg.bohr**2,
                continued_content,
            )
        return None, continued_content

    @classmethod
    def _parse_hessian_payload(cls, raw_text: str) -> dict[str, Any]:
        payload: dict[str, Any] = {}
        hessian, _remaining = cls._extract_hessian_payload(raw_text)
        if hessian is not None:
            payload["hessian"] = hessian
        return payload

    def parse_payload(
        self, full_text: str, tree_context: G16V3ComponentTreeProtocol
    ) -> dict[str, Any]:
        if tree_context.only_extract_structure:
            self.payload = {}
            return self.payload
        self.payload = self._parse_hessian_payload(self.raw_text)
        return self.payload

    def aggregate_into(self, infos: dict[str, Any], payload: Mapping[str, Any]) -> None:
        if payload.get("hessian") is not None:
            infos["hessian"] = payload["hessian"]

    @classmethod
    def build_payload_from_frame(cls, frame: Any) -> dict[str, Any]:
        if (hessian := getattr(frame, "hessian", None)) is not None:
            return {"hessian": hessian}
        return {}


class G16V3L103OptimizationComponent(G16V3BaseComponent):
    component_name = "l103.optimization"
    gaussian_block = "l103.optimization"
    boundaries = (
        G16BoundarySpec(
            component_name="l103.optimization",
            start_patterns=(r"^\s*Item\s+Value\s+Threshold\s+Converged\?",),
            contains_markers=("Item               Value     Threshold  Converged?",),
            end_patterns=(
                r"^\s*Input orientation:",
                r"^\s*Standard orientation:",
                r"^\s*Rotational constants \(GHZ\):",
                r"^\s*Population analysis using the (SCF|CC) [Dd]ensity\.",
                r"Harmonic frequencies \(cm\*\*-1\)",
                r"^\s*- Thermochemistry -",
                r"^\s*Center\s+Atomic\s+Forces \(Hartrees/Bohr\)",
                r"^\s*Center\s+Atomic\s+\s+\s+Forces \(Hartrees/Bohr\)",
                r"^\s*Force constants in Cartesian coordinates",
                r"^\s*1[\\|]1[\\|]GINC",
                r"^\s*(Job cpu time|Elapsed time):",
                r"^\s*Normal termination of Gaussian",
                r"^\s*Leave Link",
            ),
            repeatable=True,
            boundary_source="molop",
            match_mode="section",
            include_end=False,
        ),
    )
    required_frame_fields = ("geometry_optimization_status",)

    @staticmethod
    def _extract_optimization_payload(block: str) -> tuple[GeometryOptimizationStatus | None, str]:
        berny_dict: dict[str, float | bool] = {}
        start_index = block.find("Item               Value     Threshold  Converged?")
        end_index = block.find(
            "GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad",
            start_index,
        )
        if start_index == -1 or end_index == -1:
            focus_content, continued_content = (
                g16_log_patterns.BERNY_STATE_MAJOR_PART.split_content(block)
            )
            if focus_content == "":
                focus_content, continued_content = (
                    g16_log_patterns.BERNY_STATE_BACKUP_PART.split_content(block)
                )
        else:
            focus_content = block[start_index:end_index]
            continued_content = block[end_index:]
        if focus_content == "":
            return None, block
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
        return (
            GeometryOptimizationStatus.model_validate(berny_dict) if berny_dict else None,
            continued_content,
        )

    @classmethod
    def _parse_optimization_payload(cls, raw_text: str) -> dict[str, Any]:
        payload: dict[str, Any] = {}
        berny, _remaining = cls._extract_optimization_payload(raw_text)
        if berny is not None:
            payload["geometry_optimization_status"] = berny
        return payload

    def parse_payload(
        self, full_text: str, tree_context: G16V3ComponentTreeProtocol
    ) -> dict[str, Any]:
        if tree_context.only_extract_structure:
            self.payload = {}
            return self.payload
        self.payload = self._parse_optimization_payload(self.raw_text)
        return self.payload

    def aggregate_into(self, infos: dict[str, Any], payload: Mapping[str, Any]) -> None:
        if payload.get("geometry_optimization_status") is not None:
            infos["geometry_optimization_status"] = payload["geometry_optimization_status"]

    @classmethod
    def build_payload_from_frame(cls, frame: Any) -> dict[str, Any]:
        if (opt := getattr(frame, "geometry_optimization_status", None)) is not None:
            return {"geometry_optimization_status": opt}
        return {}


class G16V3L9999ArchiveComponent(G16V3BaseComponent):
    component_name = "l9999.archive"
    gaussian_block = "l9999.archive"
    allowed_child_component_names = (
        "l9999.archive.metadata",
        "l9999.archive.energies",
        "l9999.archive.thermochemistry",
        "l9999.archive.polarizability",
        "l9999.archive.hessian",
    )
    required_frame_fields = ("job_type",)
    optional_frame_fields = (
        "functional",
        "basis_set",
        "keywords",
        "title_card",
        "charge",
        "multiplicity",
        "qm_software_version",
        "atoms",
        "coords",
        "energies",
        "thermal_informations",
        "polarizability",
        "hessian",
    )
    boundaries = (
        G16BoundarySpec(
            component_name="l9999.archive",
            start_patterns=(r"^\s*1[\\|]1[\\|]GINC",),
            contains_markers=("1\\1\\GINC",),
            repeatable=True,
            boundary_source="iochem",
            match_mode="molop_pattern",
            molop_patterns=(ARCHIVE_TAIL_V2,),
        ),
    )

    @staticmethod
    def _extract_archive_metadata(block: str) -> tuple[dict[str, Any], str]:
        return parse_archive_tail(block, include_coords=True)

    @staticmethod
    def _extract_archive_energies(block: str) -> tuple[Energies | None, str]:
        energy_dict: dict[str, Any] = {}
        focus_content, remaining_block = g16_log_patterns.ENERGIES_IN_ARCHIVE_TAIL.split_content(
            block
        )
        if focus_content == "":
            return None, block
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
        return (Energies.model_validate(energy_dict) if energy_dict else None, remaining_block)

    @staticmethod
    def _extract_archive_thermal_infos(block: str) -> tuple[ThermalInformations | None, str]:
        thermal_dict: dict[str, Any] = {}
        thermal_mapping = {
            "ZeroPoint": "ZPVE",
            "Thermal": "TCE",
            "ETot": "U_T",
            "HTot": "H_T",
            "GTot": "G_T",
        }
        focus_content, remaining_block = (
            g16_log_patterns.THERMOCHEMISTRY_IN_ARCHIVE_TAIL.split_content(block)
        )
        if focus_content == "":
            return None, block
        if matches := g16_log_patterns.THERMOCHEMISTRY_IN_ARCHIVE_TAIL.get_matches(focus_content):
            for match in matches:
                if match[0] in thermal_mapping:
                    thermal_dict[thermal_mapping[match[0]]] = float(match[1]) * atom_ureg.Unit(
                        "hartree/particle"
                    )
        return (
            ThermalInformations.model_validate(thermal_dict) if thermal_dict else None,
            remaining_block,
        )

    @staticmethod
    def _extract_archive_polarizability(block: str) -> dict[str, Any] | None:
        polarizability_dict: dict[str, Any] = {}
        if matches := g16_log_patterns.DIPOLE_IN_ARCHIVE_TAIL.get_matches(block):
            polarizability_dict["dipole"] = (
                np.array([float(match[1]) for match in matches[1:]]) * atom_ureg.debye
            )
        if matches := g16_log_patterns.POLAR_IN_ARCHIVE_TAIL.get_matches(block):
            polarizability_dict["polarizability_tensor"] = (
                np.array([float(match[1]) for match in matches[1:]]) * atom_ureg.bohr**3
            )
        if matches := g16_log_patterns.QUADRUPOLE_IN_ARCHIVE_TAIL.get_matches(block):
            polarizability_dict["quadrupole"] = (
                np.array([float(match[1]) for match in matches[1:]]) * atom_ureg.bohr**3
            )
        return polarizability_dict or None

    @staticmethod
    def _extract_archive_hessian(block: str) -> tuple[NumpyQuantity | None, str]:
        focus_content, remaining_block = g16_log_patterns.HESSIAN_IN_ARCHIVE_TAIL.split_content(
            block
        )
        if matches := g16_log_patterns.HESSIAN_IN_ARCHIVE_TAIL.get_matches(focus_content):
            return (
                fill_symmetric_matrix(np.array([float(match[0]) for match in matches]))
                * atom_ureg.hartree
                / atom_ureg.bohr**2,
                remaining_block,
            )
        return None, remaining_block

    @classmethod
    def _parse_archive_payload(cls, raw_text: str, *, structure_only: bool) -> dict[str, Any]:
        metadata, remaining_tail = cls._extract_archive_metadata(raw_text)
        payload = dict(metadata)
        if structure_only:
            return payload
        energies, _remaining_tail = cls._extract_archive_energies(remaining_tail)
        if energies is not None:
            payload["energies"] = energies
        thermal, _remaining_tail = cls._extract_archive_thermal_infos(remaining_tail)
        if thermal is not None:
            payload["thermal_informations"] = thermal
        if tail_polarizability := cls._extract_archive_polarizability(remaining_tail):
            payload["polarizability"] = tail_polarizability
        hessian, _remaining_tail = cls._extract_archive_hessian(remaining_tail)
        if hessian is not None:
            payload["hessian"] = hessian
        return payload

    def parse_payload(
        self, full_text: str, tree_context: G16V3ComponentTreeProtocol
    ) -> dict[str, Any]:
        self.payload = self._parse_archive_payload(
            self.raw_text, structure_only=tree_context.only_extract_structure
        )
        return self.payload

    def aggregate_into(self, infos: dict[str, Any], payload: Mapping[str, Any]) -> None:
        for key in (
            "job_type",
            "functional",
            "basis_set",
            "keywords",
            "title_card",
            "charge",
            "multiplicity",
            "qm_software_version",
        ):
            if key not in infos and payload.get(key) is not None:
                infos[key] = payload[key]
        if payload.get("atoms") and "atoms" not in infos:
            infos["atoms"] = payload["atoms"]
        if payload.get("coords") is not None and "coords" not in infos:
            infos["coords"] = payload["coords"]
        _merge_if_present_no_force(infos, "energies", payload.get("energies"))
        _merge_if_present_no_force(
            infos, "thermal_informations", payload.get("thermal_informations")
        )
        _merge_if_present_no_force(infos, "polarizability", payload.get("polarizability"))
        if payload.get("hessian") is not None and "hessian" not in infos:
            infos["hessian"] = payload["hessian"]

    @classmethod
    def build_payload_from_frame(cls, frame: Any) -> dict[str, Any]:
        payload: dict[str, Any] = {}
        for key in (
            "job_type",
            "functional",
            "basis_set",
            "keywords",
            "title_card",
            "charge",
            "multiplicity",
            "qm_software_version",
            "atoms",
            "coords",
            "energies",
            "thermal_informations",
            "polarizability",
            "hessian",
        ):
            if (value := getattr(frame, key, None)) is not None:
                payload[key] = value
        return payload


class G16V3L9999FinalComponent(G16V3BaseComponent):
    component_name = "l9999.final"
    gaussian_block = "l9999.final"
    allowed_child_component_names = ()
    required_frame_fields = ("status",)
    boundaries = (
        G16BoundarySpec(
            component_name="l9999.final",
            start_patterns=(r"^\s*Normal termination of Gaussian",),
            literal_start_markers=("Normal termination of Gaussian",),
            contains_markers=("Normal termination of Gaussian",),
            repeatable=True,
            overlap_policy="allow_shared_end",
            boundary_source="iochem",
            match_mode="line",
        ),
    )

    @staticmethod
    def _extract_termination_status(block: str) -> dict[str, Any]:
        return _parse_termination_status_payload_impl(block)

    def parse_payload(
        self, full_text: str, tree_context: G16V3ComponentTreeProtocol
    ) -> dict[str, Any]:
        self.payload = self._extract_termination_status(self.raw_text)
        return self.payload

    def aggregate_into(self, infos: dict[str, Any], payload: Mapping[str, Any]) -> None:
        return None

    @classmethod
    def build_payload_from_frame(cls, frame: Any) -> dict[str, Any]:
        if (status := getattr(frame, "status", None)) is not None:
            return {"status": status}
        return {}


class G16V3JobCPUComponent(G16V3BaseComponent):
    component_name = "jobcpu"
    gaussian_block = "jobcpu"
    allowed_child_component_names = ()
    required_frame_fields = ("running_time",)
    boundaries = (
        G16BoundarySpec(
            component_name="jobcpu",
            start_patterns=(r"^\s*(Job cpu time|Elapsed time):",),
            literal_start_markers=("Job cpu time:", "Elapsed time:"),
            contains_markers=("Job cpu time:", "Elapsed time:"),
            end_patterns=(r"^\s*Normal termination of Gaussian.*",),
            repeatable=True,
            overlap_policy="allow_shared_end",
            boundary_source="iochem",
            match_mode="section",
            include_end=True,
        ),
    )

    @staticmethod
    def _extract_running_time(block: str) -> PlainQuantity | None:
        return _parse_running_time_impl(block)

    @classmethod
    def _parse_jobcpu_component_payload(cls, block: str) -> dict[str, Any]:
        payload = _parse_jobcpu_payload(block)
        payload.update(
            {
                key: value
                for key, value in G16V3L9999FinalComponent._extract_termination_status(
                    block
                ).items()
                if key not in {"termination_line"} or value
            }
        )
        return payload

    def parse_payload(
        self, full_text: str, tree_context: G16V3ComponentTreeProtocol
    ) -> dict[str, Any]:
        self.payload = self._parse_jobcpu_component_payload(self.raw_text)
        return self.payload

    def _render_fakeg(self, **kwargs) -> str:
        running_time = self.payload.get("job_cpu_time") or self.payload.get("running_time")
        if running_time is None:
            return self.raw_text
        seconds = float(
            running_time.to("second").m if hasattr(running_time, "to") else running_time
        )
        return _format_job_time_seconds(seconds)

    @classmethod
    def build_payload_from_frame(cls, frame: Any) -> dict[str, Any]:
        if (running_time := getattr(frame, "running_time", None)) is not None:
            return {"running_time": running_time, "job_cpu_time": running_time}
        return {}


class G16V3LinkEnterComponent(G16V3BaseComponent):
    component_name = "link.enter"
    gaussian_block = "link.enter"
    allowed_child_component_names = ()
    boundaries = (
        G16BoundarySpec(
            component_name="link.enter",
            start_patterns=(r"^\s*\(Enter\s+.+\)\s*$",),
            literal_start_markers=("(Enter ",),
            contains_markers=("(Enter ",),
            repeatable=True,
            boundary_source="iochem",
            match_mode="line",
        ),
    )


class G16V3LeaveComponent(G16V3BaseComponent):
    component_name = "leave"
    gaussian_block = "leave"
    allowed_child_component_names = ()
    required_frame_fields = ("running_time",)
    boundaries = (
        G16BoundarySpec(
            component_name="leave",
            start_patterns=(r"^\s*Leave Link",),
            literal_start_markers=("Leave Link",),
            contains_markers=("Leave Link",),
            repeatable=True,
            boundary_source="iochem",
            match_mode="line",
        ),
    )

    @staticmethod
    def _extract_leave_running_time(block: str) -> PlainQuantity | None:
        return G16V3JobCPUComponent._extract_running_time(block)

    def parse_payload(
        self, full_text: str, tree_context: G16V3ComponentTreeProtocol
    ) -> dict[str, Any]:
        payload: dict[str, Any] = {}
        if running_time := self._extract_leave_running_time(self.raw_text):
            payload["running_time"] = running_time
        self.payload = payload
        return self.payload

    def aggregate_into(self, infos: dict[str, Any], payload: Mapping[str, Any]) -> None:
        if payload.get("running_time") is not None:
            running_time = payload["running_time"]
            if hasattr(running_time, "to"):
                infos.setdefault("__running_time_accumulated_seconds__", 0.0)
                infos["__running_time_accumulated_seconds__"] += float(running_time.to("second").m)
            else:
                infos.setdefault("__running_time_accumulated_seconds__", 0.0)
                infos["__running_time_accumulated_seconds__"] += float(running_time)

    @classmethod
    def build_payload_from_frame(cls, frame: Any) -> dict[str, Any]:
        if (running_time := getattr(frame, "running_time", None)) is not None:
            return {"running_time": running_time}
        return {}


G16LOG_V3_COMPONENT_REGISTRY: tuple[type[G16V3BaseComponent], ...] = (
    G16V3L1HeaderComponent,
    G16V3L101TitleComponent,
    G16V3L1KeywordsComponent,
    G16V3L202OrientComponent,
    G16V3L202RotConstComponent,
    G16V3L502CycleComponent,
    G16V3L601PopAnalComponent,
    G16V3L716FreqComponent,
    G16V3L716ThermochemistryComponent,
    G16V3L716PolarizabilityComponent,
    G16V3L716ForcesComponent,
    G16V3L716SecondDerivComponent,
    G16V3L103OptimizationComponent,
    G16V3L9999ArchiveComponent,
    G16V3L9999FinalComponent,
    G16V3JobCPUComponent,
    G16V3LinkEnterComponent,
    G16V3LeaveComponent,
)


def get_g16log_v3_component_classes() -> tuple[type[G16V3BaseComponent], ...]:
    return G16LOG_V3_COMPONENT_REGISTRY


def get_g16log_v3_candidate_component_classes(
    full_text: str,
    *,
    context: G16ScanContext | None = None,
) -> tuple[type[G16V3BaseComponent], ...]:
    selected: list[type[G16V3BaseComponent]] = []
    for component_cls in G16LOG_V3_COMPONENT_REGISTRY:
        markers = component_cls.boundary_markers()
        if markers and not (
            context.has_any_markers(markers)
            if context is not None
            else any(marker in full_text for marker in markers)
        ):
            continue
        selected.append(component_cls)
    return tuple(selected)


class G16V3ComponentNode(BaseDataClassWithUnit):
    node_name: str = Field(default="g16.frame")
    component: G16V3BaseComponent | None = Field(default=None)
    children: list[G16V3ComponentNode] = Field(default_factory=list)
    include_in_aggregation: bool = Field(default=True, exclude=True, repr=False)
    include_in_render: bool = Field(default=True, exclude=True, repr=False)

    def iter_nodes(self) -> Iterator[G16V3ComponentNode]:
        yield self
        for child in self.children:
            yield from child.iter_nodes()

    def iter_components(self, *, include_synthetic: bool = False) -> Iterator[G16V3BaseComponent]:
        if self.component is not None and (self.include_in_aggregation or include_synthetic):
            yield self.component
        for child in self.children:
            yield from child.iter_components(include_synthetic=include_synthetic)

    def render_fakeg(self, **kwargs: Any) -> str:
        if self.component is not None:
            if not self.include_in_render:
                return ""
            return self.component.render_fakeg(**kwargs)
        return "\n".join(
            rendered
            for rendered in (child.render_fakeg(**kwargs).strip("\n") for child in self.children)
            if rendered
        )

    def render_subtree(self, **kwargs: Any) -> str:
        if self.component is not None:
            own_render = self.component.render_fakeg(**kwargs).strip("\n")
            child_renders = [
                rendered
                for rendered in (
                    child.render_subtree(**kwargs).strip("\n") for child in self.children
                )
                if rendered
            ]
            parts = [part for part in [own_render, *child_renders] if part]
            return "\n".join(parts)
        return "\n".join(
            rendered
            for rendered in (child.render_subtree(**kwargs).strip("\n") for child in self.children)
            if rendered
        )


class G16V3ComponentTree(BaseDataClassWithUnit):
    root: G16V3ComponentNode = Field(default_factory=G16V3ComponentNode)
    source_text: str = Field(default="", exclude=True, repr=False)
    only_extract_structure: bool = Field(default=False, exclude=True, repr=False)
    synthetic_children_expanded: bool = Field(default=False, exclude=True, repr=False)
    payloads_rawified: bool = Field(default=False, exclude=True, repr=False)

    def iter_nodes(self) -> Iterator[G16V3ComponentNode]:
        return self.root.iter_nodes()

    def iter_components(self, *, include_synthetic: bool = False) -> Iterator[G16V3BaseComponent]:
        return self.root.iter_components(include_synthetic=include_synthetic)

    def component_names(self, *, include_synthetic: bool = False) -> list[str]:
        return [
            component.component_name
            for component in self.iter_components(include_synthetic=include_synthetic)
        ]

    def find_nodes(self, node_name: str) -> list[G16V3ComponentNode]:
        return [node for node in self.iter_nodes() if node.node_name == node_name]

    def find_nodes_by_prefix(self, node_name_prefix: str) -> list[G16V3ComponentNode]:
        return [node for node in self.iter_nodes() if node.node_name.startswith(node_name_prefix)]

    def render_node(self, node_name: str, **kwargs: Any) -> str:
        nodes = self.find_nodes(node_name)
        return "\n".join(
            rendered
            for rendered in (node.render_subtree(**kwargs).strip("\n") for node in nodes)
            if rendered
        )

    def render_fakeg(self, **kwargs: Any) -> str:
        return self.root.render_fakeg(**kwargs)

    def validate_contracts(self) -> list[str]:
        issues: list[str] = []
        for node in self.iter_nodes():
            component = node.component
            if component is None:
                continue
            child_names = [
                child.component.component_name
                for child in node.children
                if child.component is not None
            ]
            allowed = set(component.allowed_child_component_names)
            required = set(component.required_child_component_names)
            repeatable = set(component.repeatable_child_component_names)
            if allowed:
                for child_name in child_names:
                    if child_name not in allowed:
                        issues.append(
                            f"{component.component_name}: child {child_name} not allowed by contract"
                        )
            for required_name in required:
                if required_name not in child_names:
                    issues.append(
                        f"{component.component_name}: required child {required_name} missing"
                    )
            seen_counts: dict[str, int] = {}
            for child_name in child_names:
                seen_counts[child_name] = seen_counts.get(child_name, 0) + 1
            for child_name, count in seen_counts.items():
                if count > 1 and child_name not in repeatable:
                    issues.append(
                        f"{component.component_name}: child {child_name} repeated {count} times but not repeatable"
                    )
        return issues


class G16V3ComponentTreeBuilder:
    @staticmethod
    def _append_synthetic_child(
        parent_node: G16V3ComponentNode,
        child_spec: G16SyntheticChildSpec,
    ) -> G16V3ComponentNode | None:
        normalized_payload = {
            key: value for key, value in child_spec.payload.items() if _has_meaningful_value(value)
        }
        if not normalized_payload:
            return None
        parent_component = parent_node.component
        child_component = child_spec.component_cls(
            span_start=parent_component.span_start if parent_component is not None else 0,
            span_end=parent_component.span_end if parent_component is not None else 0,
            raw_text="",
            payload=dict(normalized_payload),
        )
        child_node = G16V3ComponentNode(
            node_name=child_spec.node_name or child_component.component_name,
            component=child_component,
            include_in_aggregation=child_spec.include_in_aggregation,
            include_in_render=child_spec.include_in_render,
        )
        parent_node.children.append(child_node)
        return child_node

    @classmethod
    def expand_synthetic_children(cls, tree: G16V3ComponentTree) -> G16V3ComponentTree:
        if tree.synthetic_children_expanded:
            return tree
        for node in tree.root.children:
            component = node.component
            if component is None:
                continue
            for child_spec in component.build_synthetic_children():
                cls._append_synthetic_child(node, child_spec)
        tree.synthetic_children_expanded = True
        return tree

    @staticmethod
    def _rawify_tree_payloads(tree: G16V3ComponentTree) -> G16V3ComponentTree:
        for node in tree.iter_nodes():
            if node.component is not None:
                node.component.payload = _rawify_payload(node.component.payload)
        tree.payloads_rawified = True
        return tree

    @classmethod
    def from_g16_text(
        cls, full_text: str, *, only_extract_structure: bool = False
    ) -> G16V3ComponentTree:
        return build_g16log_v3_component_tree(
            full_text, only_extract_structure=only_extract_structure
        )

    @classmethod
    def from_frame_data(
        cls, frame: Any, *, expand_synthetic_children: bool = True
    ) -> G16V3ComponentTree:
        root = G16V3ComponentNode(node_name="g16.frame")
        for component_cls in G16LOG_V3_COMPONENT_REGISTRY:
            if not component_cls.can_build_from_frame(frame):
                continue
            payload = component_cls.build_payload_from_frame(frame)
            if not payload:
                continue
            component = component_cls(raw_text="", payload=payload)
            root.children.append(
                G16V3ComponentNode(node_name=component.component_name, component=component)
            )
        tree = G16V3ComponentTree(root=root, source_text="", only_extract_structure=False)
        if expand_synthetic_children:
            cls.expand_synthetic_children(tree)
        return tree


def build_g16log_v3_component_tree(
    full_text: str, *, only_extract_structure: bool = False
) -> G16V3ComponentTree:
    candidates: list[tuple[int, int, int, type[G16V3BaseComponent]]] = []
    context = G16ScanContext(full_text=full_text)
    component_classes = get_g16log_v3_candidate_component_classes(full_text, context=context)
    for order, component_cls in enumerate(component_classes):
        for span_start, span_end in component_cls.iter_spans(full_text, context=context):
            if span_end <= span_start:
                continue
            candidates.append((span_start, span_end, order, component_cls))
    candidates.sort(key=lambda item: (item[0], item[1], item[2]))

    filtered: list[tuple[int, int, int, type[G16V3BaseComponent]]] = []
    for candidate in candidates:
        span_start, span_end, _order, component_cls = candidate
        if filtered and span_start == filtered[-1][0] and span_end == filtered[-1][1]:
            continue
        if filtered and span_start < filtered[-1][1] and span_end <= filtered[-1][1]:
            previous_component_cls = filtered[-1][3]
            previous_policy = previous_component_cls.primary_boundary_spec().overlap_policy
            current_policy = component_cls.primary_boundary_spec().overlap_policy
            shared_end_allowed = (
                span_end == filtered[-1][1]
                and {previous_policy, current_policy} != {"forbid"}
                and (previous_policy == "allow_shared_end" or current_policy == "allow_shared_end")
            )
            if shared_end_allowed:
                filtered.append(candidate)
                continue
            continue
        filtered.append(candidate)

    root = G16V3ComponentNode(node_name="g16.frame")
    for span_start, span_end, _order, component_cls in filtered:
        component = component_cls(
            span_start=span_start,
            span_end=span_end,
            raw_text=full_text[span_start:span_end],
        )
        root.children.append(
            G16V3ComponentNode(node_name=component.component_name, component=component)
        )
    return G16V3ComponentTree(
        root=root,
        source_text=full_text,
        only_extract_structure=only_extract_structure,
    )


def aggregate_g16log_v3_tree(
    tree: G16V3ComponentTree, *, additional_data: Mapping[str, Any] | None = None
) -> dict[str, Any]:
    infos: dict[str, Any] = {"qm_software": "Gaussian"}
    if additional_data is not None:
        infos.update(dict(additional_data))

    payloads_rawified = tree.payloads_rawified
    for component in tree.iter_components():
        payload = (
            _restore_raw_payload(component.payload) if payloads_rawified else component.payload
        )
        component.aggregate_into(infos, payload)

    if "__running_time_accumulated_seconds__" in infos:
        infos["running_time"] = infos.pop("__running_time_accumulated_seconds__") * atom_ureg.second

    return infos
