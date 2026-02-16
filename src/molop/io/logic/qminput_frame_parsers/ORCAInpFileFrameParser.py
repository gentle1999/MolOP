"""
Author: TMJ
Date: 2026-02-10 00:00:00
LastEditors: TMJ
LastEditTime: 2026-02-10 00:00:00
Description: ORCA input frame parsers
"""

from __future__ import annotations

import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any, cast

import numpy as np
from pydantic import PrivateAttr
from rdkit import Chem

from molop.io.base_models.FrameParser import BaseFrameParser, _HasParseMethod
from molop.io.logic.qminput_frame_models.ORCAInpFileFrame import (
    ORCAInpFileFrameDisk,
    ORCAInpFileFrameMemory,
)
from molop.unit import atom_ureg


pt = Chem.GetPeriodicTable()

_CARTESIAN_TYPES = {"xyz", "cart", "cartesian"}
_NON_CARTESIAN_TYPES = {"int", "internal", "gzmt"}
_EXTERNAL_TYPES = {"xyzfile", "gzmtfile"}
_SUPPORTED_CTYPES = _CARTESIAN_TYPES | _NON_CARTESIAN_TYPES | _EXTERNAL_TYPES
_DUMMY_SYMBOLS = {"DA", "X", "XX"}
_ATOM_TOKEN_PATTERN = re.compile(r"^(?P<symbol>[A-Za-z]+)(?:\((?P<fragment>-?\d+)\))?$")


@dataclass(slots=True)
class _LineSpan:
    start: int
    end: int
    line: str


@dataclass(slots=True)
class _GeometrySection:
    preamble: str
    geometry: str
    postamble: str
    ctype: str
    charge: int | None
    multiplicity: int | None
    units: str | None
    coordinate_lines: list[str]


@dataclass(slots=True)
class _CartesianParseResult:
    atoms: list[int]
    coords: list[tuple[float, float, float]]
    entries: list[dict[str, Any]]
    fragments: list[list[int]]
    ghost_markers: list[int]
    dummy_markers: list[int]
    point_charges: list[dict[str, Any]]
    freeze_markers: list[int]
    isotope_tokens: dict[int, str]
    nuclear_charge_tokens: dict[int, str]


def _build_line_spans(block: str) -> list[_LineSpan]:
    spans: list[_LineSpan] = []
    cursor = 0
    for line in block.splitlines(keepends=True):
        next_cursor = cursor + len(line)
        spans.append(_LineSpan(start=cursor, end=next_cursor, line=line))
        cursor = next_cursor
    if not spans and block:
        spans.append(_LineSpan(start=0, end=len(block), line=block))
    return spans


def _parse_float(token: str) -> float | None:
    cleaned = token.strip().rstrip(",")
    if not cleaned:
        return None
    try:
        return float(cleaned.replace("D", "E").replace("d", "e"))
    except ValueError:
        return None


def _parse_int(token: str) -> int | None:
    cleaned = token.strip().strip('"').strip("'").rstrip(",")
    if not cleaned:
        return None
    try:
        return int(cleaned)
    except ValueError:
        as_float = _parse_float(cleaned)
        if as_float is None:
            return None
        if float(as_float).is_integer():
            return int(as_float)
        return None


def _strip_inline_comment(line: str) -> str:
    stripped = line.lstrip()
    if stripped.startswith("#"):
        return ""
    if "#" in line:
        return line.split("#", 1)[0]
    return line


def _split_key_value(line: str) -> tuple[str, str]:
    content = _strip_inline_comment(line).strip()
    if "=" in content:
        key, value = content.split("=", 1)
        return key.strip(), value.strip()
    tokens = content.split(maxsplit=1)
    if not tokens:
        return "", ""
    if len(tokens) == 1:
        return tokens[0], ""
    return tokens[0], tokens[1]


def _parse_symbol_token(token: str) -> tuple[str | None, int | None, bool, bool]:
    cleaned = token.strip()
    is_ghost = cleaned.endswith(":")
    if is_ghost:
        cleaned = cleaned[:-1]

    if cleaned.isdigit():
        atomic_number = int(cleaned)
        if atomic_number > 0:
            return pt.GetElementSymbol(atomic_number), None, is_ghost, False
        return None, None, is_ghost, False

    matched = _ATOM_TOKEN_PATTERN.fullmatch(cleaned)
    if matched is None:
        return None, None, is_ghost, False

    raw_symbol = matched.group("symbol")
    fragment_text = matched.group("fragment")
    symbol_upper = raw_symbol.upper()
    is_dummy = symbol_upper in _DUMMY_SYMBOLS
    if symbol_upper == "DA":
        symbol = "DA"
    elif symbol_upper == "XX":
        symbol = "Xx"
    elif len(raw_symbol) == 1:
        symbol = raw_symbol.upper()
    else:
        symbol = raw_symbol[0].upper() + raw_symbol[1:].lower()

    fragment_id = int(fragment_text) if fragment_text is not None else None
    return symbol, fragment_id, is_ghost, is_dummy


def _extract_coords_from_tokens(
    tokens: Sequence[str],
) -> tuple[tuple[float, float, float] | None, bool, list[str]]:
    coords: list[float] = []
    freeze_marker = False
    idx = 0
    while idx < len(tokens) and len(coords) < 3:
        token = tokens[idx]
        if token == "$":
            freeze_marker = True
            idx += 1
            continue
        candidate = token
        if candidate.startswith("$"):
            freeze_marker = True
            candidate = candidate[1:]
        value = _parse_float(candidate)
        if value is None:
            break
        coords.append(value)
        idx += 1

    if len(coords) != 3:
        return None, freeze_marker, list(tokens[idx:])

    tails = list(tokens[idx:])
    if any("$" in tail for tail in tails):
        freeze_marker = True
    return (coords[0], coords[1], coords[2]), freeze_marker, tails


def _is_coords_open_line(line: str) -> bool:
    lower = line.lower()
    return lower.startswith("coords") and (len(lower) == 6 or lower[6].isspace())


def _extract_star_geometry(block: str) -> _GeometrySection | None:
    spans = _build_line_spans(block)
    if not spans:
        return None

    for idx, span in enumerate(spans):
        stripped = span.line.strip()
        if not stripped.startswith("*"):
            continue

        header = stripped[1:].strip()
        if not header:
            continue

        header_tokens = header.split()
        ctype = header_tokens[0].lower()
        if ctype not in _SUPPORTED_CTYPES:
            continue

        charge = _parse_int(header_tokens[1]) if len(header_tokens) > 1 else None
        multiplicity = _parse_int(header_tokens[2]) if len(header_tokens) > 2 else None
        units = header_tokens[3] if len(header_tokens) > 3 else None

        has_closing_star = False
        if ctype in _EXTERNAL_TYPES:
            end_idx = idx + 1
        else:
            end_idx = len(spans)
            for close_idx in range(idx + 1, len(spans)):
                if spans[close_idx].line.strip() == "*":
                    has_closing_star = True
                    end_idx = close_idx + 1
                    break

        start_pos = spans[idx].start
        end_pos = spans[end_idx - 1].end if end_idx > 0 else len(block)

        if ctype in _CARTESIAN_TYPES:
            body_end = end_idx - 1 if has_closing_star else end_idx
            coordinate_lines = [
                spans[line_idx].line.rstrip("\r\n") for line_idx in range(idx + 1, body_end)
            ]
        else:
            coordinate_lines = []

        return _GeometrySection(
            preamble=block[:start_pos],
            geometry=block[start_pos:end_pos],
            postamble=block[end_pos:],
            ctype=ctype,
            charge=charge,
            multiplicity=multiplicity,
            units=units,
            coordinate_lines=coordinate_lines,
        )

    return None


def _extract_percent_coords_geometry(block: str) -> _GeometrySection | None:
    spans = _build_line_spans(block)
    if not spans:
        return None

    for idx, span in enumerate(spans):
        stripped = span.line.strip().lower()
        if not stripped.startswith("%coords"):
            continue

        nested_coords_depth = 0
        end_idx = len(spans)
        for line_idx in range(idx + 1, len(spans)):
            candidate = spans[line_idx].line.strip().lower()
            if _is_coords_open_line(candidate):
                nested_coords_depth += 1
                continue
            if candidate == "end":
                if nested_coords_depth > 0:
                    nested_coords_depth -= 1
                    continue
                end_idx = line_idx + 1
                break

        ctype: str | None = None
        charge: int | None = None
        multiplicity: int | None = None
        units: str | None = None
        coordinate_lines: list[str] = []
        in_coords = False

        for line_idx in range(idx + 1, end_idx):
            raw_line = spans[line_idx].line.rstrip("\r\n")
            stripped_line = raw_line.strip()
            if not stripped_line:
                continue
            if stripped_line.startswith("#"):
                continue

            lowered = stripped_line.lower()
            if _is_coords_open_line(lowered):
                in_coords = True
                continue
            if lowered == "end":
                if in_coords:
                    in_coords = False
                    continue
                break

            if in_coords:
                coordinate_lines.append(raw_line)
                continue

            key, value = _split_key_value(stripped_line)
            key_lower = key.lower()
            value_clean = value.strip().strip('"').strip("'")
            if key_lower == "ctyp":
                ctype = value_clean.lower()
            elif key_lower == "charge":
                charge = _parse_int(value_clean)
            elif key_lower == "mult":
                multiplicity = _parse_int(value_clean)
            elif key_lower == "units":
                units = value_clean

        start_pos = spans[idx].start
        end_pos = spans[end_idx - 1].end if end_idx > 0 else len(block)
        return _GeometrySection(
            preamble=block[:start_pos],
            geometry=block[start_pos:end_pos],
            postamble=block[end_pos:],
            ctype=(ctype or "xyz").lower(),
            charge=charge,
            multiplicity=multiplicity,
            units=units,
            coordinate_lines=coordinate_lines,
        )

    return None


def _is_bohr_unit(unit_hint: str | None) -> bool:
    if unit_hint is None:
        return False
    normalized = unit_hint.strip().strip('"').strip("'").lower()
    return normalized in {
        "bohr",
        "bohrs",
        "a0",
        "au",
        "atomic",
        "atomicunit",
        "atomicunits",
    }


def _parse_cartesian_coordinate_lines(lines: Sequence[str]) -> _CartesianParseResult:
    atoms: list[int] = []
    coords: list[tuple[float, float, float]] = []
    entries: list[dict[str, Any]] = []
    fragments: list[list[int]] = []
    ghost_markers: list[int] = []
    dummy_markers: list[int] = []
    point_charges: list[dict[str, Any]] = []
    freeze_markers: list[int] = []
    isotope_tokens: dict[int, str] = {}
    nuclear_charge_tokens: dict[int, str] = {}

    line_index = 0
    for raw_line in lines:
        parse_line = _strip_inline_comment(raw_line).strip()
        if not parse_line:
            continue
        tokens = parse_line.split()
        if not tokens:
            continue

        line_index += 1
        if tokens[0].upper() == "Q":
            if len(tokens) >= 5:
                point_charge = _parse_float(tokens[1])
                point_x = _parse_float(tokens[2])
                point_y = _parse_float(tokens[3])
                point_z = _parse_float(tokens[4])
                if (
                    point_charge is not None
                    and point_x is not None
                    and point_y is not None
                    and point_z is not None
                ):
                    point_entry: dict[str, Any] = {
                        "charge": point_charge,
                        "x": point_x,
                        "y": point_y,
                        "z": point_z,
                    }
                    entries.append(
                        {
                            "kind": "q",
                            "charge": point_charge,
                            "x": point_x,
                            "y": point_y,
                            "z": point_z,
                        }
                    )
                    point_charges.append(point_entry)
            continue

        symbol, fragment_id, is_ghost, is_dummy = _parse_symbol_token(tokens[0])
        if symbol is None:
            continue

        coord_triplet, freeze_marker, tails = _extract_coords_from_tokens(tokens[1:])
        if coord_triplet is None:
            continue

        if is_dummy:
            dummy_markers.append(line_index)
            continue

        atomic_number = pt.GetAtomicNumber(symbol)
        if atomic_number <= 0:
            continue

        atom_index = len(atoms) + 1
        atoms.append(atomic_number)
        coords.append(coord_triplet)
        entries.append(
            {
                "kind": "atom",
                "symbol": symbol,
                "x": coord_triplet[0],
                "y": coord_triplet[1],
                "z": coord_triplet[2],
            }
        )

        if fragment_id is not None:
            fragments.append([atom_index, fragment_id])
        if is_ghost:
            ghost_markers.append(atom_index)
        if freeze_marker:
            freeze_markers.append(atom_index)

        for tail in tails:
            tail_upper = tail.upper()
            if tail_upper.startswith("M="):
                isotope_tokens[atom_index] = tail.split("=", 1)[1]
            if tail_upper.startswith("Z="):
                nuclear_charge_tokens[atom_index] = tail.split("=", 1)[1]

    return _CartesianParseResult(
        atoms=atoms,
        coords=coords,
        entries=entries,
        fragments=fragments,
        ghost_markers=ghost_markers,
        dummy_markers=dummy_markers,
        point_charges=point_charges,
        freeze_markers=freeze_markers,
        isotope_tokens=isotope_tokens,
        nuclear_charge_tokens=nuclear_charge_tokens,
    )


class ORCAInpFileFrameParserMixin:
    _last_coords_ctype: str | None = PrivateAttr(default=None)
    _last_external_path: str | None = PrivateAttr(default=None)
    _last_point_charges: list[dict[str, float]] | None = PrivateAttr(default=None)

    def _parse_frame(self) -> Mapping[str, Any]:
        typed_self = cast(_HasParseMethod, self)
        block = typed_self._block

        self._last_coords_ctype = None
        self._last_external_path = None
        self._last_point_charges = None
        section = _extract_star_geometry(block)
        if section is None:
            section = _extract_percent_coords_geometry(block)

        if section is None:
            return {}

        self._last_coords_ctype = section.ctype
        if section.ctype in _EXTERNAL_TYPES:
            self._last_external_path = section.units

        metadata: dict[str, Any] = {
            "orca_raw_preamble": section.preamble,
            "orca_raw_postamble": section.postamble,
        }

        # Best-effort: ORCA keywords are typically provided via one or more leading `!` lines.
        # Preserve them as normalized text (without the leading '!') for the QM-input base fields.
        keyword_lines: list[str] = []
        for line in section.preamble.splitlines():
            stripped = line.lstrip()
            if stripped.startswith("!"):
                keyword_lines.append(stripped[1:].strip())
        if keyword_lines:
            metadata["keywords"] = "\n".join(keyword_lines)

        if section.charge is not None:
            metadata["charge"] = section.charge
        if section.multiplicity is not None:
            metadata["multiplicity"] = section.multiplicity

        if section.ctype in _CARTESIAN_TYPES and section.coordinate_lines:
            parsed = _parse_cartesian_coordinate_lines(section.coordinate_lines)
            use_bohr_units = _is_bohr_unit(section.units)
            if parsed.atoms:
                coords_array = np.array(parsed.coords, dtype=np.float32)
                if use_bohr_units:
                    bohr_coords = coords_array * atom_ureg.bohr
                    coords_quantity = cast(Any, bohr_coords).to(atom_ureg.angstrom)
                else:
                    coords_quantity = coords_array * atom_ureg.angstrom
                metadata["atoms"] = parsed.atoms
                metadata["coords"] = coords_quantity
            if parsed.point_charges:
                point_charges: list[dict[str, float]] = []
                for point_charge in parsed.point_charges:
                    charge = point_charge.get("charge")
                    x = point_charge.get("x")
                    y = point_charge.get("y")
                    z = point_charge.get("z")
                    if not all(isinstance(value, (int, float)) for value in (charge, x, y, z)):
                        continue

                    charge_f = float(cast(float, charge))
                    x_f = float(cast(float, x))
                    y_f = float(cast(float, y))
                    z_f = float(cast(float, z))

                    parsed_point_charge = {
                        "charge": charge_f,
                        "x": x_f,
                        "y": y_f,
                        "z": z_f,
                    }
                    if use_bohr_units:
                        parsed_point_charge["x"] = float(
                            cast(Any, parsed_point_charge["x"] * atom_ureg.bohr)
                            .to(atom_ureg.angstrom)
                            .m
                        )
                        parsed_point_charge["y"] = float(
                            cast(Any, parsed_point_charge["y"] * atom_ureg.bohr)
                            .to(atom_ureg.angstrom)
                            .m
                        )
                        parsed_point_charge["z"] = float(
                            cast(Any, parsed_point_charge["z"] * atom_ureg.bohr)
                            .to(atom_ureg.angstrom)
                            .m
                        )
                    point_charges.append(parsed_point_charge)

                if point_charges:
                    self._last_point_charges = point_charges

        return metadata

    def parse(self, block: str, *, additional_data: dict[str, Any] | None = None) -> Any:
        frame = cast(Any, super()).parse(block, additional_data=additional_data)
        if hasattr(frame, "_orca_coords_ctype"):
            frame._orca_coords_ctype = self._last_coords_ctype
        if hasattr(frame, "_orca_external_path"):
            frame._orca_external_path = self._last_external_path
        if hasattr(frame, "_orca_point_charges"):
            frame._orca_point_charges = self._last_point_charges
        self._last_coords_ctype = None
        self._last_external_path = None
        self._last_point_charges = None
        return frame


class ORCAInpFileFrameParserMemory(
    ORCAInpFileFrameParserMixin, BaseFrameParser[ORCAInpFileFrameMemory]
):
    _file_frame_class_ = ORCAInpFileFrameMemory


class ORCAInpFileFrameParserDisk(
    ORCAInpFileFrameParserMixin, BaseFrameParser[ORCAInpFileFrameDisk]
):
    _file_frame_class_ = ORCAInpFileFrameDisk
