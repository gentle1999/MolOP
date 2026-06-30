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
from typing import Any, Literal, cast

from rdkit import Chem

from molop.io.base_models.DataClasses import (
    AtomInInternalCoords,
    CoordinateParameter,
    CoordinateParameters,
    InternalCoords,
)
from molop.io.base_models.FrameParser import BaseFrameParser, _HasParseMethod
from molop.io.logic.qminput_frame_models.ORCAInpFileFrame import (
    ORCAAtomBasisOverride,
    ORCABlock,
    ORCABlockLine,
    ORCACommentLine,
    ORCAGeometry,
    ORCAGeometryAtom,
    ORCAInpFileFrameDisk,
    ORCAInpFileFrameMemory,
    ORCAKeywordLine,
)
from molop.unit import atom_ureg


pt = Chem.GetPeriodicTable()

_CARTESIAN_TYPES = {"xyz", "cart", "cartesian"}
_NON_CARTESIAN_TYPES = {"int", "internal", "gzmt"}
_EXTERNAL_TYPES = {"xyzfile", "gzmtfile", "pdbfile"}
_SUPPORTED_CTYPES = _CARTESIAN_TYPES | _NON_CARTESIAN_TYPES | _EXTERNAL_TYPES
_DUMMY_SYMBOLS = {"DA", "X", "XX"}
_ATOM_TOKEN_PATTERN = re.compile(r"^(?P<symbol>[A-Za-z]+)(?:\((?P<fragment>-?\d+)\))?$")
_GEOMETRY_CLOSE_TOKENS = {"*", "end", "edn"}
_NESTED_BLOCK_KEYWORDS_BY_BLOCK = {
    "compound": {"new_step"},
    "geom": {
        "connectfragments",
        "constraints",
        "constrainfragments",
        "hybrid_hess",
        "inhess",
        "hess_internal",
        "modify_internal",
        "potentials",
        "scan",
        "ts_active_atoms",
        "ts_mode",
    },
    "mrci": {"newblock", "refs"},
}


@dataclass(slots=True)
class _LineSpan:
    start: int
    end: int
    line: str


@dataclass(slots=True)
class _GeometrySection:
    ctype: Literal[
        "xyz", "cart", "cartesian", "int", "internal", "gzmt", "xyzfile", "gzmtfile", "pdbfile"
    ]
    charge: int
    multiplicity: int
    units: str | None
    external_path: str | None
    coordinate_lines: list[str]
    line_start: int
    line_end: int
    source: Literal["star", "percent_coords"]


@dataclass(slots=True)
class _BlockSpan:
    block: ORCABlock
    line_start: int
    line_end: int


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
    if "," in cleaned:
        cleaned = cleaned.split(",", 1)[0]
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
        if as_float is None or not float(as_float).is_integer():
            return None
        return int(as_float)


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


def _strip_balanced_braces(token: str) -> str:
    cleaned = token.strip()
    if cleaned.startswith("{") and cleaned.endswith("}"):
        return cleaned[1:-1].strip()
    return cleaned


def _parse_parameter_blocks(blocks: Sequence[ORCABlock]) -> CoordinateParameters:
    parameters: list[CoordinateParameter] = []
    for block in blocks:
        if block.name.lower() != "paras":
            continue
        for line in block.lines:
            key, value = _split_key_value(line.text)
            if not key or not value:
                continue
            parts = [part.strip() for part in value.split(",")]
            start = _parse_float(parts[0]) if parts else None
            stop = _parse_float(parts[1]) if len(parts) > 1 else None
            steps = _parse_int(parts[2]) if len(parts) > 2 else None
            parameters.append(
                CoordinateParameter(
                    name=key,
                    raw_value=value,
                    start=start,
                    stop=stop,
                    steps=steps,
                )
            )
    return CoordinateParameters(items=parameters)


def _resolve_float_expression(token: str, parameters: Mapping[str, float]) -> float | None:
    expression = _strip_balanced_braces(token)
    parsed = _parse_float(expression)
    if parsed is not None:
        return parsed
    if expression in parameters:
        return parameters[expression]

    matched = re.fullmatch(
        r"(?P<name>[A-Za-z_][A-Za-z0-9_]*)\s*(?P<op>[+-])\s*(?P<offset>[+-]?\d+(?:\.\d*)?)",
        expression,
    )
    if matched is None:
        return None
    base = parameters.get(matched.group("name"))
    offset = _parse_float(matched.group("offset"))
    if base is None or offset is None:
        return None
    return base + offset if matched.group("op") == "+" else base - offset


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


def _clean_quoted_token(token: str) -> str:
    return token.strip().strip('"').strip("'")


def _parse_atom_basis_overrides(tokens: Sequence[str]) -> list[ORCAAtomBasisOverride]:
    overrides: list[ORCAAtomBasisOverride] = []
    idx = 0
    while idx < len(tokens):
        token_lower = tokens[idx].lower()
        if token_lower not in {"newgto", "newauxgto"}:
            idx += 1
            continue

        directive_tokens: list[str] = []
        scan_idx = idx + 1
        while scan_idx < len(tokens):
            directive_tokens.append(tokens[scan_idx])
            if tokens[scan_idx].lower() == "end":
                scan_idx += 1
                break
            scan_idx += 1

        basis_tokens = [token for token in directive_tokens if token.lower() != "end"]
        basis_set = _clean_quoted_token(basis_tokens[0]) if len(basis_tokens) == 1 else None
        overrides.append(
            ORCAAtomBasisOverride(
                kind=cast(Literal["newgto", "newauxgto"], token_lower),
                basis_set=basis_set,
                tokens=[_clean_quoted_token(token) for token in directive_tokens],
            )
        )
        idx = scan_idx
    return overrides


def _is_coords_open_line(line: str) -> bool:
    lower = line.lower()
    return lower.startswith("coords") and (len(lower) == 6 or lower[6].isspace())


def _normalize_ctype(
    ctype: str | None,
) -> Literal[
    "xyz", "cart", "cartesian", "int", "internal", "gzmt", "xyzfile", "gzmtfile", "pdbfile"
]:
    normalized = (ctype or "xyz").strip().strip('"').strip("'").lower()
    if normalized in _SUPPORTED_CTYPES:
        return cast(
            Literal[
                "xyz",
                "cart",
                "cartesian",
                "int",
                "internal",
                "gzmt",
                "xyzfile",
                "gzmtfile",
                "pdbfile",
            ],
            normalized,
        )
    return "xyz"


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
        ctype = _normalize_ctype(header_tokens[0] if header_tokens else None)
        if ctype not in _SUPPORTED_CTYPES:
            continue

        charge = _parse_int(header_tokens[1]) if len(header_tokens) > 1 else None
        multiplicity = _parse_int(header_tokens[2]) if len(header_tokens) > 2 else None
        external_path = (
            " ".join(header_tokens[3:]).strip().strip('"').strip("'")
            if ctype in _EXTERNAL_TYPES and len(header_tokens) > 3
            else None
        )
        units = None if ctype in _EXTERNAL_TYPES else header_tokens[3] if len(header_tokens) > 3 else None

        has_closing_star = False
        if ctype in _EXTERNAL_TYPES:
            end_idx = idx + 1
        else:
            end_idx = len(spans)
            for close_idx in range(idx + 1, len(spans)):
                if spans[close_idx].line.strip().lower() in _GEOMETRY_CLOSE_TOKENS:
                    has_closing_star = True
                    end_idx = close_idx + 1
                    break

        coordinate_lines: list[str] = []
        if ctype not in _EXTERNAL_TYPES:
            body_end = end_idx - 1 if has_closing_star else end_idx
            coordinate_lines = [
                spans[line_idx].line.rstrip("\r\n") for line_idx in range(idx + 1, body_end)
            ]

        return _GeometrySection(
            ctype=ctype,
            charge=charge if charge is not None else 0,
            multiplicity=multiplicity if multiplicity is not None else 1,
            units=units,
            external_path=external_path,
            coordinate_lines=coordinate_lines,
            line_start=idx,
            line_end=end_idx,
            source="star",
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
            if not stripped_line or stripped_line.startswith("#"):
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
                ctype = value_clean
            elif key_lower == "charge":
                charge = _parse_int(value_clean)
            elif key_lower == "mult":
                multiplicity = _parse_int(value_clean)
            elif key_lower == "units":
                units = value_clean

        return _GeometrySection(
            ctype=_normalize_ctype(ctype),
            charge=charge if charge is not None else 0,
            multiplicity=multiplicity if multiplicity is not None else 1,
            units=units,
            external_path=None,
            coordinate_lines=coordinate_lines,
            line_start=idx,
            line_end=end_idx,
            source="percent_coords",
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


def _convert_length_to_angstrom(value: float, unit_hint: str | None) -> float:
    if _is_bohr_unit(unit_hint):
        return value * 0.529177210903
    return value


def _parse_cartesian_coordinate_lines(
    lines: Sequence[str], parameters: Mapping[str, float], unit_hint: str | None
) -> tuple[list[ORCAGeometryAtom], list[dict[str, float]]]:
    atoms: list[ORCAGeometryAtom] = []
    point_charges: list[dict[str, float]] = []
    atom_index = 0

    for raw_line in lines:
        parse_line = _strip_inline_comment(raw_line).strip()
        if not parse_line:
            continue
        tokens = parse_line.split()
        if not tokens:
            continue

        if tokens[0].upper() == "Q":
            if len(tokens) < 5:
                continue
            charge = _resolve_float_expression(tokens[1], parameters)
            x = _resolve_float_expression(tokens[2], parameters)
            y = _resolve_float_expression(tokens[3], parameters)
            z = _resolve_float_expression(tokens[4], parameters)
            if charge is None or x is None or y is None or z is None:
                continue
            point_charges.append(
                {
                    "charge": float(charge),
                    "x": _convert_length_to_angstrom(float(x), unit_hint),
                    "y": _convert_length_to_angstrom(float(y), unit_hint),
                    "z": _convert_length_to_angstrom(float(z), unit_hint),
                }
            )
            continue

        symbol, fragment_id, is_ghost, is_dummy = _parse_symbol_token(tokens[0])
        if symbol is None:
            continue
        resolved_coord_tokens: list[str] = []
        for token in tokens[1:]:
            resolved = _resolve_float_expression(token, parameters)
            if resolved is None:
                resolved_coord_tokens.append(token)
            else:
                resolved_coord_tokens.append(str(resolved))
        coord_triplet, freeze_marker, tails = _extract_coords_from_tokens(resolved_coord_tokens)
        if coord_triplet is None:
            continue

        atom_index += 1
        atomic_number = 0 if is_dummy else pt.GetAtomicNumber(symbol)
        if not is_dummy and atomic_number <= 0:
            continue

        isotope: str | None = None
        nuclear_charge: str | None = None
        for tail in tails:
            tail_upper = tail.upper()
            if tail_upper.startswith("M="):
                isotope = tail.split("=", 1)[1]
            elif tail_upper.startswith("Z="):
                nuclear_charge = tail.split("=", 1)[1]
        basis_overrides = _parse_atom_basis_overrides(tails)
        atom_basis_set = next(
            (
                override.basis_set
                for override in basis_overrides
                if override.kind == "newgto" and override.basis_set is not None
            ),
            None,
        )
        atom_auxiliary_basis_set = next(
            (
                override.basis_set
                for override in basis_overrides
                if override.kind == "newauxgto" and override.basis_set is not None
            ),
            None,
        )

        atoms.append(
            ORCAGeometryAtom(
                symbol=symbol,
                atomic_number=atomic_number or None,
                x=_convert_length_to_angstrom(coord_triplet[0], unit_hint),
                y=_convert_length_to_angstrom(coord_triplet[1], unit_hint),
                z=_convert_length_to_angstrom(coord_triplet[2], unit_hint),
                is_dummy=is_dummy,
                is_ghost=is_ghost,
                fragment_id=fragment_id,
                frozen=freeze_marker,
                isotope=isotope,
                nuclear_charge=nuclear_charge,
                basis_set=atom_basis_set,
                auxiliary_basis_set=atom_auxiliary_basis_set,
                basis_overrides=basis_overrides,
            )
        )

    _ = atom_index
    return atoms, point_charges


def _parse_internal_coordinate_lines(
    lines: Sequence[str], parameters: Mapping[str, float], unit_hint: str | None
) -> tuple[list[ORCAGeometryAtom], InternalCoords | None]:
    atoms: list[ORCAGeometryAtom] = []
    internal_atoms: list[AtomInInternalCoords] = []
    all_numeric = True

    for raw_line in lines:
        parse_line = _strip_inline_comment(raw_line).strip()
        if not parse_line:
            continue
        tokens = parse_line.split()
        if not tokens:
            continue
        symbol, fragment_id, is_ghost, is_dummy = _parse_symbol_token(tokens[0])
        if symbol is None:
            continue
        atomic_number = 0 if is_dummy else pt.GetAtomicNumber(symbol)
        if not is_dummy and atomic_number <= 0:
            continue

        def _ref(position: int, row_tokens: Sequence[str] = tokens) -> int:
            parsed = _parse_int(row_tokens[position]) if len(row_tokens) > position else None
            if parsed is None or parsed <= 0:
                return 0
            return parsed - 1

        distance = (
            _resolve_float_expression(tokens[4], parameters) if len(tokens) > 4 else 0.0
        )
        angle = _resolve_float_expression(tokens[5], parameters) if len(tokens) > 5 else 0.0
        dihedral = (
            _resolve_float_expression(tokens[6], parameters) if len(tokens) > 6 else 0.0
        )
        if distance is None or angle is None or dihedral is None:
            all_numeric = False
            atom_internal = None
        else:
            atom_internal = AtomInInternalCoords(
                symbol=symbol,
                distance_to_index=_ref(1),
                angle_to_index=_ref(2),
                dihedral_to_index=_ref(3),
                distance=_convert_length_to_angstrom(distance, unit_hint) * atom_ureg.angstrom,
                angle=angle * atom_ureg.degree,
                dihedral=dihedral * atom_ureg.degree,
                is_dummy=is_dummy,
                is_ghost=is_ghost,
            )
            internal_atoms.append(atom_internal)

        atoms.append(
            ORCAGeometryAtom(
                symbol=symbol,
                atomic_number=atomic_number or None,
                is_dummy=is_dummy,
                is_ghost=is_ghost,
                fragment_id=fragment_id,
                internal_coord=atom_internal,
            )
        )

    if not all_numeric or len(internal_atoms) != len(atoms):
        return atoms, None
    return atoms, InternalCoords(items=internal_atoms)


def _geometry_from_section(
    section: _GeometrySection | None, blocks: Sequence[ORCABlock]
) -> ORCAGeometry | None:
    if section is None:
        return None
    atoms: list[ORCAGeometryAtom] = []
    point_charges: list[dict[str, float]] = []
    internal_coords: InternalCoords | None = None
    coordinate_parameters = _parse_parameter_blocks(blocks)
    parameter_values = coordinate_parameters.as_value_map()
    if section.ctype in _CARTESIAN_TYPES and section.coordinate_lines:
        atoms, point_charges = _parse_cartesian_coordinate_lines(
            section.coordinate_lines, parameter_values, section.units
        )
    elif section.ctype in _NON_CARTESIAN_TYPES and section.coordinate_lines:
        atoms, internal_coords = _parse_internal_coordinate_lines(
            section.coordinate_lines, parameter_values, section.units
        )
    return ORCAGeometry(
        ctype=section.ctype,
        charge=section.charge,
        multiplicity=section.multiplicity,
        units=section.units,
        external_path=section.external_path,
        items=atoms,
        internal_coords=internal_coords,
        coordinate_parameters=coordinate_parameters if len(coordinate_parameters) > 0 else None,
        point_charges=point_charges,
        source=section.source,
    )


def _parse_percent_header(line: str) -> tuple[str, str] | None:
    stripped = line.strip()
    if not stripped.startswith("%"):
        return None
    rest = stripped[1:].strip()
    if not rest:
        return None
    tokens = rest.split(maxsplit=1)
    return tokens[0].lower(), tokens[1].strip() if len(tokens) > 1 else ""


def _is_top_level_directive(line: str) -> bool:
    stripped = line.strip()
    if not stripped:
        return False
    lowered = stripped.lower()
    return stripped.startswith(("!", "*", "%")) or lowered.startswith("$new_job")


def _first_code_token(line: str) -> str:
    content = _strip_inline_comment(line).strip()
    if not content:
        return ""
    return content.split(maxsplit=1)[0].lower()


def _code_tokens(line: str) -> list[str]:
    return _strip_inline_comment(line).strip().lower().split()


def _opens_nested_block(parent_name: str, line: str) -> bool:
    tokens = _code_tokens(line)
    if not tokens:
        return False
    return (
        tokens[0] in _NESTED_BLOCK_KEYWORDS_BY_BLOCK.get(parent_name.lower(), set())
        and "end" not in tokens[1:]
    )


def _closes_nested_block(parent_name: str, line: str) -> bool:
    return parent_name.lower() == "compound" and _first_code_token(line) == "step_end"


def _block_body_from_inline(inline: str) -> tuple[list[ORCABlockLine], bool]:
    if not inline:
        return [], False
    tokens = inline.split()
    if tokens and tokens[-1].lower() == "end":
        body = " ".join(tokens[:-1]).strip()
        return ([ORCABlockLine(text=body)] if body else []), True
    return [ORCABlockLine(text=inline)], False


def _extract_block_spans(block: str) -> list[_BlockSpan]:
    spans = _build_line_spans(block)
    result: list[_BlockSpan] = []
    idx = 0
    while idx < len(spans):
        header = _parse_percent_header(spans[idx].line)
        if header is None:
            idx += 1
            continue

        name, inline = header
        body_lines, inline_closed = _block_body_from_inline(inline)
        end_idx = idx + 1

        nested_depth = 1 if inline and not inline_closed and _opens_nested_block(name, inline) else 0
        if not inline_closed:
            scan_idx = idx + 1
            while scan_idx < len(spans):
                candidate = spans[scan_idx].line.strip()
                lowered = candidate.lower()
                if _is_coords_open_line(lowered):
                    nested_depth += 1
                    body_lines.append(ORCABlockLine(text=spans[scan_idx].line.rstrip("\r\n")))
                    scan_idx += 1
                    continue
                if _opens_nested_block(name, candidate):
                    nested_depth += 1
                    body_lines.append(ORCABlockLine(text=spans[scan_idx].line.rstrip("\r\n")))
                    scan_idx += 1
                    continue
                if _closes_nested_block(name, candidate) and nested_depth > 0:
                    nested_depth -= 1
                    body_lines.append(ORCABlockLine(text=spans[scan_idx].line.rstrip("\r\n")))
                    scan_idx += 1
                    continue
                if lowered == "end":
                    if nested_depth > 0:
                        nested_depth -= 1
                        body_lines.append(ORCABlockLine(text=spans[scan_idx].line.rstrip("\r\n")))
                        scan_idx += 1
                        continue
                    end_idx = scan_idx + 1
                    break
                if nested_depth == 0 and _is_top_level_directive(candidate):
                    end_idx = scan_idx
                    break
                body_lines.append(ORCABlockLine(text=spans[scan_idx].line.rstrip("\r\n")))
                scan_idx += 1
            else:
                end_idx = len(spans)

        result.append(
            _BlockSpan(
                block=ORCABlock(
                    name=name,
                    lines=[line for line in body_lines if line.text.strip()],
                    raw_header=spans[idx].line.rstrip("\r\n"),
                    raw_text="".join(span.line for span in spans[idx:end_idx]).rstrip("\r\n"),
                ),
                line_start=idx,
                line_end=end_idx,
            )
        )
        idx = max(end_idx, idx + 1)
    return result


def _line_in_spans(line_idx: int, spans: Sequence[tuple[int, int]]) -> bool:
    return any(start <= line_idx < end for start, end in spans)


class ORCAInpFileFrameParserMixin:
    def _parse_frame(self) -> Mapping[str, Any]:
        typed_self = cast(_HasParseMethod, self)
        block = typed_self._block
        lines = block.splitlines()

        section = _extract_star_geometry(block)
        if section is None:
            section = _extract_percent_coords_geometry(block)

        block_spans = _extract_block_spans(block)
        occupied_spans: list[tuple[int, int]] = [
            (span.line_start, span.line_end) for span in block_spans
        ]
        if section is not None:
            occupied_spans.append((section.line_start, section.line_end))

        comments: list[ORCACommentLine] = []
        keywords: list[ORCAKeywordLine] = []
        trailing_lines: list[str] = []

        for line_idx, line in enumerate(lines):
            if _line_in_spans(line_idx, occupied_spans):
                continue
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("#"):
                comments.append(ORCACommentLine(text=stripped[1:].strip()))
            elif stripped.startswith("!"):
                keywords.append(ORCAKeywordLine(text=stripped[1:].strip()))
            elif not stripped.lower().startswith("$new_job"):
                trailing_lines.append(line)

        geometry = _geometry_from_section(section, [span.block for span in block_spans])

        return {
            "comment_lines": comments,
            "keyword_lines": keywords,
            "blocks": [span.block for span in block_spans],
            "geometry": geometry,
            "trailing_lines": trailing_lines,
        }


class ORCAInpFileFrameParserMemory(
    ORCAInpFileFrameParserMixin, BaseFrameParser[ORCAInpFileFrameMemory]
):
    _file_frame_class_ = ORCAInpFileFrameMemory


class ORCAInpFileFrameParserDisk(
    ORCAInpFileFrameParserMixin, BaseFrameParser[ORCAInpFileFrameDisk]
):
    _file_frame_class_ = ORCAInpFileFrameDisk
