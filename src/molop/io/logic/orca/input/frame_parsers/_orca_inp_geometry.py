from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Literal, cast

from rdkit import Chem

from molop.io.base_models.DataClasses import (
    AtomInInternalCoords,
    CoordinateParameter,
    CoordinateParameters,
    InternalCoords,
)
from molop.io.logic.orca.common import (
    ORCAAtomBasisOverride,
    ORCABlock,
    ORCAGeometry,
    ORCAGeometryAtom,
)
from molop.io.logic.orca.input.frame_parsers._orca_inp_blocks import (
    build_line_spans,
    is_coords_open_line,
)
from molop.io.logic.orca.input.frame_parsers._orca_inp_patterns import orca_inp_patterns
from molop.io.logic.orca.input.frame_parsers._orca_inp_tokens import (
    parse_orca_float,
    parse_orca_int,
    split_orca_key_value,
    strip_orca_inline_comment,
)
from molop.unit import atom_ureg


pt = Chem.GetPeriodicTable()

_CARTESIAN_TYPES = {"xyz", "cart", "cartesian"}
_NON_CARTESIAN_TYPES = {"int", "internal", "gzmt"}
_EXTERNAL_TYPES = {"xyzfile", "gzmtfile", "pdbfile"}
_SUPPORTED_CTYPES = _CARTESIAN_TYPES | _NON_CARTESIAN_TYPES | _EXTERNAL_TYPES
_DUMMY_SYMBOLS = {"DA", "X", "XX"}
_GEOMETRY_CLOSE_TOKENS = {"*", "end", "edn"}


@dataclass(slots=True)
class ORCAGeometrySection:
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


def _strip_balanced_braces(token: str) -> str:
    cleaned = token.strip()
    if cleaned.startswith("{") and cleaned.endswith("}"):
        return cleaned[1:-1].strip()
    return cleaned


def parse_parameter_blocks(blocks: Sequence[ORCABlock]) -> CoordinateParameters:
    parameters: list[CoordinateParameter] = []
    for block in blocks:
        if block.name.lower() != "paras":
            continue
        for line in block.lines:
            key, value = split_orca_key_value(line.text)
            if not key or not value:
                continue
            parts = [part.strip() for part in value.split(",")]
            start = parse_orca_float(parts[0]) if parts else None
            stop = parse_orca_float(parts[1]) if len(parts) > 1 else None
            steps = parse_orca_int(parts[2]) if len(parts) > 2 else None
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
    parsed = parse_orca_float(expression)
    if parsed is not None:
        return parsed
    if expression in parameters:
        return parameters[expression]

    matched = orca_inp_patterns.VARIABLE_OFFSET.match(expression)
    if matched is None:
        return None
    base = parameters.get(matched.group("name"))
    offset = parse_orca_float(matched.group("offset"))
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

    matched = orca_inp_patterns.ATOM_TOKEN.match(cleaned)
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
        value = parse_orca_float(candidate)
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


def extract_star_geometry(block: str) -> ORCAGeometrySection | None:
    spans = build_line_spans(block)
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

        charge = parse_orca_int(header_tokens[1]) if len(header_tokens) > 1 else None
        multiplicity = parse_orca_int(header_tokens[2]) if len(header_tokens) > 2 else None
        external_path = (
            " ".join(header_tokens[3:]).strip().strip('"').strip("'")
            if ctype in _EXTERNAL_TYPES and len(header_tokens) > 3
            else None
        )
        units = (
            None
            if ctype in _EXTERNAL_TYPES
            else header_tokens[3]
            if len(header_tokens) > 3
            else None
        )

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

        return ORCAGeometrySection(
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


def extract_percent_coords_geometry(block: str) -> ORCAGeometrySection | None:
    spans = build_line_spans(block)
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
            if is_coords_open_line(candidate):
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
            if is_coords_open_line(lowered):
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

            key, value = split_orca_key_value(stripped_line)
            key_lower = key.lower()
            value_clean = value.strip().strip('"').strip("'")
            if key_lower == "ctyp":
                ctype = value_clean
            elif key_lower == "charge":
                charge = parse_orca_int(value_clean)
            elif key_lower == "mult":
                multiplicity = parse_orca_int(value_clean)
            elif key_lower == "units":
                units = value_clean

        return ORCAGeometrySection(
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
        parse_line = strip_orca_inline_comment(raw_line).strip()
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
        parse_line = strip_orca_inline_comment(raw_line).strip()
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
            parsed = parse_orca_int(row_tokens[position]) if len(row_tokens) > position else None
            if parsed is None or parsed <= 0:
                return 0
            return parsed - 1

        distance = _resolve_float_expression(tokens[4], parameters) if len(tokens) > 4 else 0.0
        angle = _resolve_float_expression(tokens[5], parameters) if len(tokens) > 5 else 0.0
        dihedral = _resolve_float_expression(tokens[6], parameters) if len(tokens) > 6 else 0.0
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


def geometry_from_section(
    section: ORCAGeometrySection | None, blocks: Sequence[ORCABlock]
) -> ORCAGeometry | None:
    if section is None:
        return None
    atoms: list[ORCAGeometryAtom] = []
    point_charges: list[dict[str, float]] = []
    internal_coords: InternalCoords | None = None
    coordinate_parameters = parse_parameter_blocks(blocks)
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
