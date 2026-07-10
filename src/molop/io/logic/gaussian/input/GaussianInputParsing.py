from __future__ import annotations

from collections.abc import Mapping
from typing import TypeAlias

from molop.io.logic.gaussian.input.GaussianInput import (
    GJFAtomSpecification,
    GJFGICLine,
    GJFGICOption,
    GJFGICSection,
    GJFModRedundantLine,
    GJFModRedundantSection,
    GJFMoleculeSpecifications,
    GJFMoleculeSpecificationsFragment,
    GJFNBOSection,
    GJFRouteSection,
    GJFSectionParsingDiagnostic,
    GJFTitleCard,
    GJFUnknownSection,
)
from molop.io.logic.gaussian.input.GaussianInputPatterns import g16_input_patterns
from molop.io.logic.gaussian.input.GaussianLink0 import GJFLink0, GJFLink0Commands
from molop.io.logic.gaussian.input.GaussianRouteParsing import parse_gaussian_route_semantic


GJFParsedAdditionalSection: TypeAlias = (
    GJFGICSection | GJFModRedundantSection | GJFNBOSection | GJFUnknownSection
)


def parse_gjf_link0_commands(raw: str) -> GJFLink0Commands:
    link0_keywords: list[GJFLink0] = []
    for matched in g16_input_patterns.OPTIONS.find_matches(raw):
        link0_keywords.append(GJFLink0(key=matched.group("key"), value=matched.group("value")))
    return GJFLink0Commands(
        link0_keywords=link0_keywords,
    )


def _normalize_route(raw: str) -> str:
    route = raw
    if not route.startswith("#"):
        route = f"# {route}"
    return route.replace(" \n", "").replace("\n", "")


def parse_gjf_route_section(raw: str) -> GJFRouteSection:
    route = _normalize_route(raw)
    return GJFRouteSection(
        route=route,
        semantic_route=parse_gaussian_route_semantic(route),
    )


def _normalize_title_card(raw: str) -> str:
    title = raw
    for char in GJFTitleCard.INVALID_CHARS:
        title = title.replace(char, " ").strip("\n")
    if title and len(title.splitlines()) > 5:
        raise ValueError("Title card cannot exceed 5 lines")
    return title


def parse_gjf_title_card(raw: str) -> GJFTitleCard:
    return GJFTitleCard(title_card=_normalize_title_card(raw))


def _normalize_free_separators(data: str) -> str:
    normalized_chars: list[str] = []
    depth = 0
    for ch in data:
        if ch in "([{":
            depth += 1
            normalized_chars.append(ch)
            continue
        if ch in ")]}":
            depth = max(depth - 1, 0)
            normalized_chars.append(ch)
            continue
        if depth == 0 and ch in {",", "/", "\t"}:
            normalized_chars.append(" ")
            continue
        normalized_chars.append(ch)
    return "".join(normalized_chars)


def parse_gjf_atom_specification(raw: str) -> GJFAtomSpecification:
    normalized = _normalize_free_separators(raw)
    parts = normalized.split()
    if len(parts) < 2:
        return GJFAtomSpecification(element_label=parts[0], coords_part="")

    if g16_input_patterns.CARTESIAN_COORD_START.match(parts[1]):
        coords_part = " ".join(parts[1:])
        frozen_tag = None
    else:
        frozen_tag = parts[1]
        coords_part = " ".join(parts[2:])

    atom_info = parts[0]
    matched = g16_input_patterns.ATOM_SPECIFICATION.match(atom_info)
    if matched is None:
        raise ValueError("Atom specification must match the pattern")

    element_label = matched.group("element_label")
    atom_type = matched.group("atom_type")
    charge = matched.group("charge")
    params_str = matched.group("params")
    params = (
        {
            section.split("=")[0].lower(): section.split("=")[1].lower()
            for section in params_str.split(",")
        }
        if params_str
        else {}
    )
    return GJFAtomSpecification(
        element_label=element_label,
        atom_type=atom_type,
        charge=float(charge) if charge else None,
        params=params,
        frozen_tag=frozen_tag,
        coords_part=coords_part,
    )


def _strip_comment_content(line: str) -> str:
    return line.split("!", 1)[0].rstrip()


def _normalize_molecule_free_separators(line: str) -> str:
    return line.replace("\t", " ").replace(",", " ").replace("/", " ")


def _is_zmat_variable_label(line: str) -> bool:
    return line.strip().lower() in {"variables:", "constants:"}


def _parse_zmat_variable_assignment(line: str) -> tuple[str, str] | None:
    normalized = " ".join(_normalize_molecule_free_separators(line).split())
    matched = g16_input_patterns.ZMAT_VARIABLE_ASSIGNMENT.match(normalized)
    if matched is None:
        return None
    return matched.group("name"), matched.group("value")


def _split_atom_and_variable_lines(lines: list[str]) -> tuple[list[str], dict[str, str]]:
    atom_lines: list[str] = []
    variable_map: dict[str, str] = {}
    in_variable_block = False
    blank_seen_after_atoms = False

    for line in lines:
        stripped = line.strip()
        if not stripped:
            if atom_lines:
                blank_seen_after_atoms = True
            continue

        if _is_zmat_variable_label(stripped):
            in_variable_block = True
            continue

        assignment = _parse_zmat_variable_assignment(stripped)
        if assignment is not None and (in_variable_block or blank_seen_after_atoms):
            var_name, var_value = assignment
            variable_map[var_name] = var_value
            in_variable_block = True
            continue

        atom_lines.append(stripped)

    return atom_lines, variable_map


def _replace_zmat_variables_in_line(line: str, variable_map: Mapping[str, str]) -> str:
    normalized_line = _normalize_molecule_free_separators(line)
    if not variable_map:
        return " ".join(normalized_line.split())
    tokens = normalized_line.split()
    if len(tokens) < 2:
        return " ".join(tokens)

    replaced: list[str] = [tokens[0]]
    for token in tokens[1:]:
        if token in variable_map:
            replaced.append(variable_map[token])
            continue
        if len(token) >= 2 and token[0] in "+-" and token[1:] in variable_map:
            replaced.append(f"{token[0]}{variable_map[token[1:]]}")
            continue
        replaced.append(token)
    return " ".join(replaced)


def _validate_internal_coordinate_references(
    atom_specifications: list[GJFAtomSpecification],
) -> None:
    for atom_index, atom_spec in enumerate(atom_specifications, start=1):
        if not atom_spec.coords_part.strip() or atom_spec.is_cartesian_coords():
            continue
        if not atom_spec.is_internal_coords():
            continue

        values = atom_spec.coords_part.split()
        if len(values) % 2 == 1 and len(values) >= 3 and values[-1] in {"0", "1"}:
            values = values[:-1]
        pair_count = len(values) // 2
        refs = [int(values[i]) for i in range(0, len(values), 2)]

        for ref_pos, ref in enumerate(refs, start=1):
            if ref <= 0:
                raise ValueError(
                    f"Z-matrix reference {ref_pos} on atom {atom_index} must be a positive 1-based index"
                )
            if ref >= atom_index:
                raise ValueError(
                    f"Z-matrix reference {ref_pos} on atom {atom_index} must refer to a previously defined atom"
                )

        if pair_count >= 2 and refs[0] == refs[1]:
            raise ValueError(
                f"Z-matrix distance and angle references on atom {atom_index} must be different"
            )
        if pair_count >= 3 and len(set(refs[:3])) < 3:
            raise ValueError(
                f"Z-matrix distance/angle/third references on atom {atom_index} must be distinct"
            )


def parse_gjf_molecule_specifications(raw: str) -> GJFMoleculeSpecifications:
    lines = [_strip_comment_content(line) for line in raw.splitlines()]
    if "".join(lines).strip() == "":
        return GJFMoleculeSpecifications()

    first_non_blank_idx = next((idx for idx, line in enumerate(lines) if line.strip()), None)
    if first_non_blank_idx is None:
        return GJFMoleculeSpecifications()

    body_lines = lines[first_non_blank_idx + 1 :]
    if len(body_lines) < 1:
        raise ValueError(
            "Molecule specifications must have at least 2 lines(charge and spin multiplicity line"
            " & coordinates line)"
        )

    electron_config = _normalize_molecule_free_separators(lines[first_non_blank_idx]).split()
    if len(electron_config) < 2 or len(electron_config) % 2 != 0:
        raise ValueError("charge and spin multiplicity line must have even number of values")

    atom_lines, variable_map = _split_atom_and_variable_lines(body_lines)
    normalized_atom_lines = [
        _replace_zmat_variables_in_line(line, variable_map) for line in atom_lines
    ]

    if len(electron_config) == 2:
        total_charge, spin_multiplicity = map(int, electron_config)
        atom_specifications = [parse_gjf_atom_specification(line) for line in normalized_atom_lines]
        _validate_internal_coordinate_references(atom_specifications)
        molecule_fragments = [
            GJFMoleculeSpecificationsFragment(
                total_charge=total_charge,
                spin_multiplicity=spin_multiplicity,
                atom_specifications=atom_specifications,
            )
        ]
    else:
        total_charge, spin_multiplicity = map(int, electron_config[0:2])
        fragment_charge_spin_multiplicity: dict[int, tuple[int, int]] = {}
        for fragment_id, i in enumerate(range(2, len(electron_config), 2)):
            fragment_charge_spin_multiplicity[fragment_id] = (
                int(electron_config[i]),
                int(electron_config[i + 1]),
            )
        atom_specifications = [parse_gjf_atom_specification(line) for line in normalized_atom_lines]
        _validate_internal_coordinate_references(atom_specifications)
        fragment_ids = [atom_spec.get_fragment_id() for atom_spec in atom_specifications]
        declared_fragment_count = len(fragment_charge_spin_multiplicity)

        if any(fragment_id <= 0 for fragment_id in fragment_ids):
            raise ValueError(
                "Multi-fragment molecule specifications require every atom to declare Fragment=n"
            )

        expected_fragment_ids = set(range(1, declared_fragment_count + 1))
        actual_fragment_ids = set(fragment_ids)
        if actual_fragment_ids != expected_fragment_ids:
            raise ValueError(
                "Fragment assignments must be contiguous and match declared fragment charge/spin pairs"
            )

        molecule_fragments = [
            GJFMoleculeSpecificationsFragment(
                fragment_id=fragment_id,
                total_charge=sub_charge,
                spin_multiplicity=sub_spin_multiplicity,
                atom_specifications=[
                    atom_specification
                    for atom_specification in atom_specifications
                    if atom_specification.params.get("fragment", "0") == str(fragment_id + 1)
                ],
            )
            for fragment_id, (
                sub_charge,
                sub_spin_multiplicity,
            ) in fragment_charge_spin_multiplicity.items()
        ]
        if any(len(fragment.atom_specifications) == 0 for fragment in molecule_fragments):
            raise ValueError(
                "Each declared fragment charge/spin pair must correspond to at least one atom"
            )

    return GJFMoleculeSpecifications(
        total_charge=total_charge,
        spin_multiplicity=spin_multiplicity,
        molecule_fragments=molecule_fragments,
    )


def _split_top_level_args(raw: str) -> list[str]:
    if not raw.strip():
        return []
    args: list[str] = []
    current: list[str] = []
    depth = 0
    for ch in raw:
        if ch == "," and depth == 0:
            token = "".join(current).strip()
            if token:
                args.append(token)
            current = []
            continue
        if ch in "([{":
            depth += 1
        elif ch in ")]}":
            depth = max(depth - 1, 0)
        current.append(ch)
    token = "".join(current).strip()
    if token:
        args.append(token)
    return args


def _parse_gic_function_expression(expression: str) -> tuple[str | None, list[str]]:
    matched = g16_input_patterns.GIC_FUNCTION_EXPRESSION.match(expression.strip())
    if matched is None:
        return None, []
    return matched.group("name"), _split_top_level_args(matched.group("args"))


def parse_gjf_gic_option(token: str) -> GJFGICOption:
    raw = token.strip()
    if "=" in raw:
        key, value = raw.split("=", 1)
        return GJFGICOption(
            raw=raw,
            key=key.strip().lower(),
            value=value.strip(),
            is_flag=False,
        )
    return GJFGICOption(raw=raw, key=raw.lower(), value=None, is_flag=True)


def _derive_gic_state_from_options(
    options: list[GJFGICOption], standalone_action: str | None = None
) -> tuple[str | None, dict[str, str]]:
    option_values = {opt.key: opt.value for opt in options if opt.value is not None}
    normalized_state = None

    if standalone_action is not None:
        action = standalone_action.lower()
        if action in {"freeze", "frozen"}:
            normalized_state = "frozen"
        elif action in {"remove", "inactive", "kill", "removeall"}:
            normalized_state = "inactive"
        elif action in {"active", "activate", "modify"}:
            normalized_state = "active"
        elif action in {"printonly", "print-only"}:
            normalized_state = "print-only"

    if normalized_state is None:
        option_keys = {opt.key for opt in options}
        if option_keys & {"freeze", "frozen"}:
            normalized_state = "frozen"
        elif option_keys & {"inactive", "remove", "kill", "removeall"}:
            normalized_state = "inactive"
        elif option_keys & {"active", "modify"}:
            normalized_state = "active"
        elif option_keys & {"printonly", "print-only"}:
            normalized_state = "print-only"

    return normalized_state, option_values


def _split_top_level_assignment(line: str) -> tuple[str, str] | None:
    depth = 0
    for idx, ch in enumerate(line):
        if ch in "([{":
            depth += 1
            continue
        if ch in ")]}":
            depth = max(depth - 1, 0)
            continue
        if ch == "=" and depth == 0:
            return line[:idx].strip(), line[idx + 1 :].strip()
    return None


_GIC_STANDALONE_OPTIONS = {"freezeall", "unfreezeall", "removeall"}


def parse_gjf_gic_line(raw: str) -> GJFGICLine:
    line = raw.strip()
    if not line:
        raise ValueError("GIC line cannot be empty")
    line = line.split("!", 1)[0].strip()
    lower_line = line.lower()
    if lower_line in _GIC_STANDALONE_OPTIONS or g16_input_patterns.GIC_STANDALONE_ATOM.match(
        lower_line
    ):
        parts = line.split()
        standalone_action = parts[2].lower() if len(parts) > 2 else parts[0].lower()
        normalized_state, option_values = _derive_gic_state_from_options([], standalone_action)
        return GJFGICLine(
            raw_line=raw,
            expression=line,
            is_standalone_option=True,
            expression_kind="standalone",
            standalone_keyword=parts[0],
            standalone_args=parts[1:],
            standalone_target=parts[1] if len(parts) > 1 else None,
            standalone_action=standalone_action,
            normalized_state=normalized_state,
            option_values=option_values,
        )

    assignment = _split_top_level_assignment(line)
    if assignment is None:
        function_name, function_args = _parse_gic_function_expression(line)
        return GJFGICLine(
            raw_line=raw,
            expression=line,
            expression_kind="function" if function_name else "raw",
            function_name=function_name,
            function_args=function_args,
        )

    left, right = assignment
    left_match = g16_input_patterns.GIC_LEFT_ASSIGNMENT.match(left)
    if left_match is None:
        function_name, function_args = _parse_gic_function_expression(right)
        return GJFGICLine(
            raw_line=raw,
            expression=line,
            expression_kind="function" if function_name else "raw",
            function_name=function_name,
            function_args=function_args,
        )

    opts = left_match.group("opts")
    function_name, function_args = _parse_gic_function_expression(right)
    parsed_options = (
        [parse_gjf_gic_option(opt) for opt in _split_top_level_args(opts)] if opts else []
    )
    normalized_state, option_values = _derive_gic_state_from_options(parsed_options)
    return GJFGICLine(
        raw_line=raw,
        label=left_match.group("label"),
        label_options=_split_top_level_args(opts) if opts else [],
        expression=right,
        expression_kind="function" if function_name else "assignment",
        function_name=function_name,
        function_args=function_args,
        parsed_label_options=parsed_options,
        normalized_state=normalized_state,
        option_values=option_values,
    )


_GIC_HEAD_TOKENS = (
    "r(",
    "bond(",
    "stretch(",
    "a(",
    "angle(",
    "bend(",
    "d(",
    "dihedral(",
    "torsion(",
    "l(",
    "linear(",
    "linearbend(",
    "x(",
    "y(",
    "z(",
    "cartesian(",
    "cart(",
    "dotdiff(",
    "xcntr(",
    "ycntr(",
    "zcntr(",
    "freezeall",
    "unfreezeall",
    "removeall",
    "atom ",
)


def is_gjf_gic_section(raw: str) -> bool:
    lines = [line.strip() for line in raw.splitlines() if line.strip()]
    if not lines:
        return False
    for line in lines:
        lowered = line.lower()
        if lowered.startswith(_GIC_HEAD_TOKENS) or "=" in lowered:
            return True
    return False


def parse_gjf_gic_section(raw: str) -> GJFGICSection:
    parsed_lines = [parse_gjf_gic_line(line) for line in raw.splitlines() if line.strip()]
    return GJFGICSection(raw=raw, lines=parsed_lines)


def parse_gjf_modredundant_line(raw: str) -> GJFModRedundantLine:
    line = raw.strip()
    if not line:
        raise ValueError("ModRedundant line cannot be empty")
    parts = line.split()
    if len(parts) <= 1:
        return GJFModRedundantLine(raw_line=raw, coordinate_type=parts[0])

    coordinate_type = parts[0]
    atom_refs: list[str] = []
    idx = 1
    while idx < len(parts) and g16_input_patterns.MODREDUNDANT_ATOM_REF.match(parts[idx]):
        atom_refs.append(parts[idx])
        idx += 1
    action = parts[idx] if idx < len(parts) else None
    parameters = parts[idx + 1 :] if idx + 1 < len(parts) else []
    return GJFModRedundantLine(
        raw_line=raw,
        coordinate_type=coordinate_type,
        atom_refs=atom_refs,
        action=action,
        parameters=parameters,
    )


_MODREDUNDANT_ACTIONS = {"a", "f", "b", "k", "r", "d", "h", "s"}


def is_gjf_modredundant_section(raw: str) -> bool:
    lines = [line.strip() for line in raw.splitlines() if line.strip()]
    if not lines:
        return False
    for line in lines:
        parts = line.split()
        if len(parts) < 2:
            continue
        if not g16_input_patterns.MODREDUNDANT_COORDINATE_TYPE.match(parts[0]):
            continue
        if any(part.lower() in _MODREDUNDANT_ACTIONS for part in parts[1:]):
            return True
    return False


def parse_gjf_modredundant_section(raw: str) -> GJFModRedundantSection:
    parsed_lines = [parse_gjf_modredundant_line(line) for line in raw.splitlines() if line.strip()]
    return GJFModRedundantSection(raw=raw, lines=parsed_lines)


def parse_gjf_unknown_section(raw: str) -> GJFUnknownSection:
    return GJFUnknownSection(raw=raw)


def is_gjf_nbo_section(raw: str) -> bool:
    lines = [line.strip() for line in raw.splitlines() if line.strip()]
    if len(lines) < 2:
        return False
    return lines[0].lower().startswith("$nbo") and lines[-1].lower() == "$end"


def parse_gjf_nbo_section(raw: str) -> GJFNBOSection:
    lines = [line.rstrip() for line in raw.splitlines() if line.strip()]
    if not lines:
        raise ValueError("NBO section cannot be empty")
    header = lines[0]
    footer = lines[-1] if len(lines) > 1 else "$end"
    commands = lines[1:-1] if len(lines) > 2 else []
    return GJFNBOSection(raw=raw, header=header, commands=commands, footer=footer)


def split_gjf_additional_sections(raw: str) -> list[str]:
    stripped = raw.strip("\n")
    if not stripped.strip():
        return []
    matches = g16_input_patterns.ADDITIONAL_SECTION_SEPARATOR.find_matches(stripped)
    sections: list[str] = []
    start = 0
    for matched in matches:
        sections.append(stripped[start : matched.start()])
        start = matched.end()
    sections.append(stripped[start:])
    return [section for section in sections if section.strip()]


def parse_gjf_additional_sections(
    raw: str,
) -> tuple[list[GJFParsedAdditionalSection], list[GJFSectionParsingDiagnostic]]:
    parsed_sections: list[GJFParsedAdditionalSection] = []
    diagnostics: list[GJFSectionParsingDiagnostic] = []
    for section_index, section_raw in enumerate(split_gjf_additional_sections(raw)):
        looks_like_gic = is_gjf_gic_section(section_raw)
        looks_like_modredundant = is_gjf_modredundant_section(section_raw)
        looks_like_nbo = is_gjf_nbo_section(section_raw)

        if looks_like_nbo:
            try:
                parsed_sections.append(parse_gjf_nbo_section(section_raw))
            except Exception as exc:
                parsed_sections.append(parse_gjf_unknown_section(section_raw))
                diagnostics.append(
                    GJFSectionParsingDiagnostic(
                        section_type="nbo",
                        message=f"Failed to parse as NBO section: {exc}",
                        section_index=section_index,
                    )
                )
            continue

        if looks_like_gic and looks_like_modredundant:
            parsed_sections.append(parse_gjf_unknown_section(section_raw))
            diagnostics.append(
                GJFSectionParsingDiagnostic(
                    section_type="mixed-additional-section",
                    message=(
                        "A single additional section mixes GIC and ModRedundant syntax, "
                        "which Gaussian does not allow"
                    ),
                    section_index=section_index,
                )
            )
            continue

        if looks_like_gic:
            try:
                parsed_sections.append(parse_gjf_gic_section(section_raw))
            except Exception as exc:
                parsed_sections.append(parse_gjf_unknown_section(section_raw))
                diagnostics.append(
                    GJFSectionParsingDiagnostic(
                        section_type="gic",
                        message=f"Failed to parse as GIC section: {exc}",
                        section_index=section_index,
                    )
                )
            continue

        if looks_like_modredundant:
            try:
                parsed_sections.append(parse_gjf_modredundant_section(section_raw))
            except Exception as exc:
                parsed_sections.append(parse_gjf_unknown_section(section_raw))
                diagnostics.append(
                    GJFSectionParsingDiagnostic(
                        section_type="modredundant",
                        message=f"Failed to parse as ModRedundant section: {exc}",
                        section_index=section_index,
                    )
                )
            continue

        parsed_sections.append(parse_gjf_unknown_section(section_raw))

    return parsed_sections, diagnostics
