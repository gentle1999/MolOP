from __future__ import annotations

from molop.io.logic.gaussian.input.GaussianInputPatterns import g16_input_patterns


def normalize_gjf_free_separators(line: str) -> str:
    return line.replace("\t", " ").replace(",", " ").replace("/", " ")


def strip_gjf_comment_content(line: str) -> str:
    return line.split("!", 1)[0].rstrip()


def build_gjf_context_lines(block: str) -> list[str]:
    return [strip_gjf_comment_content(line) for line in block.splitlines()]


def is_gjf_blank(line: str) -> bool:
    return strip_gjf_comment_content(line).strip() == ""


def join_gjf_lines(lines: list[str], start: int, end: int) -> str:
    if start >= end:
        return ""
    return "\n".join(lines[start:end])


def find_first_gjf_route_index(lines: list[str]) -> int | None:
    for idx, line in enumerate(lines):
        if strip_gjf_comment_content(line).lstrip().startswith("#"):
            return idx
    return None


def find_next_gjf_blank_index(lines: list[str], start: int) -> int | None:
    for idx in range(start, len(lines)):
        if is_gjf_blank(lines[idx]):
            return idx
    return None


def skip_gjf_blank_indices(lines: list[str], start: int) -> int:
    idx = start
    while idx < len(lines) and is_gjf_blank(lines[idx]):
        idx += 1
    return idx


def is_gjf_charge_multiplicity_line(line: str) -> bool:
    stripped = " ".join(normalize_gjf_free_separators(strip_gjf_comment_content(line)).split())
    if not stripped:
        return False
    parts = stripped.split()
    if len(parts) < 2 or len(parts) % 2 != 0:
        return False
    return all(g16_input_patterns.INTEGER_TOKEN.match(part) is not None for part in parts)


def find_first_gjf_charge_multiplicity_index(lines: list[str], start: int) -> int | None:
    for idx in range(start, len(lines)):
        if is_gjf_charge_multiplicity_line(lines[idx]):
            return idx
    return None


def is_gjf_zmat_variable_label(line: str) -> bool:
    stripped = strip_gjf_comment_content(line).strip().lower()
    return stripped in {"variables:", "constants:"}


def is_gjf_zmat_variable_assignment(line: str) -> bool:
    normalized = " ".join(normalize_gjf_free_separators(strip_gjf_comment_content(line)).split())
    return g16_input_patterns.ZMAT_VARIABLE_ASSIGNMENT.match(normalized) is not None


def gjf_route_flags(route_raw: str) -> tuple[bool, bool]:
    normalized = "".join(route_raw.lower().split())
    has_allcheck = "allcheck" in normalized
    has_checkpoint = any(
        key in normalized
        for key in ("geom=check", "geom=checkpoint", "geom=(check", "geom=(checkpoint")
    )
    return has_allcheck, has_checkpoint


def find_gjf_molecule_block_end(lines: list[str], charge_idx: int) -> tuple[int, int]:
    first_blank = find_next_gjf_blank_index(lines, charge_idx)
    molecule_end = first_blank if first_blank is not None else len(lines)
    additional_start = molecule_end + 1 if first_blank is not None else molecule_end

    probe_idx = skip_gjf_blank_indices(lines, additional_start)
    if probe_idx >= len(lines):
        return molecule_end, additional_start
    if not (
        is_gjf_zmat_variable_label(lines[probe_idx])
        or is_gjf_zmat_variable_assignment(lines[probe_idx])
    ):
        return molecule_end, additional_start

    scan_idx = probe_idx
    seen_assignment = False
    while scan_idx < len(lines):
        stripped = lines[scan_idx].strip()
        if not stripped:
            scan_idx += 1
            continue
        if is_gjf_zmat_variable_label(stripped):
            scan_idx += 1
            continue
        if is_gjf_zmat_variable_assignment(stripped):
            seen_assignment = True
            scan_idx += 1
            continue
        break

    if seen_assignment:
        return scan_idx, scan_idx
    return molecule_end, additional_start
