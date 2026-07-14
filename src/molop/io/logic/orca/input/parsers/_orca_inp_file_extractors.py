from __future__ import annotations

from molop.io.base_models.source import LocatedTextBlock
from molop.io.codec_exceptions import FormatMismatchError


_CARTESIAN_TYPES = {"xyz", "cart", "cartesian"}
_NON_CARTESIAN_TYPES = {"int", "internal", "gzmt"}
_INLINE_STAR_GEOMETRY_TYPES = _CARTESIAN_TYPES | _NON_CARTESIAN_TYPES
_EXTERNAL_STAR_GEOMETRY_TYPES = {"xyzfile", "gzmtfile", "pdbfile"}
_STAR_GEOMETRY_TYPES = _INLINE_STAR_GEOMETRY_TYPES | _EXTERNAL_STAR_GEOMETRY_TYPES
_ORCA_INPUT_PROBE_CHARS = 20000


def is_orca_star_geometry_header(line: str) -> bool:
    stripped = line.strip()
    if not stripped.startswith("*"):
        return False
    header = stripped[1:].strip()
    if not header:
        return False
    ctype = header.split(maxsplit=1)[0].lower()
    return ctype in _STAR_GEOMETRY_TYPES


def is_orca_coords_open_line(line: str) -> bool:
    lower = line.lower()
    return lower.startswith("coords") and (len(lower) == 6 or lower[6].isspace())


def ensure_orca_input_content(file_content: str) -> None:
    text = file_content[:_ORCA_INPUT_PROBE_CHARS]
    lowered = text.lower()
    has_simple_input = any(line.lstrip().startswith("!") for line in text.splitlines())
    has_orca_block = any(line.lstrip().startswith("%") for line in text.splitlines())
    has_geometry = any(is_orca_star_geometry_header(line) for line in text.splitlines())
    has_geometry = has_geometry or "%coords" in lowered
    if not (has_simple_input or has_orca_block):
        raise FormatMismatchError("Not an ORCA input file: missing ORCA command or block.")
    if not has_geometry:
        raise FormatMismatchError("Not an ORCA input file: missing coordinate section.")


def is_orca_new_job_delimiter(line: str) -> bool:
    lowered = line.strip().lower()
    return lowered == "$new_job" or lowered.startswith("$new_job ")


def locate_orca_input_frames(file_content: str) -> tuple[LocatedTextBlock, ...]:
    """Locate ORCA jobs while respecting geometry blocks containing delimiter-like text."""

    lines = file_content.splitlines(keepends=True)
    offsets = [0]
    for line in lines:
        offsets.append(offsets[-1] + len(line))

    frames: list[LocatedTextBlock] = []
    start_char = 0
    in_star_geometry = False
    in_percent_coords = False
    nested_coords_depth = 0

    for line_index, line in enumerate(lines):
        stripped = line.strip()
        lowered = stripped.lower()
        if not in_star_geometry and not in_percent_coords and is_orca_new_job_delimiter(stripped):
            end_char = offsets[line_index]
            if file_content[start_char:end_char].strip():
                frames.append(LocatedTextBlock(start_char, end_char))
            start_char = offsets[line_index + 1]
            continue

        if not in_percent_coords:
            if in_star_geometry:
                if stripped == "*":
                    in_star_geometry = False
            elif stripped.startswith("*"):
                header = stripped[1:].strip()
                if header:
                    ctype = header.split(maxsplit=1)[0].lower()
                    if ctype in _INLINE_STAR_GEOMETRY_TYPES:
                        in_star_geometry = True

        if not in_star_geometry:
            if not in_percent_coords:
                if lowered.startswith(r"%coords"):
                    in_percent_coords = True
                    nested_coords_depth = 0
            else:
                if is_orca_coords_open_line(lowered):
                    nested_coords_depth += 1
                elif lowered == "end":
                    if nested_coords_depth > 0:
                        nested_coords_depth -= 1
                    else:
                        in_percent_coords = False

    if file_content[start_char:].strip():
        frames.append(LocatedTextBlock(start_char, len(file_content)))
    if not frames:
        raise FormatMismatchError("Not an ORCA input file: no non-empty jobs found.")
    return tuple(frames)
