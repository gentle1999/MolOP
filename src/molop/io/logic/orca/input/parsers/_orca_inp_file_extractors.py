from __future__ import annotations

from collections.abc import Sequence

from molop.io.codec_exceptions import FormatMismatchError


_CARTESIAN_TYPES = {"xyz", "cart", "cartesian"}
_NON_CARTESIAN_TYPES = {"int", "internal", "gzmt"}
_FILE_SPLIT_STAR_GEOMETRY_CTYPES = (
    _CARTESIAN_TYPES
    | _NON_CARTESIAN_TYPES
    | {
        "xyzfile",
        "gzmtfile",
    }
)
_ORCA_INPUT_PROBE_CHARS = 20000


def is_orca_star_geometry_header(line: str) -> bool:
    stripped = line.strip()
    if not stripped.startswith("*"):
        return False
    header = stripped[1:].strip()
    if not header:
        return False
    ctype = header.split(maxsplit=1)[0].lower()
    return ctype in _FILE_SPLIT_STAR_GEOMETRY_CTYPES or ctype == "pdbfile"


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


def build_orca_input_file_lines(file_content: str) -> list[str]:
    return file_content.splitlines(keepends=True)


def is_orca_new_job_delimiter(line: str) -> bool:
    lowered = line.strip().lower()
    return lowered == "$new_job" or lowered.startswith("$new_job ")


def split_orca_input_frames(lines: Sequence[str]) -> list[str]:
    frames: list[str] = []
    current_lines: list[str] = []
    in_star_geometry = False
    in_percent_coords = False
    nested_coords_depth = 0

    for line in lines:
        stripped = line.strip()
        lowered = stripped.lower()

        if not in_star_geometry and not in_percent_coords and is_orca_new_job_delimiter(stripped):
            frames.append("".join(current_lines))
            current_lines = []
            continue

        current_lines.append(line)

        if not in_percent_coords:
            if in_star_geometry:
                if stripped == "*":
                    in_star_geometry = False
            elif stripped.startswith("*"):
                header = stripped[1:].strip()
                if header:
                    ctype = header.split(maxsplit=1)[0].lower()
                    if ctype in _FILE_SPLIT_STAR_GEOMETRY_CTYPES:
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

    frames.append("".join(current_lines))
    return [frame for frame in frames if frame]
