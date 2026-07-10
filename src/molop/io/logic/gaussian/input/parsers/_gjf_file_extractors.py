from __future__ import annotations

import os

from molop.io.base_models.SearchPattern import MolOPMatch
from molop.io.codec_exceptions import FormatMismatchError
from molop.io.logic.gaussian.input.GaussianInputPatterns import g16_input_patterns


GJF_PROBE_CHARS = 20000


def ensure_gjf_content(file_content: str) -> None:
    stripped_lines = [line.strip() for line in file_content[:GJF_PROBE_CHARS].splitlines()]
    if not stripped_lines:
        raise FormatMismatchError("Not a Gaussian input file: empty file.")
    if any(line.startswith("@") and len(line) > 1 for line in stripped_lines):
        return
    if not any(line.startswith("#") for line in stripped_lines):
        raise FormatMismatchError("Not a Gaussian input file: missing route section.")
    if not any(g16_input_patterns.CHARGE_MULTIPLICITY_LINE.match(line) for line in stripped_lines):
        raise FormatMismatchError("Not a Gaussian input file: missing charge/multiplicity line.")


def expand_gjf_at_includes(file_content: str, source_path: str | None = None) -> str:
    normalized_content = file_content.replace("\r\n", "\n").replace("\r", "\n")
    base_dir = os.path.dirname(source_path) if source_path else None

    def _expand(content: str, visited: set[str]) -> str:
        expanded_lines: list[str] = []
        for raw_line in content.splitlines():
            stripped = raw_line.strip()
            if stripped.startswith("@") and len(stripped) > 1:
                include_path_token = stripped[1:].strip()
                if include_path_token.startswith(('"', "'")) and include_path_token.endswith(
                    ('"', "'")
                ):
                    include_path_token = include_path_token[1:-1].strip()
                if not include_path_token:
                    raise ValueError("Include syntax `@filename` requires a non-empty filename")

                if base_dir is None:
                    raise ValueError(
                        "Cannot resolve `@filename` include without source file path context"
                    )

                include_path = os.path.abspath(os.path.join(base_dir, include_path_token))
                if include_path in visited:
                    raise ValueError(
                        f"Detected recursive `@filename` include cycle at {include_path}"
                    )
                if not os.path.exists(include_path):
                    raise ValueError(f"Included file does not exist: {include_path}")
                if not os.path.isfile(include_path):
                    raise ValueError(f"Included path is not a file: {include_path}")

                with open(include_path) as f:
                    include_content = f.read().replace("\r\n", "\n").replace("\r", "\n")
                expanded_lines.append(_expand(include_content, visited | {include_path}))
                continue

            expanded_lines.append(raw_line)

        return "\n".join(expanded_lines)

    initial_visited = {source_path} if source_path else set()
    return _expand(normalized_content, initial_visited)


def find_gjf_link1_matches(file_content: str) -> list[MolOPMatch]:
    return g16_input_patterns.LINK1_MARKER.find_matches(file_content)


def validate_gjf_link1_boundaries(file_content: str, link1_matches: list[MolOPMatch]) -> None:
    for match in link1_matches:
        marker_start = match.start()
        if marker_start < 2 or file_content[marker_start - 2 : marker_start] != "\n\n":
            raise ValueError("`--Link1--` must be preceded by a blank line")


def split_gjf_link1_frames(file_content: str, link1_matches: list[MolOPMatch]) -> list[str]:
    frames: list[str] = []
    start = 0
    for match in link1_matches:
        frames.append(file_content[start : match.start()])
        start = match.end()
    frames.append(file_content[start:])
    return [frame.strip() + "\n\n" for frame in frames]
