from __future__ import annotations

from molop.io.base_models.SearchPattern import MolOPMatch
from molop.io.base_models.source import LocatedTextBlock
from molop.io.codec_exceptions import FormatMismatchError, ParseError
from molop.io.logic.gaussian.input.GaussianInputPatterns import g16_input_patterns


GJF_PROBE_CHARS = 20000
GJF_INCLUDE_PROVENANCE_DIAGNOSTIC = "MOL.PARSE.GJF_INCLUDE_PROVENANCE_UNSUPPORTED"


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


def ensure_gjf_single_artifact_source(file_content: str) -> None:
    """Reject includes until provenance can identify every expanded artifact."""

    if any(line.strip().startswith("@") for line in file_content.splitlines()):
        raise ParseError(
            f"{GJF_INCLUDE_PROVENANCE_DIAGNOSTIC}: Gaussian @include expansion "
            "cannot be represented by single-artifact SourceSpan values."
        )


def find_gjf_link1_matches(file_content: str) -> list[MolOPMatch]:
    return g16_input_patterns.LINK1_MARKER.find_matches(file_content)


def validate_gjf_link1_boundaries(file_content: str, link1_matches: list[MolOPMatch]) -> None:
    for match in link1_matches:
        marker_start = match.start()
        prefix_lines = file_content[:marker_start].rstrip(" \t\f\v").splitlines(keepends=True)
        if not prefix_lines or prefix_lines[-1].strip():
            raise ValueError("`--Link1--` must be preceded by a blank line")


def locate_gjf_link1_frames(file_content: str) -> tuple[LocatedTextBlock, ...]:
    """Locate exact Link1 input ranges within one source artifact."""

    ensure_gjf_single_artifact_source(file_content)
    matches = find_gjf_link1_matches(file_content)
    validate_gjf_link1_boundaries(file_content, matches)
    starts = [0, *(matched.end() for matched in matches)]
    ends = [*(matched.start() for matched in matches), len(file_content)]
    ranges = tuple(zip(starts, ends, strict=True))
    if any(not file_content[start:end].strip() for start, end in ranges):
        raise ValueError("`--Link1--` must delimit non-empty Gaussian input jobs")
    return tuple(LocatedTextBlock(start, end) for start, end in ranges)
