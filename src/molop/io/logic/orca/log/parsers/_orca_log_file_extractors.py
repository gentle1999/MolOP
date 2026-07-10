from __future__ import annotations

from collections.abc import Sequence
from typing import Any

from molop.io.codec_exceptions import FormatMismatchError
from molop.io.logic.orca.input.parsers._orca_inp_metadata import parse_orca_input_metadata
from molop.io.logic.orca.log.parsers._orca_log_patterns import orca_log_patterns


_ORCA_PROBE_BYTES = 12000


def ensure_orca_output_content(file_content: str) -> None:
    prefix = file_content[:_ORCA_PROBE_BYTES]
    if "ORCA" not in prefix:
        raise FormatMismatchError("Not an ORCA output file.")
    if "Program Version" not in file_content:
        raise FormatMismatchError("Not an ORCA output file: missing Program Version.")
    if orca_log_patterns.BANNER.search(prefix) is None and "* O   R   C   A *" not in prefix:
        raise FormatMismatchError("Not an ORCA output file: missing ORCA banner.")


def extract_orca_output_version(file_content: str) -> str | None:
    if matched := orca_log_patterns.VERSION.search(file_content):
        return matched.group("version")
    return None


def extract_orca_printed_input(file_content: str) -> tuple[str, str]:
    blocks = orca_log_patterns.INPUT_BLOCK.find_matches(file_content)
    if not blocks:
        return "", ""
    body = blocks[-1].group("body")
    name = ""
    if matched := orca_log_patterns.INPUT_NAME.search(body):
        name = matched.group("name").strip()
    lines: list[str] = []
    for raw_line in body.splitlines():
        matched = orca_log_patterns.INPUT_LINE.match(raw_line)
        if matched is None:
            continue
        line = matched.group("line")
        if "****END OF INPUT****" in line:
            break
        lines.append(line.rstrip())
    return name, "\n".join(lines).strip() + ("\n" if lines else "")


def parse_orca_printed_input_metadata(input_text: str) -> dict[str, Any]:
    return parse_orca_input_metadata(input_text)


def split_orca_output_frames(file_content: str) -> list[str]:
    matches = orca_log_patterns.COORD_HEADER.find_matches(file_content)
    if not matches:
        return [file_content]
    frames: list[str] = []
    prefix = file_content[: matches[0].start()]
    for idx, matched in enumerate(matches):
        next_start = matches[idx + 1].start() if idx + 1 < len(matches) else len(file_content)
        start = matched.start()
        if idx == 0:
            frames.append(file_content[:next_start])
        else:
            frames.append(prefix + file_content[start:next_start])
    return [frame for frame in frames if frame.strip()]


def last_frame_value(frames: Sequence[Any], field: str) -> Any:
    for frame in reversed(frames):
        value = getattr(frame, field, None)
        if value is not None:
            return value
    return None


def first_frame_value(frames: Sequence[Any], field: str) -> Any:
    for frame in frames:
        value = getattr(frame, field, None)
        if value is not None:
            return value
    return None
