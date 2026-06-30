from __future__ import annotations

import re
from collections.abc import Sequence
from typing import TYPE_CHECKING, Any, Protocol, cast

from molop.io.base_models.DataClasses import ImplicitSolvation, Status
from molop.io.base_models.FileParser import BaseFileParserDisk, BaseFileParserMemory
from molop.io.codec_exceptions import FormatMismatchError
from molop.io.logic.QM_frame_models.ORCALogFileFrame import (
    ORCALogFileFrameDisk,
    ORCALogFileFrameMemory,
)
from molop.io.logic.QM_frame_parsers.ORCALogFileFrameParser import (
    ORCALogFileFrameParserDisk,
    ORCALogFileFrameParserMemory,
)
from molop.io.logic.QM_models.ORCALogFile import ORCALogFileDisk, ORCALogFileMemory
from molop.io.logic.qminput_frame_parsers.ORCAInpFileFrameParser import (
    ORCAInpFileFrameParserMemory,
)
from molop.unit import atom_ureg


if TYPE_CHECKING:
    from molop.io.codec_registry import Registry


_ORCA_BANNER_RE = re.compile(r"\*\s+O\s+R\s+C\s+A\s+\*", re.I)
_ORCA_PROBE_BYTES = 12000
_VERSION_RE = re.compile(r"Program Version\s+([0-9][^\s]*)")
_INPUT_BLOCK_RE = re.compile(
    r"=+\s*\n\s*INPUT FILE\s*\n=+\s*\n(?P<body>.*?)\n=+",
    re.S,
)
_INPUT_NAME_RE = re.compile(r"^\s*NAME\s*=\s*(?P<name>.+?)\s*$", re.M)
_INPUT_LINE_RE = re.compile(r"^\|\s*\d+>\s?(?P<line>.*)$")
_COORD_HEADER_RE = re.compile(
    r"^-+\s*\nCARTESIAN COORDINATES \(ANGSTROEM\)\s*\n-+\s*$",
    re.MULTILINE,
)
_RUN_TIME_RE = re.compile(
    r"TOTAL RUN TIME:\s*"
    r"(?P<days>\d+)\s+days\s+"
    r"(?P<hours>\d+)\s+hours\s+"
    r"(?P<minutes>\d+)\s+minutes?\s+"
    r"(?P<seconds>\d+)\s+seconds?\s+"
    r"(?P<msec>\d+)\s+msec",
    re.I,
)


class _HasMetadataFinalize(Protocol):
    def _update_file_metadata_from_frames(self, chem_file: Any, metadata: dict[str, Any]) -> None: ...


def _ensure_orca_output(file_content: str) -> None:
    prefix = file_content[:_ORCA_PROBE_BYTES]
    if "ORCA" not in prefix:
        raise FormatMismatchError("Not an ORCA output file.")
    if "Program Version" not in file_content:
        raise FormatMismatchError("Not an ORCA output file: missing Program Version.")
    if _ORCA_BANNER_RE.search(prefix) is None and "* O   R   C   A *" not in prefix:
        raise FormatMismatchError("Not an ORCA output file: missing ORCA banner.")


def _parse_version(file_content: str) -> str | None:
    if matched := _VERSION_RE.search(file_content):
        return matched.group(1)
    return None


def _extract_printed_input(file_content: str) -> tuple[str, str]:
    blocks = list(_INPUT_BLOCK_RE.finditer(file_content))
    if not blocks:
        return "", ""
    body = blocks[-1].group("body")
    name = ""
    if matched := _INPUT_NAME_RE.search(body):
        name = matched.group("name").strip()
    lines: list[str] = []
    for raw_line in body.splitlines():
        matched = _INPUT_LINE_RE.match(raw_line)
        if matched is None:
            continue
        line = matched.group("line")
        if "****END OF INPUT****" in line:
            break
        lines.append(line.rstrip())
    return name, "\n".join(lines).strip() + ("\n" if lines else "")


def _parse_printed_input_metadata(input_text: str) -> dict[str, Any]:
    if not input_text.strip():
        return {}
    parser = ORCAInpFileFrameParserMemory()
    frame = parser.parse(input_text)
    metadata: dict[str, Any] = {}
    for field in (
        "keywords",
        "method",
        "basis_set",
        "functional",
        "model_chemistry",
        "task_requests",
        "excited_state_requests",
        "multireference_requests",
        "resources_raw",
        "request_num_cpu",
        "request_memory",
        "resource_request",
        "charge",
        "multiplicity",
        "auxiliary_basis_set",
        "dispersion_correction",
    ):
        if hasattr(frame, field):
            metadata[field] = getattr(frame, field)
    if solvent_model := frame.model_chemistry.solvation_model:
        metadata["solvent"] = ImplicitSolvation(
            solvent_model=solvent_model,
            solvent=frame.model_chemistry.solvent,
        )
    return metadata


def _parse_running_time(file_content: str) -> Any | None:
    if matches := list(_RUN_TIME_RE.finditer(file_content)):
        matched = matches[-1]
        seconds = (
            int(matched.group("days")) * 86400
            + int(matched.group("hours")) * 3600
            + int(matched.group("minutes")) * 60
            + int(matched.group("seconds"))
            + int(matched.group("msec")) / 1000
        )
        return seconds * atom_ureg.second
    return None


def _parse_status(file_content: str) -> Status | None:
    if "****ORCA TERMINATED NORMALLY****" in file_content:
        return Status(normal_terminated=True, scf_converged=True)
    if "ORCA finished by error termination" in file_content or "ORCA TERMINATED ABNORMALLY" in file_content:
        return Status(normal_terminated=False, scf_converged=False)
    return None


def _split_orca_frames(file_content: str) -> list[str]:
    matches = list(_COORD_HEADER_RE.finditer(file_content))
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


def _last_frame_value(frames: Sequence[Any], field: str) -> Any:
    for frame in reversed(frames):
        value = getattr(frame, field, None)
        if value is not None:
            return value
    return None


def _first_frame_value(frames: Sequence[Any], field: str) -> Any:
    for frame in frames:
        value = getattr(frame, field, None)
        if value is not None:
            return value
    return None


class ORCALogFileParserMixin:
    @classmethod
    def _quick_check_file_format(cls, file_content: str) -> None:
        _ensure_orca_output(file_content)

    def _parse_metadata(self, file_content: str) -> dict[str, Any]:
        input_file_name, printed_input = _extract_printed_input(file_content)
        metadata: dict[str, Any] = {
            "qm_software": "ORCA",
            "input_file_name": input_file_name,
        }
        if version := _parse_version(file_content):
            metadata["qm_software_version"] = version
        metadata.update(_parse_printed_input_metadata(printed_input))
        if status := _parse_status(file_content):
            metadata["status"] = status
        if running_time := _parse_running_time(file_content):
            metadata["running_time"] = running_time
        return metadata

    def _split_file(self, file_content: str) -> Sequence[str]:
        return _split_orca_frames(file_content)

    def _update_file_metadata_from_frames(self, chem_file: Any, metadata: dict[str, Any]) -> None:
        base_parser = cast(_HasMetadataFinalize, super())
        base_parser._update_file_metadata_from_frames(chem_file, metadata)
        frames = list(getattr(chem_file, "frames", []))
        if not frames:
            return
        for field in ("solvent", "temperature", "pressure"):
            if getattr(chem_file, field, None) is not None:
                continue
            value = _first_frame_value(frames, field)
            if value is not None:
                setattr(chem_file, field, value)
        for field in (
            "status",
            "geometry_optimization_status",
            "electronic_states",
            "multireference_result",
        ):
            value = _last_frame_value(frames, field)
            if value is not None:
                setattr(chem_file, field, value)
        last_frame = frames[-1]
        if getattr(last_frame, "forces", None) is None:
            inherited_forces = _last_frame_value(frames[:-1], "forces")
            if inherited_forces is not None:
                last_frame.forces = inherited_forces

class ORCALogFileParserMemory(
    ORCALogFileParserMixin,
    BaseFileParserMemory[
        ORCALogFileMemory,
        ORCALogFileFrameMemory,
        ORCALogFileFrameParserMemory,
    ],
):
    _chem_file = ORCALogFileMemory
    _frame_parser = ORCALogFileFrameParserMemory


class ORCALogFileParserDisk(
    ORCALogFileParserMixin,
    BaseFileParserDisk[
        ORCALogFileDisk,
        ORCALogFileFrameDisk,
        ORCALogFileFrameParserDisk,
    ],
):
    allowed_formats = (".out", ".log", ".orcaout")
    _chem_file = ORCALogFileDisk
    _frame_parser = ORCALogFileFrameParserDisk


def register(registry: Registry) -> None:
    from typing import cast

    from molop.io.codecs._shared.reader_helpers import (
        ParserDiskReader,
        ReaderCodec,
        StructureLevel,
        extensions_for_parser,
    )

    extensions = frozenset(extensions_for_parser(ORCALogFileParserDisk))
    priority = 150

    @registry.reader_factory(format_id="orcaout", extensions=extensions, priority=priority)
    def _factory() -> ReaderCodec:
        return cast(
            ReaderCodec,
            ParserDiskReader(
                format_id="orcaout",
                extensions=extensions,
                level=StructureLevel.COORDS,
                parser_cls=ORCALogFileParserDisk,
                priority=priority,
            ),
        )
