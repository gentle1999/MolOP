from __future__ import annotations

from dataclasses import dataclass
from enum import Enum, auto
from typing import TYPE_CHECKING, Any, Protocol, cast

from molop.io.base_models.FileParser import BaseFileParserDisk, BaseFileParserMemory
from molop.io.base_models.ParseContainers import ModelParseResult, TextParseContext
from molop.io.logic.orca.log.frame_models.ORCALogFileFrame import (
    ORCALogFileFrameDisk,
    ORCALogFileFrameMemory,
)
from molop.io.logic.orca.log.frame_parsers.ORCALogFileFrameParser import (
    ORCALogFileFrameParserDisk,
    ORCALogFileFrameParserMemory,
)
from molop.io.logic.orca.log.models.ORCALogFile import ORCALogFileDisk, ORCALogFileMemory
from molop.io.logic.orca.log.parsers._orca_log_file_extractors import (
    ensure_orca_output_content,
    extract_orca_output_version,
    extract_orca_printed_input,
    first_frame_value,
    last_frame_value,
    parse_orca_printed_input_metadata,
    split_orca_output_frames,
)
from molop.io.logic.orca.log.parsers._orca_log_shared import (
    extract_orca_running_time,
    extract_orca_status,
)


if TYPE_CHECKING:
    from molop.io.codec_registry import Registry


class _HasMetadataFinalize(Protocol):
    def _update_file_metadata_from_frames(
        self, chem_file: Any, metadata: dict[str, Any]
    ) -> None: ...


class ORCALogMetadataParsePhase(Enum):
    """Explicit stages for ORCA output file metadata parsing."""

    SOFTWARE = auto()
    PRINTED_INPUT = auto()
    STATUS = auto()
    DONE = auto()


@dataclass(slots=True)
class ORCALogMetadataParseContext:
    """Mutable metadata context for ORCA output file parser phases."""

    text: TextParseContext
    printed_input: str = ""


class ORCALogFileParserMixin:
    @classmethod
    def _quick_check_file_format(cls, file_content: str) -> None:
        ensure_orca_output_content(file_content)

    def _run_software_metadata_phase(
        self, context: ORCALogMetadataParseContext, result: ModelParseResult
    ) -> ORCALogMetadataParsePhase:
        input_file_name, printed_input = extract_orca_printed_input(context.text.content)
        result.set("input_file_name", input_file_name)
        context.printed_input = printed_input
        if version := extract_orca_output_version(context.text.content):
            result.set("qm_software_version", version)
        return ORCALogMetadataParsePhase.PRINTED_INPUT

    def _run_printed_input_metadata_phase(
        self, context: ORCALogMetadataParseContext, result: ModelParseResult
    ) -> ORCALogMetadataParsePhase:
        result.update(parse_orca_printed_input_metadata(context.printed_input))
        return ORCALogMetadataParsePhase.STATUS

    def _run_status_metadata_phase(
        self, context: ORCALogMetadataParseContext, result: ModelParseResult
    ) -> ORCALogMetadataParsePhase:
        if status := extract_orca_status(context.text.content):
            result.set("status", status)
        if running_time := extract_orca_running_time(context.text.content):
            result.set("running_time", running_time)
        return ORCALogMetadataParsePhase.DONE

    def _parse_metadata_result(self, file_content: str) -> ModelParseResult:
        context = ORCALogMetadataParseContext(TextParseContext(file_content))
        result = ModelParseResult({"qm_software": "ORCA"})
        phase = ORCALogMetadataParsePhase.SOFTWARE
        while phase is not ORCALogMetadataParsePhase.DONE:
            if phase is ORCALogMetadataParsePhase.SOFTWARE:
                phase = self._run_software_metadata_phase(context, result)
            elif phase is ORCALogMetadataParsePhase.PRINTED_INPUT:
                phase = self._run_printed_input_metadata_phase(context, result)
            elif phase is ORCALogMetadataParsePhase.STATUS:
                phase = self._run_status_metadata_phase(context, result)
            else:
                raise AssertionError(f"Unexpected ORCA log metadata parse phase: {phase!r}")
        return result

    def _parse_metadata(self, file_content: str) -> dict[str, Any]:
        return self._parse_metadata_result(file_content).model_data()

    def _split_file(self, file_content: str) -> list[str]:
        return split_orca_output_frames(file_content)

    def _update_file_metadata_from_frames(self, chem_file: Any, metadata: dict[str, Any]) -> None:
        base_parser = cast(_HasMetadataFinalize, super())
        base_parser._update_file_metadata_from_frames(chem_file, metadata)
        frames = list(getattr(chem_file, "frames", []))
        if not frames:
            return
        for field in ("solvent", "temperature", "pressure"):
            if getattr(chem_file, field, None) is not None:
                continue
            value = first_frame_value(frames, field)
            if value is not None:
                setattr(chem_file, field, value)
        for field in (
            "status",
            "geometry_optimization_status",
            "electronic_states",
            "multireference_result",
        ):
            value = last_frame_value(frames, field)
            if value is not None:
                setattr(chem_file, field, value)
        last_frame = frames[-1]
        if getattr(last_frame, "forces", None) is None:
            inherited_forces = last_frame_value(frames[:-1], "forces")
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
