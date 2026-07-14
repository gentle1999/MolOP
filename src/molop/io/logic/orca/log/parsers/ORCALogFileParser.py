from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from enum import Enum, auto
from typing import TYPE_CHECKING, Any, ClassVar, Protocol, cast

from molop.io.base_models.FileParser import BaseFileParserDisk, BaseFileParserMemory
from molop.io.base_models.ParseContainers import ModelParseResult, TextParseContext
from molop.io.base_models.source import LocatedSourceSegment
from molop.io.logic.orca.log.frame_models.ORCALogFileFrame import (
    ORCALogFileFrameDisk,
    ORCALogFileFrameMemory,
)
from molop.io.logic.orca.log.frame_parsers.ORCALogFileFrameParser import (
    ORCALogFileFrameParserDisk,
    ORCALogFileFrameParserMemory,
)
from molop.io.logic.orca.log.locators import locate_orca_job_frames, locate_orca_jobs
from molop.io.logic.orca.log.models.ORCALogFile import ORCALogFileDisk, ORCALogFileMemory
from molop.io.logic.orca.log.parsers._orca_log_file_extractors import (
    ensure_orca_output_content,
    extract_orca_output_version,
    extract_orca_printed_input,
    first_frame_value,
    last_frame_value,
    parse_orca_printed_input_metadata,
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
    format_id: ClassVar[str] = "orcaout"
    _assess_segment_calculation_status = True

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

    def _parse_segment_metadata_result(self, file_content: str) -> ModelParseResult:
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

    def _locate_segments(self, file_content: str) -> Sequence[LocatedSourceSegment]:
        return tuple(
            LocatedSourceSegment(
                segment=job,
                frames=locate_orca_job_frames(file_content, job),
            )
            for job in locate_orca_jobs(file_content)
        )

    def _parse_artifact_metadata(self, file_content: str) -> dict[str, Any] | None:
        metadata: dict[str, Any] = {"qm_software": "ORCA"}
        if version := extract_orca_output_version(file_content):
            metadata["qm_software_version"] = version
        return metadata

    def _parse_segment_metadata(
        self,
        segment_content: str,
        *,
        artifact_metadata: Mapping[str, Any],
    ) -> dict[str, Any]:
        metadata = self._parse_segment_metadata_result(segment_content).model_data()
        if not metadata.get("qm_software_version") and artifact_metadata.get("qm_software_version"):
            metadata["qm_software_version"] = artifact_metadata["qm_software_version"]
        return metadata

    def _prepare_file_metadata(
        self,
        artifact_metadata: Mapping[str, Any],
        segment_metadata: Sequence[Mapping[str, Any]],
    ) -> dict[str, Any]:
        metadata = dict(artifact_metadata)
        if segment_metadata:
            metadata.update(segment_metadata[0])
        running_times = [
            running_time
            for segment in segment_metadata
            if (running_time := segment.get("running_time")) is not None
        ]
        if running_times:
            metadata["running_time"] = sum(running_times[1:], running_times[0])
        for segment in reversed(segment_metadata):
            if (status := segment.get("status")) is not None:
                metadata["status"] = status
                break
        return metadata

    def _source_frame_role(
        self,
        frame: Any,
        *,
        frame_index: int,
        frame_count: int,
    ) -> str:
        _ = self
        task_types = {task.task_type for task in frame.task_requests if task.enabled}
        if "sp" in task_types and "opt" not in task_types and frame_count == 1:
            return "single_point"
        if frame_count == 1 or frame_index == frame_count - 1:
            return "terminal"
        if frame_index == 0:
            return "initial"
        return "intermediate"

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
        if metadata.get("status") is None:
            status = last_frame_value(frames, "status")
            if status is not None:
                chem_file.status = status
        for field in (
            "geometry_optimization_status",
            "electronic_states",
            "multireference_result",
        ):
            value = last_frame_value(frames, field)
            if value is not None:
                setattr(chem_file, field, value)


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

    @registry.reader_factory(
        format_id=ORCALogFileParserDisk.format_id,
        extensions=extensions,
        priority=priority,
    )
    def _factory() -> ReaderCodec:
        return cast(
            ReaderCodec,
            ParserDiskReader(
                format_id=ORCALogFileParserDisk.format_id,
                extensions=extensions,
                level=StructureLevel.COORDS,
                parser_cls=ORCALogFileParserDisk,
                priority=priority,
            ),
        )
