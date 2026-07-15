from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import TYPE_CHECKING, Any, ClassVar, cast

from molop.io.base_models.FileParser import BaseFileParserDisk, BaseFileParserMemory
from molop.io.base_models.source import LocatedSourceSegment
from molop.io.logic.xtb.output.frame_models.XTBOutputFileFrame import (
    XTBOutputFileFrameDisk,
    XTBOutputFileFrameMemory,
)
from molop.io.logic.xtb.output.frame_parsers.XTBOutputFileFrameParser import (
    XTBOutputFileFrameParserDisk,
    XTBOutputFileFrameParserMemory,
)
from molop.io.logic.xtb.output.locators import locate_xtb_runs
from molop.io.logic.xtb.output.models.XTBOutputFile import (
    XTBOutputFileDisk,
    XTBOutputFileMemory,
)
from molop.io.logic.xtb.output.parsers._xtb_output_extractors import (
    ensure_xtb_output_content,
    extract_xtb_metadata,
    extract_xtb_running_time,
    extract_xtb_status,
    extract_xtb_version,
)


if TYPE_CHECKING:
    from molop.io.codec_registry import Registry


class XTBOutputFileParserMixin:
    format_id: ClassVar[str] = "xtbout"
    _assess_segment_calculation_status = True

    @classmethod
    def _quick_check_file_format(cls, file_content: str) -> None:
        ensure_xtb_output_content(file_content)

    def _locate_segments(self, file_content: str) -> Sequence[LocatedSourceSegment]:
        return locate_xtb_runs(file_content)

    def _parse_artifact_metadata(self, file_content: str) -> dict[str, Any]:
        return {
            "qm_software": "xTB",
            "qm_software_version": extract_xtb_version(file_content) or "",
        }

    def _parse_segment_metadata(
        self,
        segment_content: str,
        *,
        artifact_metadata: Mapping[str, Any],
    ) -> dict[str, Any]:
        metadata = extract_xtb_metadata(segment_content)
        if not metadata.get("qm_software_version"):
            metadata["qm_software_version"] = artifact_metadata.get("qm_software_version", "")
        metadata["status"] = extract_xtb_status(segment_content)
        if running_time := extract_xtb_running_time(segment_content):
            metadata["running_time"] = running_time
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
            value
            for segment in segment_metadata
            if (value := segment.get("running_time")) is not None
        ]
        if running_times:
            metadata["running_time"] = sum(running_times[1:], running_times[0])
        if segment_metadata:
            last = segment_metadata[-1]
            for field in ("status", "qm_software_version"):
                if last.get(field) is not None:
                    metadata[field] = last[field]
        return metadata

    def _source_frame_role(
        self,
        frame: Any,
        *,
        frame_index: int,
        frame_count: int,
    ) -> str:
        _ = frame_index, frame_count
        task_types = {task.task_type for task in frame.task_requests if task.enabled}
        if task_types == {"sp"}:
            return "single_point"
        return "terminal"

    def _update_file_metadata_from_frames(self, chem_file: Any, metadata: dict[str, Any]) -> None:
        base_parser = cast(Any, super())
        base_parser._update_file_metadata_from_frames(chem_file, metadata)
        frames = list(getattr(chem_file, "frames", ()))
        if not frames:
            return
        first_frame = frames[0]
        last_frame = frames[-1]
        for field in ("solvent", "temperature", "electron_temperature"):
            if getattr(chem_file, field, None) is None:
                value = getattr(first_frame, field, None)
                if value is not None:
                    setattr(chem_file, field, value)
        for field in ("status", "geometry_optimization_status", "gradient_norm"):
            value = getattr(last_frame, field, None)
            if value is not None:
                setattr(chem_file, field, value)


class XTBOutputFileParserMemory(
    XTBOutputFileParserMixin,
    BaseFileParserMemory[
        XTBOutputFileMemory,
        XTBOutputFileFrameMemory,
        XTBOutputFileFrameParserMemory,
    ],
):
    _chem_file = XTBOutputFileMemory
    _frame_parser = XTBOutputFileFrameParserMemory


class XTBOutputFileParserDisk(
    XTBOutputFileParserMixin,
    BaseFileParserDisk[
        XTBOutputFileDisk,
        XTBOutputFileFrameDisk,
        XTBOutputFileFrameParserDisk,
    ],
):
    allowed_formats = (".out", ".log", ".xtbout")
    _chem_file = XTBOutputFileDisk
    _frame_parser = XTBOutputFileFrameParserDisk


def register(registry: Registry) -> None:
    from molop.io.codecs._shared.reader_helpers import (
        ParserDiskReader,
        ReaderCodec,
        StructureLevel,
        extensions_for_parser,
    )

    extensions = frozenset(extensions_for_parser(XTBOutputFileParserDisk))
    priority = 140

    @registry.reader_factory(
        format_id=XTBOutputFileParserDisk.format_id,
        extensions=extensions,
        priority=priority,
    )
    def _factory() -> ReaderCodec:
        return cast(
            ReaderCodec,
            ParserDiskReader(
                format_id=XTBOutputFileParserDisk.format_id,
                extensions=extensions,
                level=StructureLevel.COORDS,
                parser_cls=XTBOutputFileParserDisk,
                priority=priority,
            ),
        )


__all__ = ["XTBOutputFileParserDisk", "XTBOutputFileParserMemory"]
