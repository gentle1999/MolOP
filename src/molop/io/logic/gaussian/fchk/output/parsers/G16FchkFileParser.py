from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import TYPE_CHECKING, Any, ClassVar, cast

from molop.io.base_models.FileParser import BaseFileParserDisk, BaseFileParserMemory
from molop.io.base_models.source import LocatedSourceSegment
from molop.io.logic.gaussian.fchk.output.frame_models.G16FchkFileFrame import (
    G16FchkFileFrameDisk,
    G16FchkFileFrameMemory,
)
from molop.io.logic.gaussian.fchk.output.frame_parsers.G16FchkFileFrameParser import (
    G16FchkFileFrameParserDisk,
    G16FchkFileFrameParserMemory,
)
from molop.io.logic.gaussian.fchk.output.locators import locate_fchk_content
from molop.io.logic.gaussian.fchk.output.models.G16FchkFile import (
    G16FchkFileDisk,
    G16FchkFileMemory,
)
from molop.io.logic.gaussian.fchk.output.parsers._fchk_extractors import (
    extract_fchk_metadata,
)
from molop.io.logic.gaussian.fchk.output.parsers._fchk_records import (
    ensure_fchk_content,
)


if TYPE_CHECKING:
    from molop.io.codec_registry import Registry


class G16FchkFileParserMixin:
    format_id: ClassVar[str] = "g16fchk"
    _assess_segment_calculation_status = True

    @classmethod
    def _quick_check_file_format(cls, file_content: str) -> None:
        ensure_fchk_content(file_content)

    def _locate_segments(self, file_content: str) -> Sequence[LocatedSourceSegment]:
        return locate_fchk_content(file_content)

    def _parse_artifact_metadata(self, file_content: str) -> dict[str, Any]:
        _ = file_content
        return {"qm_software": "Gaussian"}

    def _parse_segment_metadata(
        self,
        segment_content: str,
        *,
        artifact_metadata: Mapping[str, Any],
    ) -> dict[str, Any]:
        _ = artifact_metadata
        return extract_fchk_metadata(segment_content)

    def _prepare_file_metadata(
        self,
        artifact_metadata: Mapping[str, Any],
        segment_metadata: Sequence[Mapping[str, Any]],
    ) -> dict[str, Any]:
        metadata = dict(artifact_metadata)
        if segment_metadata:
            metadata.update(segment_metadata[0])
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
        return "single_point" if task_types == {"sp"} else "terminal"

    def _update_file_metadata_from_frames(self, chem_file: Any, metadata: dict[str, Any]) -> None:
        base_parser = cast(Any, super())
        base_parser._update_file_metadata_from_frames(chem_file, metadata)
        frames = list(getattr(chem_file, "frames", ()))
        if not frames:
            return
        last_frame = frames[-1]
        for field in ("status", "geometry_optimization_status"):
            value = getattr(last_frame, field, None)
            if value is not None:
                setattr(chem_file, field, value)


class G16FchkFileParserMemory(
    G16FchkFileParserMixin,
    BaseFileParserMemory[
        G16FchkFileMemory,
        G16FchkFileFrameMemory,
        G16FchkFileFrameParserMemory,
    ],
):
    _chem_file = G16FchkFileMemory
    _frame_parser = G16FchkFileFrameParserMemory


class G16FchkFileParserDisk(
    G16FchkFileParserMixin,
    BaseFileParserDisk[
        G16FchkFileDisk,
        G16FchkFileFrameDisk,
        G16FchkFileFrameParserDisk,
    ],
):
    allowed_formats = (".fchk", ".fch", ".fck")
    _chem_file = G16FchkFileDisk
    _frame_parser = G16FchkFileFrameParserDisk


def register(registry: Registry) -> None:
    from molop.io.codecs._shared.reader_helpers import (
        ParserDiskReader,
        ReaderCodec,
        StructureLevel,
        extensions_for_parser,
    )

    extensions = frozenset(extensions_for_parser(G16FchkFileParserDisk))
    priority = 180

    @registry.reader_factory(
        format_id=G16FchkFileParserDisk.format_id,
        extensions=extensions,
        priority=priority,
    )
    def _factory() -> ReaderCodec:
        return cast(
            ReaderCodec,
            ParserDiskReader(
                format_id=G16FchkFileParserDisk.format_id,
                extensions=extensions,
                level=StructureLevel.COORDS,
                parser_cls=G16FchkFileParserDisk,
                priority=priority,
                memory_parser_cls=G16FchkFileParserMemory,
            ),
        )


__all__ = ["G16FchkFileParserDisk", "G16FchkFileParserMemory"]
