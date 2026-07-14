"""
Author: TMJ
Date: 2025-07-30 14:30:03
LastEditors: TMJ
LastEditTime: 2026-03-23 19:01:43
Description: 请填写简介
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING, Any, ClassVar

from molop.io.base_models.FileParser import BaseFileParserDisk, BaseFileParserMemory
from molop.io.base_models.ParseContainers import ModelParseResult
from molop.io.base_models.source import LocatedSourceSegment
from molop.io.logic.gaussian.input.frame_models.GJFFileFrame import (
    GJFFileFrameDisk,
    GJFFileFrameMemory,
)
from molop.io.logic.gaussian.input.frame_parsers.GJFFileFrameParser import (
    GJFFileFrameParserDisk,
    GJFFileFrameParserMemory,
)
from molop.io.logic.gaussian.input.models.GJFFile import GJFFileDisk, GJFFileMemory
from molop.io.logic.gaussian.input.parsers._gjf_file_extractors import (
    ensure_gjf_content,
    locate_gjf_link1_frames,
)


if TYPE_CHECKING:
    from molop.io.codec_registry import Registry


class GJFFileParserMixin:
    format_id: ClassVar[str] = "gjf"

    @classmethod
    def _quick_check_file_format(cls, file_content: str) -> None:
        ensure_gjf_content(file_content)

    def _parse_artifact_metadata(self, file_content: str) -> dict[str, Any]:
        _ = file_content
        return ModelParseResult(
            {
                "qm_software": "Gaussian",
                "qm_software_version": "Any",
            }
        ).model_data()

    def _locate_segments(self, file_content: str) -> Sequence[LocatedSourceSegment]:
        frames = locate_gjf_link1_frames(file_content)
        return tuple(LocatedSourceSegment(segment=frame, frames=(frame,)) for frame in frames)


class GJFFileParserMemory(
    GJFFileParserMixin,
    BaseFileParserMemory[GJFFileMemory, GJFFileFrameMemory, GJFFileFrameParserMemory],
):
    _frame_parser = GJFFileFrameParserMemory
    _chem_file = GJFFileMemory


class GJFFileParserDisk(
    GJFFileParserMixin,
    BaseFileParserDisk[GJFFileDisk, GJFFileFrameDisk, GJFFileFrameParserDisk],
):
    allowed_formats = ("gjf", "gif", "com", ".gau", ".gjc")
    _frame_parser = GJFFileFrameParserDisk
    _chem_file = GJFFileDisk


def register(registry: Registry) -> None:
    """Register this file parser as a reader codec.

    Called by lazy activation via `molop.io.codecs.catalog`.
    """

    from typing import cast

    from molop.io.codecs._shared.reader_helpers import (
        ParserDiskReader,
        ReaderCodec,
        StructureLevel,
        extensions_for_parser,
    )

    extensions = frozenset(extensions_for_parser(GJFFileParserDisk))
    priority = 100

    @registry.reader_factory(
        format_id=GJFFileParserDisk.format_id,
        extensions=extensions,
        priority=priority,
    )
    def _factory() -> ReaderCodec:
        return cast(
            ReaderCodec,
            ParserDiskReader(
                format_id=GJFFileParserDisk.format_id,
                extensions=extensions,
                level=StructureLevel.COORDS,
                parser_cls=GJFFileParserDisk,
                priority=priority,
            ),
        )
