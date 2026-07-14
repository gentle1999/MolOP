"""
Author: TMJ
Date: 2025-07-30 10:30:16
LastEditors: TMJ
LastEditTime: 2026-02-05 19:52:55
Description: 请填写简介
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING, ClassVar

from molop.io.base_models.FileParser import BaseFileParserDisk, BaseFileParserMemory
from molop.io.base_models.source import LocatedSourceSegment, LocatedTextBlock
from molop.io.logic.coords.frame_models.SDFFileFrame import SDFFileFrameDisk, SDFFileFrameMemory
from molop.io.logic.coords.frame_parsers.SDFFileFrameParser import (
    SDFFileFrameParserDisk,
    SDFFileFrameParserMemory,
)
from molop.io.logic.coords.models.SDFFile import SDFFileDisk, SDFFileMemory
from molop.io.logic.coords.parsers._coords_file_extractors import locate_sdf_frames


if TYPE_CHECKING:
    from molop.io.codec_registry import Registry


class SDFFileParserMixin:
    format_id: ClassVar[str] = "sdf"

    @classmethod
    def _quick_check_file_format(cls, file_content: str) -> None:
        _ = file_content

    def _locate_segments(self, file_content: str) -> Sequence[LocatedSourceSegment]:
        frames = locate_sdf_frames(file_content)
        return (
            LocatedSourceSegment(
                segment=LocatedTextBlock(0, len(file_content)),
                frames=frames,
            ),
        )


class SDFFileParserMemory(
    SDFFileParserMixin,
    BaseFileParserMemory[SDFFileMemory, SDFFileFrameMemory, SDFFileFrameParserMemory],
):
    _frame_parser = SDFFileFrameParserMemory
    _chem_file = SDFFileMemory


class SDFFileParserDisk(
    SDFFileParserMixin,
    BaseFileParserDisk[SDFFileDisk, SDFFileFrameDisk, SDFFileFrameParserDisk],
):
    allowed_formats = ("sdf", "sd", "mol")
    _frame_parser = SDFFileFrameParserDisk
    _chem_file = SDFFileDisk


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

    extensions = frozenset(extensions_for_parser(SDFFileParserDisk))
    priority = 100

    @registry.reader_factory(
        format_id=SDFFileParserDisk.format_id,
        extensions=extensions,
        priority=priority,
    )
    def _factory() -> ReaderCodec:
        return cast(
            ReaderCodec,
            ParserDiskReader(
                format_id=SDFFileParserDisk.format_id,
                extensions=extensions,
                level=StructureLevel.COORDS,
                parser_cls=SDFFileParserDisk,
                priority=priority,
            ),
        )
