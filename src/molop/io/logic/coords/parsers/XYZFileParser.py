"""
Author: TMJ
Date: 2025-07-29 22:53:30
LastEditors: TMJ
LastEditTime: 2026-02-05 19:53:43
Description: 请填写简介
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING, ClassVar

from molop.io.base_models.FileParser import BaseFileParserDisk, BaseFileParserMemory
from molop.io.base_models.source import LocatedSourceSegment, LocatedTextBlock
from molop.io.logic.coords.frame_models.XYZFileFrame import XYZFileFrameDisk, XYZFileFrameMemory
from molop.io.logic.coords.frame_parsers.XYZFileFrameParser import (
    XYZFileFrameParserDisk,
    XYZFileFrameParserMemory,
)
from molop.io.logic.coords.models.XYZFile import XYZFileDisk, XYZFileMemory
from molop.io.logic.coords.parsers._coords_file_extractors import locate_xyz_frames


if TYPE_CHECKING:
    from molop.io.codec_registry import Registry


class XYZFileParserMixin:
    format_id: ClassVar[str] = "xyz"

    @classmethod
    def _quick_check_file_format(cls, file_content: str) -> None:
        _ = file_content

    def _locate_segments(self, file_content: str) -> Sequence[LocatedSourceSegment]:
        frames = locate_xyz_frames(file_content)
        return (
            LocatedSourceSegment(
                segment=LocatedTextBlock(0, len(file_content)),
                frames=frames,
            ),
        )


class XYZFileParserMemory(
    XYZFileParserMixin,
    BaseFileParserMemory[XYZFileMemory, XYZFileFrameMemory, XYZFileFrameParserMemory],
):
    _frame_parser = XYZFileFrameParserMemory
    _chem_file = XYZFileMemory


class XYZFileParserDisk(
    XYZFileParserMixin,
    BaseFileParserDisk[XYZFileDisk, XYZFileFrameDisk, XYZFileFrameParserDisk],
):
    allowed_formats = ("xyz",)
    _frame_parser = XYZFileFrameParserDisk
    _chem_file = XYZFileDisk


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

    extensions = frozenset(extensions_for_parser(XYZFileParserDisk))
    priority = 100

    @registry.reader_factory(
        format_id=XYZFileParserDisk.format_id,
        extensions=extensions,
        priority=priority,
    )
    def _factory() -> ReaderCodec:
        return cast(
            ReaderCodec,
            ParserDiskReader(
                format_id=XYZFileParserDisk.format_id,
                extensions=extensions,
                level=StructureLevel.COORDS,
                parser_cls=XYZFileParserDisk,
                priority=priority,
            ),
        )
