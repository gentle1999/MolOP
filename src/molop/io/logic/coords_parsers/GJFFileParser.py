"""
Author: TMJ
Date: 2025-07-30 14:30:03
LastEditors: TMJ
LastEditTime: 2026-02-05 19:23:34
Description: 请填写简介
"""

from __future__ import annotations

import re
from collections.abc import Sequence
from typing import TYPE_CHECKING

from molop.io.base_models.FileParser import BaseFileParserDisk, BaseFileParserMemory
from molop.io.logic.coords_frame_models.GJFFileFrame import GJFFileFrameDisk, GJFFileFrameMemory
from molop.io.logic.coords_frame_parsers.GJFFileFrameParser import (
    GJFFileFrameParserDisk,
    GJFFileFrameParserMemory,
)
from molop.io.logic.coords_models.GJFFile import GJFFileDisk, GJFFileMemory


if TYPE_CHECKING:
    from molop.io.codec_registry import Registry


class GJFFileParserMixin:
    def _parse_metadata(self, file_content: str): ...

    def _split_file(self, file_content: str) -> Sequence[str]:
        return [
            frame.strip() + "\n\n" for frame in re.split(r"--[lL][iI][nN][kK]1--", file_content)
        ]


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

    @registry.reader_factory(format_id="gjf", extensions=extensions, priority=priority)
    def _factory() -> ReaderCodec:
        return cast(
            ReaderCodec,
            ParserDiskReader(
                format_id="gjf",
                extensions=extensions,
                level=StructureLevel.COORDS,
                parser_cls=GJFFileParserDisk,
                priority=priority,
            ),
        )
