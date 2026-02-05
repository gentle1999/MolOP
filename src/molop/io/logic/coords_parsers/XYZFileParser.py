"""
Author: TMJ
Date: 2025-07-29 22:53:30
LastEditors: TMJ
LastEditTime: 2026-02-05 19:53:43
Description: 请填写简介
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING, Any

from molop.io.base_models.FileParser import BaseFileParserDisk, BaseFileParserMemory
from molop.io.logic.coords_frame_models.XYZFileFrame import XYZFileFrameDisk, XYZFileFrameMemory
from molop.io.logic.coords_frame_parsers.XYZFileFrameParser import (
    XYZFileFrameParserDisk,
    XYZFileFrameParserMemory,
)
from molop.io.logic.coords_models.XYZFile import XYZFileDisk, XYZFileMemory


if TYPE_CHECKING:
    from molop.io.codec_registry import Registry


class XYZFileParserMixin:
    def _parse_metadata(self, file_content: str) -> dict[str, Any] | None: ...

    def _split_file(self, file_content: str) -> Sequence[str]:
        lines = file_content.splitlines()
        anchor = 0
        frames = []
        while anchor < len(lines):
            try:
                num_atoms = int(lines[anchor])
            except ValueError:
                anchor += 1
                continue
            frames.append("\n".join(lines[anchor : anchor + num_atoms + 2]))
            anchor += num_atoms + 2
        return frames


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

    @registry.reader_factory(format_id="xyz", extensions=extensions, priority=priority)
    def _factory() -> ReaderCodec:
        return cast(
            ReaderCodec,
            ParserDiskReader(
                format_id="xyz",
                extensions=extensions,
                level=StructureLevel.COORDS,
                parser_cls=XYZFileParserDisk,
                priority=priority,
            ),
        )
