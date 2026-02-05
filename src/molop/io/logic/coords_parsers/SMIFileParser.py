"""
Author: TMJ
Date: 2025-12-14 23:30:46
LastEditors: TMJ
LastEditTime: 2025-12-14 23:41:30
Description: 请填写简介
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING, Any

from molop.io.base_models.FileParser import BaseFileParserDisk, BaseFileParserMemory
from molop.io.logic.coords_frame_models.SMIFileFrame import SMIFileFrameDisk, SMIFileFrameMemory
from molop.io.logic.coords_frame_parsers.SMIFileFrameParser import (
    SMIFileFrameParserDisk,
    SMIFileFrameParserMemory,
)
from molop.io.logic.coords_models.SMIFile import SMIFileDisk, SMIFileMemory


if TYPE_CHECKING:
    from molop.io.codec_registry import Registry


class SMIFileParserMixin:
    def _parse_metadata(self, file_content: str) -> dict[str, Any] | None: ...

    def _split_file(self, file_content: str) -> Sequence[str]:
        return file_content.splitlines()


class SMIFileParserMemory(
    SMIFileParserMixin,
    BaseFileParserMemory[SMIFileMemory, SMIFileFrameMemory, SMIFileFrameParserMemory],
):
    _frame_parser = SMIFileFrameParserMemory
    _chem_file = SMIFileMemory


class SMIFileParserDisk(
    SMIFileParserMixin,
    BaseFileParserDisk[SMIFileDisk, SMIFileFrameDisk, SMIFileFrameParserDisk],
):
    allowed_formats = ("smi", "txt")
    _frame_parser = SMIFileFrameParserDisk
    _chem_file = SMIFileDisk


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

    extensions = frozenset(extensions_for_parser(SMIFileParserDisk))
    priority = 100

    @registry.reader_factory(format_id="smi", extensions=extensions, priority=priority)
    def _factory() -> ReaderCodec:
        return cast(
            ReaderCodec,
            ParserDiskReader(
                format_id="smi",
                extensions=extensions,
                level=StructureLevel.COORDS,
                parser_cls=SMIFileParserDisk,
                priority=priority,
            ),
        )
