"""
Author: TMJ
Date: 2025-07-30 10:30:16
LastEditors: TMJ
LastEditTime: 2026-02-05 19:52:55
Description: 请填写简介
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING, Any

from molop.io.base_models.FileParser import BaseFileParserDisk, BaseFileParserMemory
from molop.io.base_models.ParseContainers import ModelParseResult
from molop.io.logic.coords.frame_models.SDFFileFrame import SDFFileFrameDisk, SDFFileFrameMemory
from molop.io.logic.coords.frame_parsers.SDFFileFrameParser import (
    SDFFileFrameParserDisk,
    SDFFileFrameParserMemory,
)
from molop.io.logic.coords.models.SDFFile import SDFFileDisk, SDFFileMemory
from molop.io.logic.coords.parsers._coords_file_extractors import split_sdf_frames


if TYPE_CHECKING:
    from molop.io.codec_registry import Registry


class SDFFileParserMixin:
    @classmethod
    def _quick_check_file_format(cls, file_content: str) -> None:
        _ = file_content

    def _parse_metadata_result(self, file_content: str) -> ModelParseResult:
        _ = file_content
        return ModelParseResult()

    def _parse_metadata(self, file_content: str) -> dict[str, Any]:
        return self._parse_metadata_result(file_content).model_data()

    def _split_file(self, file_content: str) -> Sequence[str]:
        return split_sdf_frames(file_content)


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

    @registry.reader_factory(format_id="sdf", extensions=extensions, priority=priority)
    def _factory() -> ReaderCodec:
        return cast(
            ReaderCodec,
            ParserDiskReader(
                format_id="sdf",
                extensions=extensions,
                level=StructureLevel.COORDS,
                parser_cls=SDFFileParserDisk,
                priority=priority,
            ),
        )
