"""
Author: TMJ
Date: 2026-02-10 00:00:00
LastEditors: TMJ
LastEditTime: 2026-02-10 00:00:00
Description: ORCA input file parsers
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING, Any, ClassVar

from molop.io.base_models.FileParser import BaseFileParserDisk, BaseFileParserMemory
from molop.io.base_models.ParseContainers import ModelParseResult
from molop.io.base_models.source import LocatedSourceSegment
from molop.io.logic.orca.input.frame_models.ORCAInpFileFrame import (
    ORCAInpFileFrameDisk,
    ORCAInpFileFrameMemory,
)
from molop.io.logic.orca.input.frame_parsers.ORCAInpFileFrameParser import (
    ORCAInpFileFrameParserDisk,
    ORCAInpFileFrameParserMemory,
)
from molop.io.logic.orca.input.models.ORCAInpFile import ORCAInpFileDisk, ORCAInpFileMemory
from molop.io.logic.orca.input.parsers._orca_inp_file_extractors import (
    ensure_orca_input_content,
    locate_orca_input_frames,
)


if TYPE_CHECKING:
    from molop.io.codec_registry import Registry


class ORCAInpFileParserMixin:
    format_id: ClassVar[str] = "orcainp"

    @classmethod
    def _quick_check_file_format(cls, file_content: str) -> None:
        ensure_orca_input_content(file_content)

    def _parse_artifact_metadata(self, file_content: str) -> dict[str, Any]:
        _ = file_content
        return ModelParseResult(
            {
                "qm_software": "ORCA",
                "qm_software_version": "Any",
            }
        ).model_data()

    def _locate_segments(self, file_content: str) -> Sequence[LocatedSourceSegment]:
        frames = locate_orca_input_frames(file_content)
        return tuple(LocatedSourceSegment(segment=frame, frames=(frame,)) for frame in frames)


class ORCAInpFileParserMemory(
    ORCAInpFileParserMixin,
    BaseFileParserMemory[
        ORCAInpFileMemory,
        ORCAInpFileFrameMemory,
        ORCAInpFileFrameParserMemory,
    ],
):
    _frame_parser = ORCAInpFileFrameParserMemory
    _chem_file = ORCAInpFileMemory


class ORCAInpFileParserDisk(
    ORCAInpFileParserMixin,
    BaseFileParserDisk[
        ORCAInpFileDisk,
        ORCAInpFileFrameDisk,
        ORCAInpFileFrameParserDisk,
    ],
):
    allowed_formats = ("inp",)
    _frame_parser = ORCAInpFileFrameParserDisk
    _chem_file = ORCAInpFileDisk


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

    extensions = frozenset(extensions_for_parser(ORCAInpFileParserDisk))
    priority = 100

    @registry.reader_factory(
        format_id=ORCAInpFileParserDisk.format_id,
        extensions=extensions,
        priority=priority,
    )
    def _factory() -> ReaderCodec:
        return cast(
            ReaderCodec,
            ParserDiskReader(
                format_id=ORCAInpFileParserDisk.format_id,
                extensions=extensions,
                level=StructureLevel.COORDS,
                parser_cls=ORCAInpFileParserDisk,
                priority=priority,
                memory_parser_cls=ORCAInpFileParserMemory,
            ),
        )
