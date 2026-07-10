"""
Author: TMJ
Date: 2026-02-10 00:00:00
LastEditors: TMJ
LastEditTime: 2026-02-10 00:00:00
Description: ORCA input file parsers
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field
from enum import Enum, auto
from typing import TYPE_CHECKING, Any

from molop.io.base_models.FileParser import BaseFileParserDisk, BaseFileParserMemory
from molop.io.base_models.ParseContainers import ModelParseResult
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
    build_orca_input_file_lines,
    ensure_orca_input_content,
    split_orca_input_frames,
)


if TYPE_CHECKING:
    from molop.io.codec_registry import Registry


class ORCAInpFileSplitPhase(Enum):
    """Explicit stages for ORCA input file-level splitting."""

    PREPARE_LINES = auto()
    SCAN_LINES = auto()
    DONE = auto()


@dataclass(slots=True)
class ORCAInpFileSplitContext:
    """Mutable file-splitting context for ORCA input files."""

    file_content: str
    lines: list[str] = field(default_factory=list)
    frames: list[str] = field(default_factory=list)


class ORCAInpFileParserMixin:
    @classmethod
    def _quick_check_file_format(cls, file_content: str) -> None:
        ensure_orca_input_content(file_content)

    def _parse_metadata_result(self, file_content: str) -> ModelParseResult:
        _ = file_content
        return ModelParseResult(
            {
                "qm_software": "ORCA",
                "qm_software_version": "Any",
            }
        )

    def _parse_metadata(self, file_content: str) -> dict[str, Any]:
        return self._parse_metadata_result(file_content).model_data()

    def _run_prepare_lines_split_phase(
        self, context: ORCAInpFileSplitContext
    ) -> ORCAInpFileSplitPhase:
        context.lines = build_orca_input_file_lines(context.file_content)
        return ORCAInpFileSplitPhase.SCAN_LINES

    def _run_scan_lines_split_phase(
        self, context: ORCAInpFileSplitContext
    ) -> ORCAInpFileSplitPhase:
        context.frames = split_orca_input_frames(context.lines)
        return ORCAInpFileSplitPhase.DONE

    def _split_file(self, file_content: str) -> Sequence[str]:
        context = ORCAInpFileSplitContext(file_content)
        phase = ORCAInpFileSplitPhase.PREPARE_LINES
        while phase is not ORCAInpFileSplitPhase.DONE:
            if phase is ORCAInpFileSplitPhase.PREPARE_LINES:
                phase = self._run_prepare_lines_split_phase(context)
            elif phase is ORCAInpFileSplitPhase.SCAN_LINES:
                phase = self._run_scan_lines_split_phase(context)
            else:
                raise AssertionError(f"Unexpected ORCA input file split phase: {phase!r}")
        return context.frames


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

    @registry.reader_factory(format_id="orcainp", extensions=extensions, priority=priority)
    def _factory() -> ReaderCodec:
        return cast(
            ReaderCodec,
            ParserDiskReader(
                format_id="orcainp",
                extensions=extensions,
                level=StructureLevel.COORDS,
                parser_cls=ORCAInpFileParserDisk,
                priority=priority,
            ),
        )
