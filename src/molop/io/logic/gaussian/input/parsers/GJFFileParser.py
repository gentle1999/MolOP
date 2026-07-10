"""
Author: TMJ
Date: 2025-07-30 14:30:03
LastEditors: TMJ
LastEditTime: 2026-03-23 19:01:43
Description: 请填写简介
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field
from enum import Enum, auto
from typing import TYPE_CHECKING, Any

from molop.io.base_models.FileParser import BaseFileParserDisk, BaseFileParserMemory
from molop.io.base_models.ParseContainers import ModelParseResult
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
    expand_gjf_at_includes,
    find_gjf_link1_matches,
    split_gjf_link1_frames,
    validate_gjf_link1_boundaries,
)


if TYPE_CHECKING:
    from molop.io.codec_registry import Registry


class GJFFileSplitPhase(Enum):
    """Explicit stages for Gaussian input file-level splitting."""

    EXPAND_INCLUDES = auto()
    VALIDATE_LINK1 = auto()
    COLLECT_FRAMES = auto()
    DONE = auto()


@dataclass(slots=True)
class GJFFileSplitContext:
    """Mutable file-splitting context for Gaussian input files."""

    file_content: str
    normalized_content: str = ""
    link1_matches: list[Any] = field(default_factory=list)
    frames: list[str] = field(default_factory=list)


class GJFFileParserMixin:
    @classmethod
    def _quick_check_file_format(cls, file_content: str) -> None:
        ensure_gjf_content(file_content)

    def _parse_metadata_result(self, file_content: str) -> ModelParseResult:
        _ = file_content
        return ModelParseResult(
            {
                "qm_software": "Gaussian",
                "qm_software_version": "Any",
            }
        )

    def _parse_metadata(self, file_content: str) -> dict[str, Any]:
        return self._parse_metadata_result(file_content).model_data()

    def _run_expand_includes_split_phase(self, context: GJFFileSplitContext) -> GJFFileSplitPhase:
        context.normalized_content = expand_gjf_at_includes(
            context.file_content,
            getattr(self, "_file_path", None),
        )
        context.link1_matches = find_gjf_link1_matches(context.normalized_content)
        return GJFFileSplitPhase.VALIDATE_LINK1

    def _run_validate_link1_split_phase(self, context: GJFFileSplitContext) -> GJFFileSplitPhase:
        validate_gjf_link1_boundaries(context.normalized_content, context.link1_matches)
        return GJFFileSplitPhase.COLLECT_FRAMES

    def _run_collect_frames_split_phase(self, context: GJFFileSplitContext) -> GJFFileSplitPhase:
        context.frames = split_gjf_link1_frames(context.normalized_content, context.link1_matches)
        return GJFFileSplitPhase.DONE

    def _split_file(self, file_content: str) -> Sequence[str]:
        context = GJFFileSplitContext(file_content)
        phase = GJFFileSplitPhase.EXPAND_INCLUDES
        while phase is not GJFFileSplitPhase.DONE:
            if phase is GJFFileSplitPhase.EXPAND_INCLUDES:
                phase = self._run_expand_includes_split_phase(context)
            elif phase is GJFFileSplitPhase.VALIDATE_LINK1:
                phase = self._run_validate_link1_split_phase(context)
            elif phase is GJFFileSplitPhase.COLLECT_FRAMES:
                phase = self._run_collect_frames_split_phase(context)
            else:
                raise AssertionError(f"Unexpected GJF file split phase: {phase!r}")
        return context.frames


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
