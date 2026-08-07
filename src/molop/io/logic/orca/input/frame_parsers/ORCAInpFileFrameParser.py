"""
Author: TMJ
Date: 2026-02-10 00:00:00
LastEditors: TMJ
LastEditTime: 2026-02-10 00:00:00
Description: ORCA input frame parsers
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from enum import Enum, auto
from typing import Any

from molop.io.base_models.FrameParser import BaseFrameParser, FrameParseContext
from molop.io.base_models.ParseContainers import ModelParseResult
from molop.io.logic.orca.common import (
    ORCABlock,
    ORCACommentLine,
    ORCAGeometry,
    ORCAKeywordLine,
)
from molop.io.logic.orca.input.frame_models.ORCAInpFileFrame import (
    ORCAInpFileFrameDisk,
    ORCAInpFileFrameMemory,
)
from molop.io.logic.orca.input.frame_parsers._orca_inp_blocks import (
    ORCABlockSpan,
    extract_block_spans,
    line_in_spans,
)
from molop.io.logic.orca.input.frame_parsers._orca_inp_geometry import (
    ORCAGeometrySection,
    extract_percent_coords_geometry,
    extract_star_geometry,
    geometry_from_section,
)
from molop.io.logic.orca.input.frame_parsers._orca_inp_resources import (
    parse_output_print_settings,
    parse_request_memory,
    parse_request_num_cpu,
    render_resource_blocks,
)
from molop.io.logic.orca.input.frame_parsers._orca_inp_semantics import (
    build_orca_input_semantic_payload,
)


class ORCAInpParsePhase(Enum):
    """Explicit stages for ORCA input parsing."""

    STRUCTURE = auto()
    LINES = auto()
    GEOMETRY = auto()
    RESOURCES = auto()
    SEMANTICS = auto()
    DONE = auto()


@dataclass(slots=True)
class ORCAInpParseContext:
    """Mutable parse context shared by ORCA input parser phases."""

    block: str
    lines: list[str]
    geometry_section: ORCAGeometrySection | None = None
    block_spans: list[ORCABlockSpan] = field(default_factory=list)
    blocks: list[ORCABlock] = field(default_factory=list)
    comments: list[ORCACommentLine] = field(default_factory=list)
    keywords: list[ORCAKeywordLine] = field(default_factory=list)
    trailing_lines: list[str] = field(default_factory=list)
    geometry: ORCAGeometry | None = None
    keyword_text: str = ""


class ORCAInpFileFrameParserMixin:
    def _run_structure_phase(
        self, context: ORCAInpParseContext, result: ModelParseResult
    ) -> ORCAInpParsePhase:
        _ = result
        context.geometry_section = extract_star_geometry(context.block)
        if context.geometry_section is None:
            context.geometry_section = extract_percent_coords_geometry(context.block)

        context.block_spans = extract_block_spans(context.block)
        context.blocks = [span.block for span in context.block_spans]
        return ORCAInpParsePhase.LINES

    def _run_line_phase(
        self, context: ORCAInpParseContext, result: ModelParseResult
    ) -> ORCAInpParsePhase:
        _ = result
        occupied_spans: list[tuple[int, int]] = [
            (span.line_start, span.line_end) for span in context.block_spans
        ]
        if context.geometry_section is not None:
            occupied_spans.append(
                (context.geometry_section.line_start, context.geometry_section.line_end)
            )

        for line_idx, line in enumerate(context.lines):
            if line_in_spans(line_idx, occupied_spans):
                continue
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("#"):
                context.comments.append(ORCACommentLine(text=stripped[1:].strip()))
            elif stripped.startswith("!"):
                context.keywords.append(ORCAKeywordLine(text=stripped[1:].strip()))
            elif not stripped.lower().startswith("$new_job"):
                context.trailing_lines.append(line)
        return ORCAInpParsePhase.GEOMETRY

    def _run_geometry_phase(
        self, context: ORCAInpParseContext, result: ModelParseResult
    ) -> ORCAInpParsePhase:
        context.geometry = geometry_from_section(context.geometry_section, context.blocks)
        context.keyword_text = "\n".join(
            line.text for line in context.keywords if line.text.strip()
        )

        result.set("comment_lines", context.comments)
        result.set("keyword_lines", context.keywords)
        result.set("blocks", context.blocks)
        result.set("geometry", context.geometry)
        result.set("trailing_lines", context.trailing_lines)
        return ORCAInpParsePhase.RESOURCES

    def _run_resource_phase(
        self, context: ORCAInpParseContext, result: ModelParseResult
    ) -> ORCAInpParsePhase:
        result.set("resources_raw", render_resource_blocks(context.blocks))
        result.set("request_num_cpu", parse_request_num_cpu(context.blocks))
        result.set("request_memory", parse_request_memory(context.blocks))
        result.set("output_print_settings", parse_output_print_settings(context.blocks))
        return ORCAInpParsePhase.SEMANTICS

    def _run_semantic_phase(
        self, context: ORCAInpParseContext, result: ModelParseResult
    ) -> ORCAInpParsePhase:
        result.update(
            build_orca_input_semantic_payload(
                context.keyword_text, context.blocks, context.geometry
            )
        )
        return ORCAInpParsePhase.DONE

    def _parse_block_to_result(self, block: str) -> ModelParseResult:
        context = ORCAInpParseContext(block=block, lines=block.splitlines())
        result = ModelParseResult()
        phase = ORCAInpParsePhase.STRUCTURE
        while phase is not ORCAInpParsePhase.DONE:
            if phase is ORCAInpParsePhase.STRUCTURE:
                phase = self._run_structure_phase(context, result)
            elif phase is ORCAInpParsePhase.LINES:
                phase = self._run_line_phase(context, result)
            elif phase is ORCAInpParsePhase.GEOMETRY:
                phase = self._run_geometry_phase(context, result)
            elif phase is ORCAInpParsePhase.RESOURCES:
                phase = self._run_resource_phase(context, result)
            elif phase is ORCAInpParsePhase.SEMANTICS:
                phase = self._run_semantic_phase(context, result)
            else:
                raise AssertionError(f"Unexpected ORCA input parse phase: {phase!r}")
        return result

    def _parse_frame(self, block: str, *, context: FrameParseContext) -> Mapping[str, Any]:
        _ = context
        return self._parse_block_to_result(block).model_data()


class ORCAInpFileFrameParserMemory(
    ORCAInpFileFrameParserMixin, BaseFrameParser[ORCAInpFileFrameMemory]
):
    _file_frame_class_ = ORCAInpFileFrameMemory


class ORCAInpFileFrameParserDisk(
    ORCAInpFileFrameParserMixin, BaseFrameParser[ORCAInpFileFrameDisk]
):
    _file_frame_class_ = ORCAInpFileFrameDisk


def parse_orca_input_frame_result(block: str) -> ModelParseResult:
    parser = ORCAInpFileFrameParserMemory()
    return parser._parse_block_to_result(block)
