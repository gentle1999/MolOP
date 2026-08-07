"""
Author: TMJ
Date: 2025-10-09 15:23:26
LastEditors: TMJ
LastEditTime: 2026-03-23 22:44:13
Description: 请填写简介
"""

from collections.abc import Mapping
from dataclasses import dataclass
from enum import Enum, auto
from typing import Any

from molop.io.base_models.FrameParser import BaseFrameParser, FrameParseContext
from molop.io.base_models.ParseContainers import ModelParseResult
from molop.io.logic.gaussian.input.frame_models.GJFFileFrame import (
    GJFFileFrameDisk,
    GJFFileFrameMemory,
)
from molop.io.logic.gaussian.input.frame_parsers._gjf_extractors import (
    build_gjf_context_lines,
    find_first_gjf_charge_multiplicity_index,
    find_first_gjf_route_index,
    find_gjf_molecule_block_end,
    find_next_gjf_blank_index,
    gjf_route_flags,
    join_gjf_lines,
    skip_gjf_blank_indices,
)
from molop.io.logic.gaussian.input.GaussianInputParsing import (
    parse_gjf_additional_sections,
    parse_gjf_link0_commands,
    parse_gjf_molecule_specifications,
    parse_gjf_route_section,
    parse_gjf_title_card,
)


class GJFParsePhase(Enum):
    """Explicit stages for Gaussian input frame parsing."""

    PREAMBLE = auto()
    TITLE = auto()
    MOLECULE = auto()
    ADDITIONAL = auto()
    DONE = auto()


@dataclass(slots=True)
class GJFParseContext:
    """Mutable parse context shared by Gaussian input parser phases."""

    block: str
    lines: list[str]
    route_start: int | None = None
    route_blank_idx: int | None = None
    route_raw: str = ""
    has_allcheck: bool = False
    has_checkpoint: bool = False
    section_start: int = 0
    charge_idx: int | None = None
    molecule_end: int = 0
    additional_start: int = 0


class GJFFileFrameParserMixin:
    @classmethod
    def _set_additional_sections(cls, result: ModelParseResult, raw: str) -> None:
        parsed_sections, diagnostics = parse_gjf_additional_sections(raw)
        result.set("additional_sections", raw)
        result.set("parsed_additional_sections", parsed_sections)
        result.set("additional_section_diagnostics", diagnostics)

    def _run_preamble_phase(
        self, context: GJFParseContext, result: ModelParseResult
    ) -> GJFParsePhase:
        context.route_start = find_first_gjf_route_index(context.lines)
        if context.route_start is None:
            self._set_additional_sections(result, context.block.strip("\n"))
            return GJFParsePhase.DONE

        link0_raw = join_gjf_lines(context.lines, 0, context.route_start)
        result.set("link0_commands", parse_gjf_link0_commands(link0_raw))

        context.route_blank_idx = find_next_gjf_blank_index(context.lines, context.route_start)
        route_end = (
            context.route_blank_idx if context.route_blank_idx is not None else len(context.lines)
        )
        context.route_raw = join_gjf_lines(context.lines, context.route_start, route_end)
        result.set("route_section", parse_gjf_route_section(context.route_raw))
        context.has_allcheck, context.has_checkpoint = gjf_route_flags(context.route_raw)

        if context.route_blank_idx is None:
            self._set_additional_sections(result, "")
            return GJFParsePhase.DONE

        context.section_start = skip_gjf_blank_indices(context.lines, context.route_blank_idx + 1)
        if context.section_start >= len(context.lines):
            self._set_additional_sections(result, "")
            return GJFParsePhase.DONE

        return GJFParsePhase.TITLE

    def _run_title_phase(self, context: GJFParseContext, result: ModelParseResult) -> GJFParsePhase:
        context.charge_idx = find_first_gjf_charge_multiplicity_index(
            context.lines, context.section_start
        )
        if context.charge_idx is None:
            if context.has_allcheck:
                self._set_additional_sections(
                    result,
                    join_gjf_lines(context.lines, context.section_start, len(context.lines)),
                )
                return GJFParsePhase.DONE

            if context.has_checkpoint:
                title_end = find_next_gjf_blank_index(context.lines, context.section_start)
                if title_end is None:
                    title_raw = join_gjf_lines(
                        context.lines,
                        context.section_start,
                        len(context.lines),
                    )
                    if title_raw.strip():
                        result.set("title_card", parse_gjf_title_card(title_raw))
                    self._set_additional_sections(result, "")
                    return GJFParsePhase.DONE

                title_raw = join_gjf_lines(context.lines, context.section_start, title_end)
                if title_raw.strip():
                    result.set("title_card", parse_gjf_title_card(title_raw))
                context.additional_start = skip_gjf_blank_indices(context.lines, title_end + 1)
                self._set_additional_sections(
                    result,
                    join_gjf_lines(context.lines, context.additional_start, len(context.lines)),
                )
                return GJFParsePhase.DONE

            title_raw = join_gjf_lines(context.lines, context.section_start, len(context.lines))
            result.set("title_card", parse_gjf_title_card(title_raw))
            self._set_additional_sections(result, "")
            return GJFParsePhase.DONE

        title_raw = join_gjf_lines(context.lines, context.section_start, context.charge_idx)
        if not title_raw.strip() and not context.has_allcheck:
            raise ValueError("GJF title card section is required and cannot be empty")
        if title_raw.strip():
            result.set("title_card", parse_gjf_title_card(title_raw))
        return GJFParsePhase.MOLECULE

    def _run_molecule_phase(
        self, context: GJFParseContext, result: ModelParseResult
    ) -> GJFParsePhase:
        assert context.charge_idx is not None
        context.molecule_end, context.additional_start = find_gjf_molecule_block_end(
            context.lines, context.charge_idx
        )
        molecule_raw = join_gjf_lines(context.lines, context.charge_idx, context.molecule_end)

        if molecule_raw.strip():
            try:
                result.set(
                    "molecule_specifications",
                    parse_gjf_molecule_specifications(molecule_raw),
                )
            except Exception as exc:
                raise ValueError(f"Failed to parse GJF molecule specifications: {exc}") from exc

        return GJFParsePhase.ADDITIONAL

    def _run_additional_phase(
        self, context: GJFParseContext, result: ModelParseResult
    ) -> GJFParsePhase:
        self._set_additional_sections(
            result,
            join_gjf_lines(context.lines, context.additional_start, len(context.lines)),
        )
        return GJFParsePhase.DONE

    def _parse_block_to_result(self, block: str) -> ModelParseResult:
        context = GJFParseContext(
            block=block,
            lines=build_gjf_context_lines(block),
        )
        result = ModelParseResult()
        phase = GJFParsePhase.PREAMBLE
        while phase is not GJFParsePhase.DONE:
            if phase is GJFParsePhase.PREAMBLE:
                phase = self._run_preamble_phase(context, result)
            elif phase is GJFParsePhase.TITLE:
                phase = self._run_title_phase(context, result)
            elif phase is GJFParsePhase.MOLECULE:
                phase = self._run_molecule_phase(context, result)
            elif phase is GJFParsePhase.ADDITIONAL:
                phase = self._run_additional_phase(context, result)
            else:
                raise AssertionError(f"Unexpected GJF parse phase: {phase!r}")
        return result

    def _parse_frame(self, block: str, *, context: FrameParseContext) -> Mapping[str, Any]:
        _ = context
        return self._parse_block_to_result(block).model_data()


class GJFFileFrameParserMemory(GJFFileFrameParserMixin, BaseFrameParser[GJFFileFrameMemory]):
    _file_frame_class_ = GJFFileFrameMemory


class GJFFileFrameParserDisk(GJFFileFrameParserMixin, BaseFrameParser[GJFFileFrameDisk]):
    _file_frame_class_ = GJFFileFrameDisk
