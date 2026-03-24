"""
Author: TMJ
Date: 2025-10-09 15:23:26
LastEditors: TMJ
LastEditTime: 2026-03-23 22:44:13
Description: 请填写简介
"""

import re
from collections.abc import Mapping
from typing import Any, cast

from rdkit import Chem

from molop.io.base_models.FrameParser import BaseFrameParser, _HasParseMethod
from molop.io.logic.qminput_frame_models.GJFFileFrame import (
    GJFFileFrameDisk,
    GJFFileFrameMemory,
    GJFLink0Commands,
    GJFMoleculeSpecifications,
    GJFRouteSection,
    GJFTitleCard,
)


pt = Chem.GetPeriodicTable()


class GJFFileFrameParserMixin:
    _charge_multiplicity_pattern = re.compile(r"^\s*[+-]?\d+(?:\s+[+-]?\d+)+\s*$")
    _zmat_variable_assignment_pattern = re.compile(
        r"^[A-Za-z][A-Za-z0-9_]*\s+[-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?$"
    )

    @classmethod
    def _normalize_free_separators(cls, line: str) -> str:
        return line.replace("\t", " ").replace(",", " ").replace("/", " ")

    @classmethod
    def _strip_comment_content(cls, line: str) -> str:
        return line.split("!", 1)[0].rstrip()

    @classmethod
    def _is_blank(cls, line: str) -> bool:
        return cls._strip_comment_content(line).strip() == ""

    @classmethod
    def _join_lines(cls, lines: list[str], start: int, end: int) -> str:
        if start >= end:
            return ""
        return "\n".join(lines[start:end])

    @classmethod
    def _find_first_route_index(cls, lines: list[str]) -> int | None:
        for idx, line in enumerate(lines):
            if cls._strip_comment_content(line).lstrip().startswith("#"):
                return idx
        return None

    @classmethod
    def _find_next_blank_index(cls, lines: list[str], start: int) -> int | None:
        for idx in range(start, len(lines)):
            if cls._is_blank(lines[idx]):
                return idx
        return None

    @classmethod
    def _skip_blank_indices(cls, lines: list[str], start: int) -> int:
        idx = start
        while idx < len(lines) and cls._is_blank(lines[idx]):
            idx += 1
        return idx

    @classmethod
    def _is_charge_multiplicity_line(cls, line: str) -> bool:
        stripped = " ".join(
            cls._normalize_free_separators(cls._strip_comment_content(line)).split()
        )
        if not stripped:
            return False
        parts = stripped.split()
        if len(parts) < 2 or len(parts) % 2 != 0:
            return False
        return all(re.match(r"^[+-]?\d+$", part) is not None for part in parts)

    @classmethod
    def _find_first_charge_multiplicity_index(cls, lines: list[str], start: int) -> int | None:
        for idx in range(start, len(lines)):
            if cls._is_charge_multiplicity_line(lines[idx]):
                return idx
        return None

    @classmethod
    def _is_zmat_variable_label(cls, line: str) -> bool:
        stripped = cls._strip_comment_content(line).strip().lower()
        return stripped in {"variables:", "constants:"}

    @classmethod
    def _is_zmat_variable_assignment(cls, line: str) -> bool:
        normalized = " ".join(
            cls._normalize_free_separators(cls._strip_comment_content(line)).split()
        )
        return cls._zmat_variable_assignment_pattern.match(normalized) is not None

    @classmethod
    def _route_flags(cls, route_raw: str) -> tuple[bool, bool]:
        normalized = re.sub(r"\s+", "", route_raw.lower())
        has_allcheck = "allcheck" in normalized
        has_checkpoint = any(
            key in normalized
            for key in ("geom=check", "geom=checkpoint", "geom=(check", "geom=(checkpoint")
        )
        return has_allcheck, has_checkpoint

    @classmethod
    def _find_molecule_block_end(cls, lines: list[str], charge_idx: int) -> tuple[int, int]:
        first_blank = cls._find_next_blank_index(lines, charge_idx)
        molecule_end = first_blank if first_blank is not None else len(lines)
        additional_start = molecule_end + 1 if first_blank is not None else molecule_end

        probe_idx = cls._skip_blank_indices(lines, additional_start)
        if probe_idx >= len(lines):
            return molecule_end, additional_start
        if not (
            cls._is_zmat_variable_label(lines[probe_idx])
            or cls._is_zmat_variable_assignment(lines[probe_idx])
        ):
            return molecule_end, additional_start

        scan_idx = probe_idx
        seen_assignment = False
        while scan_idx < len(lines):
            stripped = lines[scan_idx].strip()
            if not stripped:
                scan_idx += 1
                continue
            if cls._is_zmat_variable_label(stripped):
                scan_idx += 1
                continue
            if cls._is_zmat_variable_assignment(stripped):
                seen_assignment = True
                scan_idx += 1
                continue
            break

        if seen_assignment:
            return scan_idx, scan_idx
        return molecule_end, additional_start

    def _parse_frame(self) -> Mapping[str, Any]:
        typed_self = cast(_HasParseMethod, self)
        block = typed_self._block
        metadata: dict[str, Any] = {}
        lines = [self._strip_comment_content(line) for line in block.splitlines()]

        route_start = self._find_first_route_index(lines)
        if route_start is None:
            metadata["additional_sections"] = block.strip("\n")
            return metadata

        link0_raw = self._join_lines(lines, 0, route_start)
        metadata["link0_commands"] = GJFLink0Commands.from_str(link0_raw)

        route_blank_idx = self._find_next_blank_index(lines, route_start)
        route_end = route_blank_idx if route_blank_idx is not None else len(lines)
        route_raw = self._join_lines(lines, route_start, route_end)
        metadata["route_section"] = GJFRouteSection.from_str(route_raw)
        has_allcheck, has_checkpoint = self._route_flags(route_raw)

        if route_blank_idx is None:
            metadata["additional_sections"] = ""
            return metadata

        section_start = self._skip_blank_indices(lines, route_blank_idx + 1)
        if section_start >= len(lines):
            metadata["additional_sections"] = ""
            return metadata

        charge_idx = self._find_first_charge_multiplicity_index(lines, section_start)
        if charge_idx is None:
            if has_allcheck:
                metadata["additional_sections"] = self._join_lines(lines, section_start, len(lines))
                return metadata

            if has_checkpoint:
                title_end = self._find_next_blank_index(lines, section_start)
                if title_end is None:
                    title_raw = self._join_lines(lines, section_start, len(lines))
                    if title_raw.strip():
                        metadata["title_card"] = GJFTitleCard.from_str(title_raw)
                    metadata["additional_sections"] = ""
                    return metadata

                title_raw = self._join_lines(lines, section_start, title_end)
                if title_raw.strip():
                    metadata["title_card"] = GJFTitleCard.from_str(title_raw)
                additional_start = self._skip_blank_indices(lines, title_end + 1)
                metadata["additional_sections"] = self._join_lines(
                    lines, additional_start, len(lines)
                )
                return metadata

            title_raw = self._join_lines(lines, section_start, len(lines))
            metadata["title_card"] = GJFTitleCard.from_str(title_raw)
            metadata["additional_sections"] = ""
            return metadata

        title_raw = self._join_lines(lines, section_start, charge_idx)
        if not title_raw.strip() and not has_allcheck:
            raise ValueError("GJF title card section is required and cannot be empty")
        if title_raw.strip():
            metadata["title_card"] = GJFTitleCard.from_str(title_raw)

        molecule_end, additional_start = self._find_molecule_block_end(lines, charge_idx)
        molecule_raw = self._join_lines(lines, charge_idx, molecule_end)

        if molecule_raw.strip():
            try:
                metadata["molecule_specifications"] = GJFMoleculeSpecifications.from_str(
                    molecule_raw
                )
            except Exception as exc:
                raise ValueError(f"Failed to parse GJF molecule specifications: {exc}") from exc

        metadata["additional_sections"] = self._join_lines(lines, additional_start, len(lines))
        return metadata


class GJFFileFrameParserMemory(GJFFileFrameParserMixin, BaseFrameParser[GJFFileFrameMemory]):
    _file_frame_class_ = GJFFileFrameMemory


class GJFFileFrameParserDisk(GJFFileFrameParserMixin, BaseFrameParser[GJFFileFrameDisk]):
    _file_frame_class_ = GJFFileFrameDisk
