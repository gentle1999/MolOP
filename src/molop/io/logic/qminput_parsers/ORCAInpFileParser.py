"""
Author: TMJ
Date: 2026-02-10 00:00:00
LastEditors: TMJ
LastEditTime: 2026-02-10 00:00:00
Description: ORCA input file parsers
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING, Any

from molop.io.base_models.FileParser import BaseFileParserDisk, BaseFileParserMemory
from molop.io.logic.qminput_frame_models.ORCAInpFileFrame import (
    ORCAInpFileFrameDisk,
    ORCAInpFileFrameMemory,
)
from molop.io.logic.qminput_frame_parsers.ORCAInpFileFrameParser import (
    ORCAInpFileFrameParserDisk,
    ORCAInpFileFrameParserMemory,
)
from molop.io.logic.qminput_models.ORCAInpFile import ORCAInpFileDisk, ORCAInpFileMemory


if TYPE_CHECKING:
    from molop.io.codec_registry import Registry


_CARTESIAN_TYPES = {"xyz", "cart", "cartesian", "int", "internal", "gzmt", "xyzfile", "gzmtfile"}


def _is_coords_open_line(line: str) -> bool:
    lowered = line.lower()
    return lowered.startswith("coords") and (len(lowered) == 6 or lowered[6].isspace())


def _is_new_job_delimiter(line: str) -> bool:
    lowered = line.lower()
    return lowered == "$new_job" or lowered.startswith("$new_job ")


class ORCAInpFileParserMixin:
    def _parse_metadata(self, file_content: str) -> dict[str, Any] | None:
        _ = file_content
        return None

    def _split_file(self, file_content: str) -> Sequence[str]:
        lines = file_content.splitlines(keepends=True)
        frames: list[str] = []
        current_lines: list[str] = []

        in_star_geometry = False
        in_percent_coords = False
        nested_coords_depth = 0

        for line in lines:
            stripped = line.strip()
            lowered = stripped.lower()

            if not in_star_geometry and not in_percent_coords and _is_new_job_delimiter(lowered):
                frames.append("".join(current_lines))
                current_lines = []
                continue

            current_lines.append(line)

            if not in_percent_coords:
                if in_star_geometry:
                    if stripped == "*":
                        in_star_geometry = False
                elif stripped.startswith("*"):
                    header = stripped[1:].strip()
                    if header:
                        ctype = header.split(maxsplit=1)[0].lower()
                        if ctype in _CARTESIAN_TYPES:
                            in_star_geometry = True

            if not in_star_geometry:
                if not in_percent_coords:
                    if lowered.startswith(r"%coords"):
                        in_percent_coords = True
                        nested_coords_depth = 0
                else:
                    if _is_coords_open_line(lowered):
                        nested_coords_depth += 1
                    elif lowered == "end":
                        if nested_coords_depth > 0:
                            nested_coords_depth -= 1
                        else:
                            in_percent_coords = False

        frames.append("".join(current_lines))
        return [frame for frame in frames if frame]


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
