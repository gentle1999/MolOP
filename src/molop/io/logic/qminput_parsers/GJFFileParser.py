"""
Author: TMJ
Date: 2025-07-30 14:30:03
LastEditors: TMJ
LastEditTime: 2026-03-23 19:01:43
Description: 请填写简介
"""

from __future__ import annotations

import os
import re
from collections.abc import Sequence
from typing import TYPE_CHECKING

from molop.io.base_models.FileParser import BaseFileParserDisk, BaseFileParserMemory
from molop.io.logic.qminput_frame_models.GJFFileFrame import GJFFileFrameDisk, GJFFileFrameMemory
from molop.io.logic.qminput_frame_parsers.GJFFileFrameParser import (
    GJFFileFrameParserDisk,
    GJFFileFrameParserMemory,
)
from molop.io.logic.qminput_models.GJFFile import GJFFileDisk, GJFFileMemory


if TYPE_CHECKING:
    from molop.io.codec_registry import Registry


class GJFFileParserMixin:
    def _expand_at_includes(self, file_content: str) -> str:
        normalized_content = file_content.replace("\r\n", "\n").replace("\r", "\n")

        source_path = getattr(self, "_file_path", None)
        base_dir: str | None = os.path.dirname(source_path) if source_path else None

        def _expand(content: str, visited: set[str]) -> str:
            expanded_lines: list[str] = []
            for raw_line in content.splitlines():
                stripped = raw_line.strip()
                if stripped.startswith("@") and len(stripped) > 1:
                    include_path_token = stripped[1:].strip()
                    if include_path_token.startswith(('"', "'")) and include_path_token.endswith(
                        ('"', "'")
                    ):
                        include_path_token = include_path_token[1:-1].strip()
                    if not include_path_token:
                        raise ValueError("Include syntax `@filename` requires a non-empty filename")

                    if base_dir is None:
                        raise ValueError(
                            "Cannot resolve `@filename` include without source file path context"
                        )

                    include_path = os.path.abspath(os.path.join(base_dir, include_path_token))
                    if include_path in visited:
                        raise ValueError(
                            f"Detected recursive `@filename` include cycle at {include_path}"
                        )
                    if not os.path.exists(include_path):
                        raise ValueError(f"Included file does not exist: {include_path}")
                    if not os.path.isfile(include_path):
                        raise ValueError(f"Included path is not a file: {include_path}")

                    with open(include_path) as f:
                        include_content = f.read().replace("\r\n", "\n").replace("\r", "\n")
                    expanded_lines.append(_expand(include_content, visited | {include_path}))
                    continue

                expanded_lines.append(raw_line)

            return "\n".join(expanded_lines)

        initial_visited = {source_path} if source_path else set()
        return _expand(normalized_content, initial_visited)

    def _parse_metadata(self, file_content: str):
        return {
            "qm_software": "Gaussian",
            "qm_software_version": "Any",
        }

    def _split_file(self, file_content: str) -> Sequence[str]:
        normalized_content = self._expand_at_includes(file_content)
        for match in re.finditer(r"--[lL][iI][nN][kK]1--", normalized_content):
            marker_start = match.start()
            if marker_start < 2 or normalized_content[marker_start - 2 : marker_start] != "\n\n":
                raise ValueError("`--Link1--` must be preceded by a blank line")
        return [
            frame.strip() + "\n\n"
            for frame in re.split(r"--[lL][iI][nN][kK]1--", normalized_content)
        ]


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
