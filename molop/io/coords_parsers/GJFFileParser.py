"""
Author: TMJ
Date: 2025-07-30 14:30:03
LastEditors: TMJ
LastEditTime: 2025-08-20 15:11:32
Description: 请填写简介
"""

from typing import Any, Sequence

from molop.io.base_models.FileParser import BaseFileParserDisk, BaseFileParserMemory
from molop.io.coords_models.GJFFile import GJFFileDisk, GJFFileMemory
from molop.io.coords_models.GJFFileFrame import GJFFileFrameDisk, GJFFileFrameMemory
from molop.io.coords_parsers.GJFFileFrameParser import (
    GJFFileFrameParserDisk,
    GJFFileFrameParserMemory,
)
from molop.io.patterns.G16Patterns import G16InputPatterns

g16_input_patterns = G16InputPatterns()


class GJFFileParserMixin:
    def _parse_metadata(self, file_content: str) -> dict[str, Any]:
        block = file_content
        metadata: dict[str, Any] = {}
        if matches := g16_input_patterns.OPTIONS.match_content(block):
            options = "\n".join([f"{match[0]}={match[1]}" for match in matches])
            metadata.update({"options": options})
        if indexes := g16_input_patterns.ROUTE.locate_content(block):
            start_start, start_end, end_start, end_end = indexes
            if g16_input_patterns.ROUTE.content_pattern_compiled:
                if matches := g16_input_patterns.ROUTE.content_pattern_compiled.search(
                    block[start_start:end_end]
                ):
                    route = matches.groups()[0].strip()
                    metadata.update({"route": route})
            block = block[end_end:]
        if matches := g16_input_patterns.TITLE.match_content(block):
            title = matches[0][0]
            metadata.update({"title_card": title})
        if matches := g16_input_patterns.CHARGE_MULTIPLICITY.match_content(block):
            charge, multiplicity = matches[0]
            metadata.update({"charge": charge, "multiplicity": multiplicity})
        if g16_input_patterns.ATOMS.content_pattern_compiled:
            while match := g16_input_patterns.ATOMS.content_pattern_compiled.search(
                block
            ):
                block = block[match.end() :]
            suffix = block.strip()
            metadata.update({"suffix": suffix})
        return metadata

    def _split_file(self, file_content: str) -> Sequence[str]:
        return [file_content]


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
    _frame_parser = GJFFileFrameParserDisk
    _chem_file = GJFFileDisk
