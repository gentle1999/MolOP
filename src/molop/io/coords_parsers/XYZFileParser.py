"""
Author: TMJ
Date: 2025-07-29 22:53:30
LastEditors: TMJ
LastEditTime: 2025-08-02 19:58:01
Description: 请填写简介
"""

from collections.abc import Sequence
from typing import Any

from molop.io.base_models.FileParser import BaseFileParserDisk, BaseFileParserMemory
from molop.io.coords_models.XYZFile import XYZFileDisk, XYZFileMemory
from molop.io.coords_models.XYZFileFrame import XYZFileFrameDisk, XYZFileFrameMemory
from molop.io.coords_parsers.XYZFileFrameParser import (
    XYZFileFrameParserDisk,
    XYZFileFrameParserMemory,
)


class XYZFileParserMixin:
    def _parse_metadata(self, file_content: str) -> dict[str, Any] | None: ...

    def _split_file(self, file_content: str) -> Sequence[str]:
        lines = file_content.splitlines()
        anchor = 0
        frames = []
        while anchor < len(lines):
            try:
                num_atoms = int(lines[anchor])
            except ValueError:
                anchor += 1
                continue
            frames.append("\n".join(lines[anchor : anchor + num_atoms + 2]))
            anchor += num_atoms + 2
        return frames


class XYZFileParserMemory(
    XYZFileParserMixin,
    BaseFileParserMemory[XYZFileMemory, XYZFileFrameMemory, XYZFileFrameParserMemory],
):
    _frame_parser = XYZFileFrameParserMemory
    _chem_file = XYZFileMemory


class XYZFileParserDisk(
    XYZFileParserMixin,
    BaseFileParserDisk[XYZFileDisk, XYZFileFrameDisk, XYZFileFrameParserDisk],
):
    allowed_formats = ("xyz",)
    _frame_parser = XYZFileFrameParserDisk
    _chem_file = XYZFileDisk
