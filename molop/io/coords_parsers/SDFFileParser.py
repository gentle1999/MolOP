"""
Author: TMJ
Date: 2025-07-30 10:30:16
LastEditors: TMJ
LastEditTime: 2025-08-02 19:57:46
Description: 请填写简介
"""

from typing import Any, Optional, Sequence

from molop.io.base_models.FileParser import BaseFileParserDisk, BaseFileParserMemory
from molop.io.coords_models.SDFFile import SDFFileDisk, SDFFileMemory
from molop.io.coords_models.SDFFileFrame import SDFFileFrameDisk, SDFFileFrameMemory
from molop.io.coords_parsers.SDFFileFrameParser import (
    SDFFileFrameParserDisk,
    SDFFileFrameParserMemory,
)


class SDFFileParserMixin:
    def _parse_metadata(self, file_content: str) -> Optional[dict[str, Any]]: ...

    def _split_file(self, file_content: str) -> Sequence[str]:
        return [block for block in file_content.split("$$$$\n") if block.strip()]


class SDFFileParserMemory(
    SDFFileParserMixin,
    BaseFileParserMemory[SDFFileMemory, SDFFileFrameMemory, SDFFileFrameParserMemory],
):
    _frame_parser = SDFFileFrameParserMemory
    _chem_file = SDFFileMemory


class SDFFileParserDisk(
    SDFFileParserMixin,
    BaseFileParserDisk[SDFFileDisk, SDFFileFrameDisk, SDFFileFrameParserDisk],
):
    _frame_parser = SDFFileFrameParserDisk
    _chem_file = SDFFileDisk
