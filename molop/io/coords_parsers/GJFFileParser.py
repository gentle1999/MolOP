"""
Author: TMJ
Date: 2025-07-30 14:30:03
LastEditors: TMJ
LastEditTime: 2025-12-14 20:24:50
Description: 请填写简介
"""
import re
from typing import Sequence

from molop.io.base_models.FileParser import BaseFileParserDisk, BaseFileParserMemory
from molop.io.coords_models.GJFFile import GJFFileDisk, GJFFileMemory
from molop.io.coords_models.GJFFileFrame import GJFFileFrameDisk, GJFFileFrameMemory
from molop.io.coords_parsers.GJFFileFrameParser import (
    GJFFileFrameParserDisk,
    GJFFileFrameParserMemory,
)


class GJFFileParserMixin:
    def _parse_metadata(self, file_content: str): ...

    def _split_file(self, file_content: str) -> Sequence[str]:
        return [frame.strip()+"\n\n" for frame in re.split(r"--[lL][iI][nN][kK]1--", file_content)]


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
    _allowed_formats_ = ("gjf", "gif", "com", ".gau", ".gjc")
    _frame_parser = GJFFileFrameParserDisk
    _chem_file = GJFFileDisk
