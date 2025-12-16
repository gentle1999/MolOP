"""
Author: TMJ
Date: 2025-12-14 23:30:46
LastEditors: TMJ
LastEditTime: 2025-12-14 23:41:30
Description: 请填写简介
"""

from collections.abc import Sequence
from typing import Any

from molop.io.base_models.FileParser import BaseFileParserDisk, BaseFileParserMemory
from molop.io.coords_models.SMIFile import SMIFileDisk, SMIFileMemory
from molop.io.coords_models.SMIFileFrame import SMIFileFrameDisk, SMIFileFrameMemory
from molop.io.coords_parsers.SMIFileFrameParser import (
    SMIFileFrameParserDisk,
    SMIFileFrameParserMemory,
)


class SMIFileParserMixin:
    def _parse_metadata(self, file_content: str) -> dict[str, Any] | None: ...

    def _split_file(self, file_content: str) -> Sequence[str]:
        return file_content.splitlines()


class SMIFileParserMemory(
    SMIFileParserMixin,
    BaseFileParserMemory[SMIFileMemory, SMIFileFrameMemory, SMIFileFrameParserMemory],
):
    _frame_parser = SMIFileFrameParserMemory
    _chem_file = SMIFileMemory


class SMIFileParserDisk(
    SMIFileParserMixin,
    BaseFileParserDisk[SMIFileDisk, SMIFileFrameDisk, SMIFileFrameParserDisk],
):
    _allowed_formats_ = ("smi", "txt")
    _frame_parser = SMIFileFrameParserDisk
    _chem_file = SMIFileDisk
