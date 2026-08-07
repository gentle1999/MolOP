"""
Author: TMJ
Date: 2025-12-14 23:30:32
LastEditors: TMJ
LastEditTime: 2026-02-04 15:15:05
Description: 请填写简介
"""

from collections.abc import Mapping
from typing import Any

from molop.io.base_models.FrameParser import BaseFrameParser, FrameParseContext
from molop.io.base_models.ParseContainers import ModelParseResult
from molop.io.logic.coords.frame_models.SMIFileFrame import SMIFileFrameDisk, SMIFileFrameMemory
from molop.io.logic.coords.frame_parsers._coords_extractors import extract_smi_frame_payload


class SMIFileFrameParserMixin:
    def _parse_block_to_result(self, block: str) -> ModelParseResult:
        return ModelParseResult(extract_smi_frame_payload(block))

    def _parse_frame(self, block: str, *, context: FrameParseContext) -> Mapping[str, Any]:
        _ = context
        return self._parse_block_to_result(block).model_data()


class SMIFileFrameParserMemory(SMIFileFrameParserMixin, BaseFrameParser[SMIFileFrameMemory]):
    _file_frame_class_ = SMIFileFrameMemory


class SMIFileFrameParserDisk(SMIFileFrameParserMixin, BaseFrameParser[SMIFileFrameDisk]):
    _file_frame_class_ = SMIFileFrameDisk
