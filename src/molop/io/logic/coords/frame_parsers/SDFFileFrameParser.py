"""
Author: TMJ
Date: 2025-07-30 10:30:03
LastEditors: TMJ
LastEditTime: 2026-02-04 15:13:04
Description: 请填写简介
"""

from collections.abc import Mapping
from typing import Any

from molop.io.base_models.FrameParser import BaseFrameParser, FrameParseContext
from molop.io.base_models.ParseContainers import ModelParseResult
from molop.io.logic.coords.frame_models.SDFFileFrame import SDFFileFrameDisk, SDFFileFrameMemory
from molop.io.logic.coords.frame_parsers._coords_extractors import extract_sdf_frame_payload


class SDFFileFrameParserMixin:
    def _parse_block_to_result(self, block: str) -> ModelParseResult:
        return ModelParseResult(extract_sdf_frame_payload(block))

    def _parse_frame(self, block: str, *, context: FrameParseContext) -> Mapping[str, Any]:
        _ = context
        return self._parse_block_to_result(block).model_data()


class SDFFileFrameParserMemory(SDFFileFrameParserMixin, BaseFrameParser[SDFFileFrameMemory]):
    _file_frame_class_ = SDFFileFrameMemory


class SDFFileFrameParserDisk(SDFFileFrameParserMixin, BaseFrameParser[SDFFileFrameDisk]):
    _file_frame_class_ = SDFFileFrameDisk
