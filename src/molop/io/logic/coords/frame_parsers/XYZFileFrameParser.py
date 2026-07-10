"""
Author: TMJ
Date: 2025-10-09 15:33:21
LastEditors: TMJ
LastEditTime: 2026-02-05 19:55:18
Description: 请填写简介
"""

from collections.abc import Mapping
from typing import Any, cast

from molop.io.base_models.FrameParser import BaseFrameParser, _HasParseMethod
from molop.io.base_models.ParseContainers import ModelParseResult
from molop.io.logic.coords.frame_models.XYZFileFrame import XYZFileFrameDisk, XYZFileFrameMemory
from molop.io.logic.coords.frame_parsers._coords_extractors import extract_xyz_frame_payload


class XYZFileFrameParserMixin:
    def _parse_block_to_result(self, block: str) -> ModelParseResult:
        return ModelParseResult(extract_xyz_frame_payload(block))

    def _parse_frame(self) -> Mapping[str, Any]:
        typed_self = cast(_HasParseMethod, self)
        return self._parse_block_to_result(typed_self._block).model_data()


class XYZFileFrameParserMemory(XYZFileFrameParserMixin, BaseFrameParser[XYZFileFrameMemory]):
    _file_frame_class_ = XYZFileFrameMemory


class XYZFileFrameParserDisk(XYZFileFrameParserMixin, BaseFrameParser[XYZFileFrameDisk]):
    _file_frame_class_ = XYZFileFrameDisk
