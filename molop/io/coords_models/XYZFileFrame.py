"""
Author: TMJ
Date: 2025-07-28 23:05:56
LastEditors: TMJ
LastEditTime: 2025-11-21 16:43:16
Description: 请填写简介
"""

from pydantic import Field

from molop.io.base_models.Bases import BaseDataClassWithUnit
from molop.io.base_models.ChemFileFrame import BaseCoordsFrame
from molop.io.base_models.Mixins import DiskStorageWithFrameMixin, MemoryStorageMixin


class XYZFileFrameMixin(BaseDataClassWithUnit):
    comment: str = Field(default="", description="comment")


class XYZFileFrameMemory(
    MemoryStorageMixin, XYZFileFrameMixin, BaseCoordsFrame["XYZFileFrameMemory"]
): ...


class XYZFileFrameDisk(
    DiskStorageWithFrameMixin, XYZFileFrameMixin, BaseCoordsFrame["XYZFileFrameDisk"]
):
    _allowed_formats_ = ("xyz",)
