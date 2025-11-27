"""
Author: TMJ
Date: 2025-07-28 23:09:31
LastEditors: TMJ
LastEditTime: 2025-11-21 16:35:26
Description: 请填写简介
"""

from pydantic import Field

from molop.io.base_models.Bases import BaseDataClassWithUnit
from molop.io.base_models.ChemFileFrame import BaseCoordsFrame
from molop.io.base_models.Mixins import DiskStorageWithFrameMixin, MemoryStorageMixin


class GJFFileFrameMixin(BaseDataClassWithUnit):
    options: str = Field(default="", description="options comment")
    route: str = Field(default="", description="route comment")
    title_card: str = Field(default="", description="title card")
    suffix: str = Field(default="", description="suffix comment")


class GJFFileFrameMemory(
    MemoryStorageMixin, GJFFileFrameMixin, BaseCoordsFrame["GJFFileFrameMemory"]
): ...


class GJFFileFrameDisk(
    DiskStorageWithFrameMixin, GJFFileFrameMixin, BaseCoordsFrame["GJFFileFrameDisk"]
):
    _allowed_formats_ = ("gjf", "gif", "com", ".gau", ".gjc")
