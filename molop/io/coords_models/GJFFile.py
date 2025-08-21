'''
Author: TMJ
Date: 2025-07-29 17:00:29
LastEditors: TMJ
LastEditTime: 2025-07-29 23:21:54
Description: 请填写简介
'''
from pydantic import Field

from molop.io.base_models.Bases import BaseDataClassWithUnit
from molop.io.base_models.ChemFile import CoordsFileMixin
from molop.io.base_models.Mixins import DiskStorageMixin, MemoryStorageMixin
from molop.io.coords_models.GJFFileFrame import GJFFileFrameDisk, GJFFileFrameMemory


class GJFFileMixin(BaseDataClassWithUnit):
    options: str = Field(default="", description="options comment")
    route: str = Field(default="", description="route comment")
    title_card: str = Field(default="", description="title card")
    suffix: str = Field(default="", description="suffix comment")


class GJFFileMemory(GJFFileMixin, MemoryStorageMixin, CoordsFileMixin[GJFFileFrameMemory]): ...


class GJFFileDisk(GJFFileMixin, DiskStorageMixin, CoordsFileMixin[GJFFileFrameDisk]):
    _allowed_formats_ = ("gjf", "gif", "com", ".gau", ".gjc")
