"""
Author: TMJ
Date: 2025-07-29 16:54:42
LastEditors: TMJ
LastEditTime: 2025-09-12 10:55:48
Description: 请填写简介
"""

from molop.io.base_models.ChemFile import BaseCoordsFile
from molop.io.base_models.Mixins import DiskStorageMixin, MemoryStorageMixin
from molop.io.coords_models.XYZFileFrame import XYZFileFrameDisk, XYZFileFrameMemory


class XYZFileMemory(MemoryStorageMixin, BaseCoordsFile[XYZFileFrameMemory]): ...


class XYZFileDisk(DiskStorageMixin, BaseCoordsFile[XYZFileFrameDisk]):
    _allowed_formats_ = ("xyz",)
