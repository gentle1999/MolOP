"""
Author: TMJ
Date: 2025-07-29 16:54:42
LastEditors: TMJ
LastEditTime: 2025-07-29 16:58:23
Description: 请填写简介
"""

from molop.io.base_models.ChemFile import CoordsFileMixin
from molop.io.base_models.Mixins import DiskStorageMixin, MemoryStorageMixin
from molop.io.coords_models.XYZFileFrame import XYZFileFrameDisk, XYZFileFrameMemory


class XYZFileMemory(MemoryStorageMixin, CoordsFileMixin[XYZFileFrameMemory]): ...


class XYZFileDisk(DiskStorageMixin, CoordsFileMixin[XYZFileFrameDisk]):
    _allowed_formats_ = ("xyz",)
