"""
Author: TMJ
Date: 2025-07-29 16:59:36
LastEditors: TMJ
LastEditTime: 2025-07-29 17:00:45
Description: 请填写简介
"""

from molop.io.base_models.ChemFile import BaseCoordsFile
from molop.io.base_models.Mixins import DiskStorageMixin, MemoryStorageMixin
from molop.io.coords_models.SDFFileFrame import SDFFileFrameDisk, SDFFileFrameMemory


class SDFFileMemory(MemoryStorageMixin, BaseCoordsFile[SDFFileFrameMemory]): ...


class SDFFileDisk(DiskStorageMixin, BaseCoordsFile[SDFFileFrameDisk]):
    _allowed_formats_ = ("sdf", "sd", "mol")
