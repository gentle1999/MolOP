"""
Author: TMJ
Date: 2025-07-29 16:53:34
LastEditors: TMJ
LastEditTime: 2025-11-21 16:41:26
Description: 请填写简介
"""

from molop.io.base_models.Bases import BaseDataClassWithUnit
from molop.io.base_models.ChemFileFrame import BaseCoordsFrame
from molop.io.base_models.Mixins import DiskStorageWithFrameMixin, MemoryStorageMixin


class SDFFileFrameMixin(BaseDataClassWithUnit): ...


class SDFFileFrameMemory(
    MemoryStorageMixin, SDFFileFrameMixin, BaseCoordsFrame["SDFFileFrameMemory"]
): ...


class SDFFileFrameDisk(
    DiskStorageWithFrameMixin, SDFFileFrameMixin, BaseCoordsFrame["SDFFileFrameDisk"]
):
    _allowed_formats_ = ("sdf", "sd", "mol")
