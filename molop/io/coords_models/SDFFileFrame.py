"""
Author: TMJ
Date: 2025-07-29 16:53:34
LastEditors: TMJ
LastEditTime: 2025-08-01 14:18:35
Description: 请填写简介
"""

from molop.io.base_models.ChemFileFrame import BaseCoordsFrame
from molop.io.base_models.Mixins import DiskStorageMixin, MemoryStorageMixin
from molop.unit import atom_ureg


class SDFFileFrameMemory(MemoryStorageMixin, BaseCoordsFrame["SDFFileFrameMemory"]):

    def _add_default_units(self) -> None:
        super()._add_default_units()
        self._default_units.update({"coords": atom_ureg.angstrom})


class SDFFileFrameDisk(DiskStorageMixin, BaseCoordsFrame["SDFFileFrameDisk"]):
    _allowed_formats_ = ("sdf", "sd", "mol")

    def _add_default_units(self) -> None:
        super()._add_default_units()
        self._default_units.update({"coords": atom_ureg.angstrom})
