"""
Author: TMJ
Date: 2025-07-28 23:05:56
LastEditors: TMJ
LastEditTime: 2025-09-12 10:55:49
Description: 请填写简介
"""

from pydantic import Field

from molop.io.base_models.ChemFileFrame import BaseCoordsFrame
from molop.io.base_models.Mixins import DiskStorageMixin, MemoryStorageMixin
from molop.io.base_models.Molecule import Molecule
from molop.unit import atom_ureg


class XYZFileFrameMixin(Molecule):
    comment: str = Field(default="", description="comment")


class XYZFileFrameMemory(
    MemoryStorageMixin, XYZFileFrameMixin, BaseCoordsFrame["XYZFileFrameMemory"]
):

    def _add_default_units(self) -> None:
        super()._add_default_units()
        self._default_units.update({"coords": atom_ureg.angstrom})


class XYZFileFrameDisk(
    DiskStorageMixin, XYZFileFrameMixin, BaseCoordsFrame["XYZFileFrameDisk"]
):
    _allowed_formats_ = ("xyz",)

    def _add_default_units(self) -> None:
        super()._add_default_units()
        self._default_units.update({"coords": atom_ureg.angstrom})
