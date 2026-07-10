"""
Author: TMJ
Date: 2025-12-14 23:26:19
LastEditors: TMJ
LastEditTime: 2026-02-04 11:46:52
Description: 请填写简介
"""

from typing import cast

from molop.io.base_models.ChemFileFrame import BaseCoordsFrame, _HasCoords
from molop.io.base_models.Mixins import DiskStorageMixin, MemoryStorageMixin
from molop.io.logic.coords.frame_models._coords_renderers import render_smi_frame


class SMIFileFrameMixin:
    def _render(self, **kwargs) -> str:
        typed_self = cast(_HasCoords, self)
        return render_smi_frame(typed_self)


class SMIFileFrameMemory(
    MemoryStorageMixin, SMIFileFrameMixin, BaseCoordsFrame["SMIFileFrameMemory"]
): ...


class SMIFileFrameDisk(
    DiskStorageMixin, SMIFileFrameMixin, BaseCoordsFrame["SMIFileFrameDisk"]
): ...
