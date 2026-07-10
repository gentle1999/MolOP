"""
Author: TMJ
Date: 2025-07-29 16:53:34
LastEditors: TMJ
LastEditTime: 2026-02-12 20:33:42
Description: 请填写简介
"""

from typing import Literal, cast

from molop.io.base_models.ChemFileFrame import BaseCoordsFrame, _HasCoords
from molop.io.base_models.Mixins import DiskStorageMixin, MemoryStorageMixin
from molop.io.logic.coords.frame_models._coords_renderers import render_sdf_frame


class SDFFileFrameMixin:
    def _render(self, engine: Literal["rdkit", "openbabel"] = "rdkit", **kwargs) -> str:
        """
        Render the SDFFileFrame as a string.

        Args:
            engine (Literal["rdkit", "openbabel"], optional): The engine to use for rendering. Defaults to "rdkit".

        Returns:
            str: The rendered SDFFileFrame.
        """
        typed_self = cast(_HasCoords, self)
        return render_sdf_frame(typed_self, engine=engine)


class SDFFileFrameMemory(
    MemoryStorageMixin, SDFFileFrameMixin, BaseCoordsFrame["SDFFileFrameMemory"]
): ...


class SDFFileFrameDisk(
    DiskStorageMixin, SDFFileFrameMixin, BaseCoordsFrame["SDFFileFrameDisk"]
): ...
