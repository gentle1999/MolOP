"""
Author: TMJ
Date: 2025-07-29 16:53:34
LastEditors: TMJ
LastEditTime: 2026-02-04 11:00:34
Description: 请填写简介
"""

from typing import Literal, cast

from rdkit import Chem

from molop.io.base_models.ChemFileFrame import BaseCoordsFrame, _HasCoords
from molop.io.base_models.Mixins import DiskStorageMixin, MemoryStorageMixin


class SDFFileFrameMixin:
    def _render(self, engine: Literal["rdkit", "openbabel"] = "rdkit", **kwargs) -> str:
        """
        Render the SDFFileFrame as a string.

        Returns:
            str: The rendered SDFFileFrame.
        """
        typed_self = cast(_HasCoords, self)
        if engine == "rdkit":
            return Chem.MolToMolBlock(typed_self.rdmol) if typed_self.rdmol else ""
        elif engine == "openbabel":
            return cast(str, typed_self.omol.write("sdf")) if typed_self.omol else ""
        else:
            raise ValueError(f"Unsupported engine: {engine}")


class SDFFileFrameMemory(
    MemoryStorageMixin, SDFFileFrameMixin, BaseCoordsFrame["SDFFileFrameMemory"]
): ...


class SDFFileFrameDisk(
    DiskStorageMixin, SDFFileFrameMixin, BaseCoordsFrame["SDFFileFrameDisk"]
): ...
