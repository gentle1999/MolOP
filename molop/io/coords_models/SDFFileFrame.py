"""
Author: TMJ
Date: 2025-07-29 16:53:34
LastEditors: TMJ
LastEditTime: 2025-12-14 21:30:00
Description: 请填写简介
"""

from typing import TYPE_CHECKING, Literal, Optional, Protocol

from rdkit import Chem

from molop.io.base_models.Bases import BaseDataClassWithUnit
from molop.io.base_models.ChemFileFrame import BaseCoordsFrame
from molop.io.base_models.Mixins import DiskStorageMixin, MemoryStorageMixin
from molop.utils.types import OMol, RdMol


class SDFFileFrameProtocol(Protocol):
    rdmol: Optional[RdMol]
    omol: Optional[OMol]


if TYPE_CHECKING:

    class _SDFFileFrameProtocol(SDFFileFrameProtocol, BaseDataClassWithUnit): ...
else:

    class _SDFFileFrameProtocol(BaseDataClassWithUnit): ...


class SDFFileFrameMixin(_SDFFileFrameProtocol):
    def _render(self, engine: Literal["rdkit", "openbabel"] = "rdkit", **kwargs) -> str:
        """
        Render the SDFFileFrame as a string.

        Returns:
            str: The rendered SDFFileFrame.
        """
        if engine == "rdkit":
            return Chem.MolToMolBlock(self.rdmol) if self.rdmol else ""
        elif engine == "openbabel":
            return self.omol.write("sdf") if self.omol else ""  # type: ignore
        else:
            raise ValueError(f"Unsupported engine: {engine}")


class SDFFileFrameMemory(
    MemoryStorageMixin, SDFFileFrameMixin, BaseCoordsFrame["SDFFileFrameMemory"]
): ...


class SDFFileFrameDisk(
    DiskStorageMixin, SDFFileFrameMixin, BaseCoordsFrame["SDFFileFrameDisk"]
): ...
