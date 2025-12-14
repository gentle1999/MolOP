"""
Author: TMJ
Date: 2025-07-28 23:05:56
LastEditors: TMJ
LastEditTime: 2025-12-14 21:30:22
Description: 请填写简介
"""

from typing import TYPE_CHECKING, Protocol

from pint.facets.numpy.quantity import NumpyQuantity
from pydantic import Field

from molop.io.base_models.Bases import BaseDataClassWithUnit
from molop.io.base_models.ChemFileFrame import BaseCoordsFrame
from molop.io.base_models.Mixins import DiskStorageMixin, MemoryStorageMixin


class XYZFileFrameProtocol(Protocol):
    atoms: list[int]
    atom_symbols: list[str]
    coords: NumpyQuantity
    charge: int
    multiplicity: int


if TYPE_CHECKING:

    class _XYZFileFrameProtocol(XYZFileFrameProtocol, BaseDataClassWithUnit): ...
else:

    class _XYZFileFrameProtocol(BaseDataClassWithUnit): ...


class XYZFileFrameMixin(_XYZFileFrameProtocol):
    comment: str = Field(default="", description="comment")

    def _render(self, **kwargs) -> str:
        """Render the XYZ file frame as a string.

        Returns:
            str: The rendered XYZ file frame.
        """
        return (
            f"{len(self.atoms)}\n"
            + (
                f"{self.comment}\n"
                if self.comment
                else f"charge {self.charge} multiplicity {self.multiplicity}\n"
            )
            + "\n".join(
                [
                    f"{atom:10s}{x:18.10f}{y:18.10f}{z:18.10f}"
                    for atom, (x, y, z) in zip(
                        self.atom_symbols, self.coords.m, strict=True
                    )
                ]
            )
        )


class XYZFileFrameMemory(
    MemoryStorageMixin, XYZFileFrameMixin, BaseCoordsFrame["XYZFileFrameMemory"]
): ...


class XYZFileFrameDisk(
    DiskStorageMixin, XYZFileFrameMixin, BaseCoordsFrame["XYZFileFrameDisk"]
): ...
