"""
Author: TMJ
Date: 2025-07-28 23:05:56
LastEditors: TMJ
LastEditTime: 2026-02-04 11:49:07
Description: 请填写简介
"""

from typing import cast

from pydantic import Field

from molop.io.base_models.ChemFileFrame import BaseCoordsFrame, _HasCoords
from molop.io.base_models.Mixins import DiskStorageMixin, MemoryStorageMixin


class XYZFileFrameMixin:
    comment: str = Field(default="", description="comment")

    def _render(self, **kwargs) -> str:
        """Render the XYZ file frame as a string.

        Returns:
            str: The rendered XYZ file frame.
        """
        typed_self = cast(_HasCoords, self)
        return (
            f"{len(typed_self.atoms)}\n"
            + (
                f"{self.comment}\n"
                if self.comment
                else f"charge {typed_self.charge} multiplicity {typed_self.multiplicity}\n"
            )
            + "\n".join(
                [
                    f"{atom:10s}{x:18.10f}{y:18.10f}{z:18.10f}"
                    for atom, (x, y, z) in zip(
                        typed_self.atom_symbols, typed_self.coords.m, strict=True
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
