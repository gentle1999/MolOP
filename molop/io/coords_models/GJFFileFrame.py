'''
Author: TMJ
Date: 2025-07-28 23:09:31
LastEditors: TMJ
LastEditTime: 2025-08-03 21:09:24
Description: 请填写简介
'''

from pydantic import Field
from rdkit import Chem

from molop.io.base_models.ChemFileFrame import BaseCoordsFrame
from molop.io.base_models.Mixins import DiskStorageMixin, MemoryStorageMixin
from molop.io.base_models.Molecule import Molecule
from molop.io.patterns.G16Patterns import options_parser


class GJFFileFrameMixin(Molecule):
    options: str = Field(default="", description="options comment")
    route: str = Field(default="", description="route comment")
    title_card: str = Field(default="", description="title card")
    suffix: str = Field(default="", description="suffix comment")

    def to_GJF_block(self) -> str:
        """
        Convert the GJFFileFrame to a GJF block string.
        """
        return (
            self.options
            + self.route
            + "\n\n"
            + f"{self.title_card}\n\n"
            + f"{self.charge} {self.multiplicity}\n"
            + "\n".join(
                [
                    f"{Chem.Atom(atom).GetSymbol():10s}{x:14.8f}{y:14.8f}{z:14.8f}"
                    for atom, (x, y, z) in zip(self.atoms, self.coords.m, strict=True)
                ]
            )
            + "\n\n"
            + self.suffix
            + "\n\n"
        )


class GJFFileFrameMemory(
    GJFFileFrameMixin, MemoryStorageMixin, BaseCoordsFrame["GJFFileFrameMemory"]
): ...


class GJFFileFrameDisk(
    GJFFileFrameMixin, DiskStorageMixin, BaseCoordsFrame["GJFFileFrameDisk"]
):
    _allowed_formats_ = ("gjf", "gif", "com", ".gau", ".gjc")
    def to_GJF_block(
        self,
        chk: bool = True,
        oldchk: bool = False,
    ) -> str:
        """
        Convert the GJFFileFrame to a GJF block string.
        """
        _options = options_parser(self.options)
        if chk:
            _options[r"%chk"] = f"{self.pure_filename}.chk"
        if oldchk:
            _options[r"%oldchk"] = f"{self.pure_filename}.chk"
        options_lines = (
            "\n".join([f"{key}={val}" for key, val in _options.items()]) + "\n"
            if _options
            else ""
        )
        return (
            options_lines
            + self.route
            + "\n\n"
            + f"{self.title_card}\n\n"
            + f"{self.charge} {self.multiplicity}\n"
            + "\n".join(
                [
                    f"{Chem.Atom(atom).GetSymbol():10s}{x:14.8f}{y:14.8f}{z:14.8f}"
                    for atom, (x, y, z) in zip(self.atoms, self.coords.m, strict=True)
                ]
            )
            + "\n\n"
            + self.suffix
            + "\n\n"
        )
