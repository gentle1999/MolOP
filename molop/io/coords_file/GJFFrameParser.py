"""
Author: TMJ
Date: 2024-06-20 22:47:26
LastEditors: TMJ
LastEditTime: 2024-06-25 13:46:55
Description: 请填写简介
"""

import re
from typing import Literal

import numpy as np
from pydantic import Field, computed_field
from rdkit import Chem

from molop.io.bases.BaseMolFrameParser import BaseMolFrameParser
from molop.unit import atom_ureg
from molop.utils.g16patterns import (
    get_solvent,
    get_solvent_model,
    link0_parser,
    parameter_comment_parser,
)


class GJFFrameParser(BaseMolFrameParser):
    _frame_type = "G16 GJF"
    keywords: str = Field(
        default="",
        description="Parameter comment",
    )

    @property
    def link(self) -> str:
        """
        Get the link.
        """
        return self.keywords.split("#")[0]

    @property
    def route(self):
        return "#" + self.keywords.split("#")[1].replace("\n", " ")

    @property
    def link0(self) -> dict:
        """
        Get the link0.
        """
        return link0_parser(self.link)

    @property
    def route_params(self) -> dict:
        """
        Get the route parameters.
        """
        return parameter_comment_parser(self.route)[0]

    @property
    def dieze_tag(self) -> Literal["#N", "#P", "#T"]:
        return parameter_comment_parser(self.route)[1]

    @property
    def solvent_model(self) -> str:
        return get_solvent_model(self.route_params)

    @property
    def solvent(self) -> str:
        return get_solvent(self.route_params)

    def _parse(self):
        """
        Parse the block.
        """
        lines = self.frame_content.split("\n")
        temp_coords = []
        for line in lines[1:]:
            if re.search(
                r"[a-zA-z0-9]+\s+([\s\-]\d+\.\d+)\s+([\s\-]\d+\.\d+)\s+([\s\-]\d+\.\d+)",
                line,
            ):
                atom, x, y, z = line.split()
                self.atoms.append(Chem.Atom(atom).GetAtomicNum())
                temp_coords.append((float(x), float(y), float(z)))
        self.coords = np.array(temp_coords) * atom_ureg.angstrom
