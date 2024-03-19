"""
Author: TMJ
Date: 2024-02-17 15:17:37
LastEditors: TMJ
LastEditTime: 2024-02-21 21:19:52
Description: 请填写简介
"""
import re
from typing import Literal

import numpy as np

from molop.io.bases.molblock_base import BaseBlockParser
from molop.unit import atom_ureg
from molop.utils import (
    get_solvent,
    get_solvent_model,
    link0_parser,
    parameter_comment_parser,
)


class GJFBlockParser(BaseBlockParser):
    """
    Parser for XYZ Blocks.
    """

    _block_type = "G16 GJF"

    def __init__(
        self,
        block: str,
        charge=0,
        multiplicity=1,
        file_path="",
        parameter_comment: str = None,
    ):
        super().__init__(block)
        self._file_path = file_path
        self._charge = charge
        self._multiplicity = multiplicity
        self._parameter_comment = parameter_comment
        link, route = self._parameter_comment.split("#")
        link = link.replace("\n", " ")
        route = "#" + route.replace("\n", " ")

        self._link0 = link0_parser(link)
        (
            self._route_params,
            self._dieze_tag,
        ) = parameter_comment_parser(route)
        self._parse()

    @property
    def link0(self) -> dict:
        return self._link0

    @property
    def route_params(self) -> dict:
        return self._route_params

    @property
    def dieze_tag(self) -> Literal["#N", "#P", "#T"]:
        return self._dieze_tag

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
        lines = self._block.split("\n")
        temp_coords = []
        for line in lines[1:]:
            if re.search(
                r"[a-zA-z0-9]+\s+([\s\-]\d+\.\d+)\s+([\s\-]\d+\.\d+)\s+([\s\-]\d+\.\d+)",
                line,
            ):
                atom, x, y, z = line.split()
                self._atoms.append(atom)
                temp_coords.append((float(x), float(y), float(z)))
        self._coords = np.array(temp_coords) * atom_ureg.angstrom

    @property
    def parameter_comment(self) -> str:
        return self._parameter_comment
