"""
Author: TMJ
Date: 2024-10-19 09:57:26
LastEditors: TMJ
LastEditTime: 2025-01-14 20:35:34
Description: 请填写简介
"""

import re
from typing import Literal

from pydantic import Field

from molop.io.bases.BaseMolFileParser import BaseMolFileParser
from molop.io.coords_file.GJFFrameParser import GJFFrameParser
from molop.utils.g16patterns import (
    get_solvent,
    get_solvent_model,
    link0_parser,
    parameter_comment_parser,
)


class GJFFileParser(BaseMolFileParser[GJFFrameParser]):
    _allowed_formats = (".gjf", ".com", ".gau", ".gjc")

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
        with open(self.file_path, "r") as fr:
            lines = fr.readlines()
        fr.close()
        if self.charge is not None:
            charge = self.charge
        if self.multiplicity is not None:
            multi = self.multiplicity
        for idx, line in enumerate(lines):
            if re.match(r"^\s*[\+\-\d]+\s+\d+$", line):
                block_start = idx
                charge, multi = map(int, line.split())
                self.charge = charge
                self.multiplicity = multi
                self.keywords = "".join(lines[: idx - 2])

            if re.search(
                r"[a-zA-z0-9]+\s+([\s\-]\d+\.\d+)\s+([\s\-]\d+\.\d+)\s+([\s\-]\d+\.\d+)",
                line,
            ):
                block_end = idx
        self.__append__(
            GJFFrameParser(
                frame_content="".join(lines[block_start : block_end + 1]),
                charge=charge,
                multiplicity=multi,
                file_path=self.file_path,
                keywords=self.keywords,
            )
        )
