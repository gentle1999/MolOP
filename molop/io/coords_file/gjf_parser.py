"""
Author: TMJ
Date: 2024-02-17 15:17:37
LastEditors: TMJ
LastEditTime: 2024-03-18 16:53:46
Description: 请填写简介
"""
import os
import re
from typing import Literal

from molop.io.bases.file_base import BaseFileParser
from molop.io.coords_file.GJFBlockParser import GJFBlockParser
from molop.utils import (
    get_solvent,
    get_solvent_model,
    link0_parser,
    parameter_comment_parser,
)


class GJFParser(BaseFileParser):
    """
    Parser for GJF files.
    """

    _allowed_formats = (".gjf", ".com", ".gau", ".gjc")

    def __init__(self, file_path: str, charge=0, multiplicity=1):
        self._check_formats(file_path)
        super().__init__(file_path)
        self.__force_charge = charge
        self.__force_multiplicity = multiplicity
        self._parse()

    def _parse(self):
        """
        Parse the file.
        """
        with open(self.file_path, "r") as fr:
            lines = fr.readlines()
        fr.close()
        for idx, line in enumerate(lines):
            if re.match(r"^\s*[\+\-\d]+\s+\d+$", line):
                block_start = idx
                charge, multi = map(int, line.split())
                if self.__force_charge is not None:
                    charge = self.__force_charge
                if self.__force_multiplicity is not None:
                    multi = self.__force_multiplicity
                self._parameter_comment = "".join(lines[: idx - 2])

            if re.search(
                r"[a-zA-z0-9]+\s+([\s\-]\d+\.\d+)\s+([\s\-]\d+\.\d+)\s+([\s\-]\d+\.\d+)",
                line,
            ):
                block_end = idx
        link, route = self._parameter_comment.split("#")
        link = link.replace("\n", " ")
        route = "#" + route.replace("\n", " ")

        self._link0 = link0_parser(link)
        (
            self._route_params,
            self._dieze_tag,
        ) = parameter_comment_parser(route)

        self.append(
            GJFBlockParser(
                "".join(lines[block_start : block_end + 1]),
                charge,
                multi,
                file_path=self._file_path,
                parameter_comment=self.parameter_comment,
            )
        )

    @property
    def parameter_comment(self) -> str:
        link, route = self._parameter_comment.split("#")
        link = link.replace("\n", " ")
        route = "#" + route.replace("\n", " ")
        return "\n".join([link, route])

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
