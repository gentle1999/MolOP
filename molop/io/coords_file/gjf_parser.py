"""
Author: TMJ
Date: 2024-01-09 19:54:01
LastEditors: TMJ
LastEditTime: 2024-01-11 10:41:27
Description: 请填写简介
"""
import os
import re

from molop.io.bases.file_base import BaseFileParser
from molop.io.coords_file.GJFBlockParser import GJFBlockParser


class GJFParser(BaseFileParser):
    """
    Parser for GJF files.
    """

    _allowed_formats = (".gjf",)

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
                self._parameter_comment = "".join(lines[:idx])

            if re.match(r"^\s*[A-Z][a-z]?(\s+\-?\d+(\.\d+)?){3}$", line):
                block_end = idx
        self.append(
            GJFBlockParser(
                "".join(lines[block_start : block_end + 1]),
                charge,
                multi,
                parameter_comment=self._parameter_comment,
            )
        )

    @property
    def parameter_comment(self) -> str:
        return self._parameter_comment
