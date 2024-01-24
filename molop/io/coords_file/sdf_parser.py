"""
Author: TMJ
Date: 2023-10-30 18:21:31
LastEditors: TMJ
LastEditTime: 2024-01-24 22:12:28
Description: 请填写简介
"""

import os

from openbabel import pybel

from molop.io.bases.file_base import BaseFileParser
from molop.io.coords_file.SDFBlockParser import SDFBlockParser


class SDFParser(BaseFileParser):
    """
    Parser for SDF files.
    """

    _allowed_formats = (".sdf",)

    def __init__(self, file_path: str):
        self._check_formats(file_path)
        super().__init__(file_path)
        self._parse()

    def _parse(self):
        for mol in pybel.readfile("sdf", self._file_path):
            self.append(
                SDFBlockParser(
                    mol.write("sdf"),
                    file_path=self._file_path,
                )
            )
