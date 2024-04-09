'''
Author: TMJ
Date: 2024-02-17 15:17:37
LastEditors: TMJ
LastEditTime: 2024-04-09 10:25:21
Description: 请填写简介
'''

import os

from openbabel import pybel

from molop.io.bases.file_base import BaseFileParser
from molop.io.coords_file.SDFBlockParser import SDFBlockParser


class SDFParser(BaseFileParser):
    """
    Parser for SDF files.
    """

    _allowed_formats = (".sdf", ".mol")

    def __init__(self, file_path: str, only_last_frame=False):
        self._check_formats(file_path)
        super().__init__(file_path)
        self._only_last_frame = only_last_frame
        self._parse()
        self._post_parse()

    def _parse(self):
        for mol in pybel.readfile("sdf", self._file_path):
            if not self._only_last_frame:
                self.append(
                    SDFBlockParser(
                        mol.write("sdf"),
                        file_path=self._file_path,
                    )
                )
        if self._only_last_frame:
            self.append(
                SDFBlockParser(
                    mol.write("sdf"),
                    file_path=self._file_path,
                )
            )
