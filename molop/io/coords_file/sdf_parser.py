'''
Author: TMJ
Date: 2023-10-30 18:21:31
LastEditors: TMJ
LastEditTime: 2024-01-11 09:58:19
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

    def __init__(self, file_path: str):
        super().__init__(file_path)
        _, file_format = os.path.splitext(file_path)
        if file_format != ".sdf":
            raise ValueError("File format must be .sdf")
        self._parse()

    def _parse(self):
        for mol in pybel.readfile("sdf", self._file_path):
            self.append(SDFBlockParser(mol.write("sdf")))
