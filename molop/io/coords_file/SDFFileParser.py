"""
Author: TMJ
Date: 2024-06-20 22:47:13
LastEditors: TMJ
LastEditTime: 2024-06-20 22:59:19
Description: 请填写简介
"""
from openbabel import pybel

from molop.io.bases.BaseMolFileParser import BaseMolFileParser
from molop.io.coords_file.SDFFrameParser import SDFFrameParser


class SDFFileParser(BaseMolFileParser[SDFFrameParser]):
    _allowed_formats = (".sdf", ".mol")

    def _parse(self):
        for mol in pybel.readfile("sdf", self.file_path):
            if not self.only_last_frame:
                self.__append__(
                    SDFFrameParser(
                        frame_content=mol.write("sdf"),
                        file_path=self.file_path,
                    )
                )
        if self.only_last_frame:
            self.__append__(
                SDFFrameParser(
                    frame_content=mol.write("sdf"),
                    file_path=self.file_path,
                )
            )
