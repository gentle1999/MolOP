"""
Author: TMJ
Date: 2024-06-20 20:52:50
LastEditors: TMJ
LastEditTime: 2024-08-04 21:01:32
Description: 请填写简介
"""

from typing import Tuple

from molop.io.bases.BaseMolFileParser import BaseMolFileParser
from molop.io.coords_file.XYZFrameParser import XYZFrameParser


class XYZFileParser(BaseMolFileParser[XYZFrameParser]):
    """
    XYZ file parser
    """

    _allowed_formats: Tuple[str] = (".xyz",)

    def _parse(self):
        with open(self.file_path, "r") as fr:
            lines = fr.readlines()
        fr.close()
        anchor = 0
        frames = []
        while anchor < len(lines):
            try:
                num_atoms = int(lines[anchor])
            except:
                anchor += 1
                continue
            frames.append("".join(lines[anchor : anchor + num_atoms + 2]))
            anchor += num_atoms + 2
        if self.only_last_frame:
            self.__append__(
                XYZFrameParser(
                    frame_content=frames[-1],
                    charge=self.charge,
                    multiplicity=self.multiplicity,
                    file_path=self.file_path,
                )
            )
        else:
            for frame in frames:
                self.__append__(
                    XYZFrameParser(
                        frame_content=frame,
                        charge=self.charge,
                        multiplicity=self.multiplicity,
                        file_path=self.file_path,
                    )
                )
        if len(frames) > 0:
            self.charge = self.frames[-1].charge
            self.multiplicity = self.frames[-1].multiplicity
