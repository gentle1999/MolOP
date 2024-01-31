"""
Author: TMJ
Date: 2023-10-30 15:41:19
LastEditors: TMJ
LastEditTime: 2024-01-31 22:24:15
Description: 请填写简介
"""
import os

from molop.io.bases.file_base import BaseFileParser
from molop.io.coords_file.XYZBlockParser import XYZBlockParser


class XYZParser(BaseFileParser):
    """
    Parser for XYZ files.
    Supports one or more molecules in one file.
    Make sure molecules in the file have same charge and multiplicity.
    """

    _allowed_formats = (".xyz",)

    def __init__(
        self, file_path: str, charge=None, multiplicity=None, only_last_frame=False
    ):
        self._check_formats(file_path)
        super().__init__(file_path)
        self._only_last_frame = only_last_frame
        self._charge = charge if charge else 0
        self._multiplicity = multiplicity if multiplicity else 1
        self._parse()

    def _parse(self):
        """
        Parse the file.
        """
        with open(self.file_path, "r") as fr:
            lines = fr.readlines()
        fr.close()
        anchor = 0
        blocks = []
        while anchor < len(lines):
            num_atoms = int(lines[anchor])
            blocks.append("".join(lines[anchor : anchor + num_atoms + 2]))
            anchor += num_atoms + 2
        if self._only_last_frame:
            self.append(
                XYZBlockParser(
                    blocks[-1],
                    charge=self._charge,
                    multiplicity=self._multiplicity,
                    file_path=self._file_path,
                )
            )
        else:
            for block in blocks:
                block = XYZBlockParser(
                    block,
                    charge=self._charge,
                    multiplicity=self._multiplicity,
                    file_path=self._file_path,
                )
                self.append(block)
