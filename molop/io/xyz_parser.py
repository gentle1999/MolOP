"""
Author: TMJ
Date: 2023-10-30 15:41:19
LastEditors: TMJ
LastEditTime: 2023-11-02 10:21:19
Description: 请填写简介
"""
import os
import re
from molop.io.base_parser import BaseParser, MultiFrameBaseParser


class XYZBlockParser(BaseParser):
    """
    Parser for XYZ Blocks.
    """

    def __init__(self, block):
        super().__init__(block)
        self._parse()

    def _parse(self):
        """
        Parse the file.
        """
        lines = self._block.split("\n")
        num_atoms = int(lines[0])
        if re.search("Created by MolOP", lines[1]):
            extra_infos = lines[1].split("; ")[1:]
            for extra_info in extra_infos:
                prop, value = extra_info.split("==")
                self.__setattr__(f"_{prop}_attach_info", value)
        for line in lines[2 : 2 + num_atoms]:
            self._atom_coord_match(line)

    def _atom_coord_match(self, line):
        if re.match(r"^\s*[A-Z][a-z]?(\s+\-?\d+(\.\d+)?){3}$", line):
            atom, x, y, z = line.split()
            self._atoms_attach_info.append(atom)
            self._coords_attach_info.extend([float(x), float(y), float(z)])


class XYZParser(MultiFrameBaseParser):
    """
    Parser for XYZ files.
    """

    def __init__(self, file_path):
        super().__init__(file_path)
        _, file_format = os.path.splitext(file_path)
        if file_format != ".xyz":
            raise ValueError("File format must be .xyz")
        self._parse_meta()
        self._split_frames()
        self._parse()

    def _parse_meta(self):
        pass

    def _split_frames(self):
        with open(self.file_path, "r") as fr:
            lines = fr.readlines()
        fr.close()
        anchor = 0
        while anchor < len(lines):
            num_atoms = int(lines[anchor])
            self._blocks.append("".join(lines[anchor : anchor + num_atoms + 2]))
            anchor += num_atoms + 2

    def _parse(self):
        """
        Parse the file.
        """
        for block in self._blocks:
            frame = XYZBlockParser(block)
            self._append(frame)
