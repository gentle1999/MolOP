"""
Author: TMJ
Date: 2023-10-30 15:38:58
LastEditors: TMJ
LastEditTime: 2023-11-01 21:49:27
Description: 请填写简介
"""
import os
import re

import parse

from molop.io.base_parser import BaseParser, MultiFrameBaseParser


class Gaussian16GJFBlockParser(BaseParser):
    """
    Parser for Gaussian16 gjf files.
    """

    def __init__(self, block):
        super().__init__(block)
        self._parse()

    def _parse(self):
        """
        Parse the file.
        """
        lines = self._block.split("\n")
        for line in lines:
            self._charge_multi_match(line)
            self._atom_coord_match(line)

    def _charge_multi_match(self, line: str):
        if re.match(r"^\s*\d+\s+\d+$", line):
            charge, multi = map(int, line.strip().split())
            self._charge_attach_info, self._multi_attach_info = charge, multi

    def _atom_coord_match(self, line: str):
        if re.match(r"^\s*[A-Z][a-z]?(\s+\-?\d+(\.\d+)?){3}$", line):
            atom, x, y, z = line.split()
            self._atoms_attach_info.append(atom)
            self._coords_attach_info.extend([float(x), float(y), float(z)])

    def _to_format_coords(self) -> str:
        block = f"{self.charge} {self.multi}\n"
        coords = [self.coords[i : i + 3] for i in range(0, len(self.coords), 3)]
        for atom, coord in zip(self.atoms, coords):
            block += f"{atom:3} {coord[0]:>15.8f} {coord[1]:>15.8f} {coord[2]:>15.8f}\n"
        return block

    def to_GJF_block(self):
        block = "".join(
            [
                f"%{attr[4:-12]}={self.__getattribute__(attr)}\n"
                for attr in self._get_attach_infos()
                if "G16" in attr
            ]
        )
        block += f"#{self._RouteSection_attach_info}\n{self._PreSetting}"
        block += self._to_format_coords()
        block += self._PostSetting
        return block


class Gaussian16GJFParser(MultiFrameBaseParser):
    """
    Parser for Gaussian16 gjf files.
    """

    def __init__(self, file_path):
        super().__init__(file_path)
        _, file_format = os.path.splitext(file_path)
        if file_format != ".gjf":
            raise ValueError("File format must be .gjf")
        self._parse_meta()
        self._split_frames()
        self._parse()

    def _parse_meta(self):
        with open(self.file_path, "r") as fr:
            lines = fr.readlines()
        fr.close()
        for i, line in enumerate(lines):
            if re.match(r"%[^=]+", line):
                attr, value = line.strip().split("%")[1].split("=")
                self.__setattr__(f"_G16{attr}_attach_info", value)
            if re.match(r"^\s*\#", line):
                self._RouteSection_attach_info = line.strip().split("#")[1]
                pre_setting_start = i
            if re.match(r"^\s*\d+\s+\d+$", line):
                self._block_start = i
            if re.match(r"^\s*[A-Z][a-z]?(\s+\-?\d+(\.\d+)?){3}$", line):
                self._block_end = i
        self._PreSetting = "".join(lines[pre_setting_start + 1 : self._block_start])
        self._PostSetting = "".join(lines[self._block_end + 1 :])

    def _split_frames(self):
        with open(self.file_path, "r") as fr:
            lines = fr.readlines()
        fr.close()
        self._blocks.append("".join(lines[self._block_start : self._block_end + 1]))

    def _parse(self):
        """
        Parse the file.
        """
        for block in self._blocks:
            frame = Gaussian16GJFBlockParser(block)
            frame._PreSetting = self._PreSetting
            frame._PostSetting = self._PostSetting
            self._append(frame)

    def to_GJF_block(self):
        return self[0].to_GJF_block()


class Gaussian16LOGParser(BaseParser):
    """
    Parser for Gaussian16 log files.
    """

    def __init__(self, file_path):
        super().__init__(file_path)
        _, file_format = os.path.splitext(file_path)
        if file_format != ".log":
            raise ValueError("File format must be .log")
        self._parse()

    def _parse(self):
        """
        Parse the file.
        """
