"""
Author: TMJ
Date: 2023-10-30 15:41:19
LastEditors: TMJ
LastEditTime: 2023-10-30 17:45:24
Description: 请填写简介
"""
import os
import re
from molop.io.base_parser import BaseParser


class XYZParser(BaseParser):
    """
    Parser for XYZ files.
    """

    def __init__(self, file_path):
        super().__init__(file_path)
        _, file_format = os.path.splitext(file_path)
        if file_format != ".xyz":
            raise ValueError("File format must be .xyz")
        self.parse()

    def parse(self):
        """
        Parse the file.
        """
        with open(self.file_path, "r") as fr:
            lines = fr.readlines()
        fr.close()
        num_atoms = int(lines[0])
        if re.search("Created by MolOP", lines[1]):
            extra_infos = lines[1].split("; ")[1:]
            for extra_info in extra_infos:
                prop, value = extra_info.split("=")
                self.__setattr__(f"_{prop}", eval(value))
        for line in lines[2 : 2 + num_atoms]:
            self.atom_coord_match(line)

    def atom_coord_match(self, line):
        if re.match(r"^\s*[A-Z][a-z]?(\s+\-?\d+(\.\d+)?){3}$", line):
            atom, x, y, z = line.split()
            self._atoms.append(atom)
            self._coords.extend([float(x), float(y), float(z)])
