'''
Author: TMJ
Date: 2023-10-30 15:38:58
LastEditors: TMJ
LastEditTime: 2023-10-30 17:31:47
Description: 请填写简介
'''
"""
Author: TMJ
Date: 2023-10-30 15:38:58
LastEditors: TMJ
LastEditTime: 2023-10-30 15:42:58
Description: 请填写简介
"""
import os
import re
from molop.io.base_parser import BaseParser


class Gaussian16GJFParser(BaseParser):
    """
    Parser for Gaussian16 gjf files.
    """

    def __init__(self, file_path):
        super().__init__(file_path)
        _, file_format = os.path.splitext(file_path)
        if file_format != ".gjf":
            raise ValueError("File format must be .gjf")
        self.parse()

    def parse(self):
        """
        Parse the file.
        """
        with open(self.file_path, "r") as fr:
            lines = fr.readlines()
        fr.close()
        for line in lines:
            self.charge_multi_match(line)
            self.atom_coord_match(line)

    def charge_multi_match(self, line):
        if re.match(r"^\s*\d+\s+\d+$", line):
            charge, multi = map(int, line.split())
            self._charge, self._multi = charge, multi

    def atom_coord_match(self, line):
        if re.match(r"^\s*[A-Z][a-z]?(\s+\-?\d+(\.\d+)?){3}$", line):
            atom, x, y, z = line.split()
            self._atoms.append(atom)
            self._coords.extend([float(x), float(y), float(z)])

    #TODO gjf文件中还有很多其他信息，例如基组，计算参数等


class Gaussian16LOGParser(BaseParser):
    """
    Parser for Gaussian16 log files.
    """

    def __init__(self, file_path):
        super().__init__(file_path)
        _, file_format = os.path.splitext(file_path)
        if file_format != ".log":
            raise ValueError("File format must be .log")
        self.parse()

    def parse(self):
        """
        Parse the file.
        """
