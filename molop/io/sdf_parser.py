'''
Author: TMJ
Date: 2023-10-30 18:21:31
LastEditors: TMJ
LastEditTime: 2023-10-31 09:19:28
Description: 请填写简介
'''

import os
import re
from molop.io.base_parser import BaseParser

class SDFParser(BaseParser):
    """
    Parser for SDF files.
    """

    def __init__(self, file_path):
        super().__init__(file_path)
        _, file_format = os.path.splitext(file_path)
        if file_format != ".sdf":
            raise ValueError("File format must be .sdf")
        self.parse()

    def parse(self):
        pass