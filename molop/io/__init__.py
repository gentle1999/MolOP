"""
Author: TMJ
Date: 2023-10-30 15:40:03
LastEditors: TMJ
LastEditTime: 2024-01-07 14:02:21
Description: 请填写简介
"""
import os

from molop.io.gaussian_parser import Gaussian16GJFParser, Gaussian16LOGParser
from molop.io.xyz_parser import XYZParser

parsers = {
    ".gjf": Gaussian16GJFParser,
    ".log": Gaussian16LOGParser,
    ".xyz": XYZParser,
}


def AutoParser(file_path):
    _, file_format = os.path.splitext(file_path)
    try:
        parser = parsers[file_format]
    except:
        raise NotImplementedError("Unknown file format: {}".format(file_format))
    return parser(file_path)
