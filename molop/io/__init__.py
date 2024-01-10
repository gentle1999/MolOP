'''
Author: TMJ
Date: 2023-10-30 15:40:03
LastEditors: TMJ
LastEditTime: 2024-01-10 19:07:34
Description: 请填写简介
'''
import os

from molop.io.bases.file_base import BaseFileParser
from molop.io.coords_file.xyz_parser import XYZParser
from molop.io.coords_file.sdf_parser import SDFParser
from molop.io.coords_file.gjf_parser import GJFParser
from molop.io.qm_file.g16log_parser import G16LOGParser


parsers = {
    ".gjf": GJFParser,
    ".log": G16LOGParser,
    ".xyz": XYZParser,
    ".sdf": SDFParser,
}


def AutoParser(file_path, charge=0, multiplicity=1) -> BaseFileParser:
    _, file_format = os.path.splitext(file_path)
    try:
        parser = parsers[file_format]
    except:
        raise NotImplementedError("Unknown file format: {}".format(file_format))
    if parser == XYZParser:
        return parser(file_path, charge=charge, multiplicity=multiplicity)
    else:
        return parser(file_path)
