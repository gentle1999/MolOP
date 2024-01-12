"""
Author: TMJ
Date: 2023-10-30 15:40:03
LastEditors: TMJ
LastEditTime: 2024-01-10 20:39:40
Description: 请填写简介
"""
import os
from typing import Union

from molop.io.bases.file_base import BaseFileParser
from molop.io.coords_file.gjf_parser import GJFParser
from molop.io.coords_file.sdf_parser import SDFParser
from molop.io.coords_file.xyz_parser import XYZParser
from molop.io.qm_file.g16irc_parser import G16IRCParser
from molop.io.qm_file.g16log_parser import G16LOGParser
from molop.io.qm_file.xtbout_parser import XTBOUTParser

from molop.logger.logger import logger

parsers = {
    ".gjf": (GJFParser,),
    ".log": (G16LOGParser,),
    ".xyz": (XYZParser,),
    ".sdf": (SDFParser,),
    ".out": (XTBOUTParser, G16IRCParser),
    ".irc": (G16IRCParser,),
}


def AutoParser(
    file_path, charge=None, multiplicity=None, only_extract_structure=False
) -> Union[BaseFileParser, XYZParser, SDFParser, GJFParser, G16LOGParser, XTBOUTParser]:
    _, file_format = os.path.splitext(file_path)
    if file_format not in parsers:
        raise NotImplementedError("Unknown file format: {}".format(file_format))
    for idx, parser in enumerate(parsers[file_format]):
        try:
            if parser in (G16LOGParser, XTBOUTParser):
                return parser(
                    file_path,
                    charge=charge,
                    multiplicity=multiplicity,
                    only_extract_structure=only_extract_structure,
                )
            elif parser in (XYZParser, GJFParser):
                return parser(file_path, charge=charge, multiplicity=multiplicity)
            else:
                return parser(file_path)
        except:
            if idx == len(parsers[file_format]) - 1:
                logger.error(f"Failed to parse file {file_path} with {parser.__name__}")
                raise Exception(f"Failed to parse file {file_path}")
            logger.warning(
                f"Failed to parse file {file_path} with {parser.__name__}, try {parsers[file_format][idx+1].__name__} instead"
            )
