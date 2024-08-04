"""
Author: TMJ
Date: 2024-02-25 14:51:27
LastEditors: TMJ
LastEditTime: 2024-08-04 21:02:19
Description: 请填写简介
"""

from typing import TypeVar, Union

from molop.io.bases.BaseMolFileParser import BaseMolFileParser, BaseQMMolFileParser
from molop.io.coords_file.GJFFileParser import GJFFileParser
from molop.io.coords_file.SDFFileParser import SDFFileParser
from molop.io.coords_file.XYZFileParser import XYZFileParser
from molop.io.qm_file.G16FchkFileParser import G16FchkFileParser
from molop.io.qm_file.G16LogFileParser import G16LogFileParser
from molop.io.qm_file.XTBFileParser import XTBFileParser

MolFileParserType = TypeVar(
    "MolFileParserType", bound=BaseMolFileParser, covariant=True
)
QMMolFileParserType = TypeVar(
    "QMMolFileParserType", bound=BaseQMMolFileParser, covariant=True
)
PARSERTYPES = Union[
    GJFFileParser,
    XYZFileParser,
    SDFFileParser,
    G16LogFileParser,
    G16FchkFileParser,
    XTBFileParser,
]
