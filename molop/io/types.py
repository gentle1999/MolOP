"""
Author: TMJ
Date: 2024-02-25 14:51:27
LastEditors: TMJ
LastEditTime: 2024-03-18 21:27:47
Description: 请填写简介
"""
from typing import Literal, Union

from molop.io.bases.file_base import BLOCKTYPES, QMBLOCKTYPES
from molop.io.coords_file.gjf_parser import GJFParser
from molop.io.coords_file.sdf_parser import SDFParser
from molop.io.coords_file.xyz_parser import XYZParser
from molop.io.qm_file.g16fchk_parser import G16FCHKParser
from molop.io.qm_file.g16log_parser import G16LOGParser
from molop.io.qm_file.xtbout_parser import XTBOUTParser

PARSERTYPES = Union[
    GJFParser,
    SDFParser,
    XYZParser,
    G16FCHKParser,
    G16LOGParser,
    XTBOUTParser,
]
