"""
Author: TMJ
Date: 2025-09-10 23:03:06
LastEditors: TMJ
LastEditTime: 2025-09-12 13:38:07
Description: 请填写简介
"""

from typing import Dict, Tuple, Type

from molop.io.coords_models import (
    GJFFileDisk,
    GJFFileFrameDisk,
    GJFFileFrameMemory,
    GJFFileMemory,
    SDFFileDisk,
    SDFFileFrameDisk,
    SDFFileFrameMemory,
    SDFFileMemory,
    XYZFileDisk,
    XYZFileFrameDisk,
    XYZFileFrameMemory,
    XYZFileMemory,
)
from molop.io.coords_parsers import (
    GJFFileFrameParserDisk,
    GJFFileFrameParserMemory,
    GJFFileParserDisk,
    GJFFileParserMemory,
    SDFFileFrameParserDisk,
    SDFFileFrameParserMemory,
    SDFFileParserDisk,
    SDFFileParserMemory,
    XYZFileFrameParserDisk,
    XYZFileFrameParserMemory,
    XYZFileParserDisk,
    XYZFileParserMemory,
)
from molop.io.QM_models import (
    G16LogFileDisk,
    G16LogFileFrameDisk,
    G16LogFileFrameMemory,
    G16LogFileMemory,
)
from molop.io.QM_parsers import (
    G16LogFileFrameParserDisk,
    G16LogFileFrameParserMemory,
    G16LogFileParserDisk,
    G16LogFileParserMemory,
)

FILEDISK = G16LogFileDisk | GJFFileDisk | XYZFileDisk | SDFFileDisk
PARSERDISK = (
    GJFFileParserDisk | XYZFileParserDisk | SDFFileParserDisk | G16LogFileParserDisk
)
PARSERTYPEDISK = Type[PARSERDISK]
PARSERS_DICT: Dict[str, Tuple[PARSERTYPEDISK, ...]] = {
    ".gjf": (GJFFileParserDisk,),
    ".gau": (G16LogFileParserDisk,),
    ".com": (GJFFileParserDisk,),
    ".gjc": (GJFFileParserDisk,),
    ".log": (G16LogFileParserDisk,),
    ".g16": (G16LogFileParserDisk,),
    ".gal": (G16LogFileParserDisk,),
    ".xyz": (XYZFileParserDisk,),
    ".sdf": (SDFFileParserDisk,),
    ".mol": (SDFFileParserDisk,),
    ".out": (G16LogFileParserDisk,),
    ".irc": (G16LogFileParserDisk,),
}
FILEMEMORY = G16LogFileMemory | GJFFileMemory | XYZFileMemory | SDFFileMemory
PARSERMEMORY = (
    GJFFileParserMemory
    | XYZFileParserMemory
    | SDFFileParserMemory
    | G16LogFileParserMemory
)
PARSERTYPEMEMORY = Type[PARSERMEMORY]


__all__ = [
    "FILEDISK",
    "FILEMEMORY",
    "PARSERDISK",
    "PARSERMEMORY",
    "PARSERTYPEDISK",
    "PARSERTYPEMEMORY",
    "PARSERS_DICT",
    "G16LogFileFrameDisk",
    "G16LogFileFrameMemory",
    "GJFFileFrameDisk",
    "GJFFileFrameMemory",
    "SDFFileFrameDisk",
    "SDFFileFrameMemory",
    "XYZFileFrameDisk",
    "XYZFileFrameMemory",
    "GJFFileFrameParserDisk",
    "GJFFileFrameParserMemory",
    "SDFFileFrameParserDisk",
    "SDFFileFrameParserMemory",
    "XYZFileFrameParserDisk",
    "XYZFileFrameParserMemory",
    "G16LogFileFrameParserDisk",
    "G16LogFileFrameParserMemory",
]
