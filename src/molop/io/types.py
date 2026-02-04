"""
Author: TMJ
Date: 2025-09-10 23:03:06
LastEditors: TMJ
LastEditTime: 2026-02-04 15:25:27
Description: 请填写简介
"""

from typing import TypeAlias

from molop.io.coords_models import (
    GJFFileDisk,
    GJFFileFrameDisk,
    GJFFileFrameMemory,
    GJFFileMemory,
    SDFFileDisk,
    SDFFileFrameDisk,
    SDFFileFrameMemory,
    SDFFileMemory,
    SMIFileDisk,
    SMIFileFrameDisk,
    SMIFileFrameMemory,
    SMIFileMemory,
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
    SMIFileFrameParserDisk,
    SMIFileFrameParserMemory,
    SMIFileParserDisk,
    SMIFileParserMemory,
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


FILEDISK: TypeAlias = G16LogFileDisk | GJFFileDisk | XYZFileDisk | SDFFileDisk | SMIFileDisk
PARSERDISK: TypeAlias = (
    GJFFileParserDisk
    | XYZFileParserDisk
    | SDFFileParserDisk
    | G16LogFileParserDisk
    | SMIFileParserDisk
)
PARSERTYPEDISK: TypeAlias = (
    type[GJFFileParserDisk]
    | type[XYZFileParserDisk]
    | type[SDFFileParserDisk]
    | type[G16LogFileParserDisk]
    | type[SMIFileParserDisk]
)
PARSERS_DICT: dict[str, tuple[PARSERTYPEDISK, ...]] = {
    ".gjf": (GJFFileParserDisk,),
    ".gau": (G16LogFileParserDisk,),
    ".com": (GJFFileParserDisk,),
    ".gjc": (GJFFileParserDisk,),
    ".log": (G16LogFileParserDisk,),
    ".g16": (G16LogFileParserDisk,),
    ".gal": (G16LogFileParserDisk,),
    ".xyz": (XYZFileParserDisk,),
    ".smi": (SMIFileParserDisk,),
    ".txt": (SMIFileParserDisk,),
    ".sdf": (SDFFileParserDisk,),
    ".mol": (SDFFileParserDisk,),
    ".out": (G16LogFileParserDisk,),
    ".irc": (G16LogFileParserDisk,),
}
FILEMEMORY: TypeAlias = (
    G16LogFileMemory | GJFFileMemory | XYZFileMemory | SDFFileMemory | SMIFileMemory
)
PARSERMEMORY: TypeAlias = (
    GJFFileParserMemory
    | XYZFileParserMemory
    | SDFFileParserMemory
    | G16LogFileParserMemory
    | SMIFileParserMemory
)
PARSERTYPEMEMORY: TypeAlias = (
    type[GJFFileParserMemory]
    | type[XYZFileParserMemory]
    | type[SDFFileParserMemory]
    | type[G16LogFileParserMemory]
    | type[SMIFileParserMemory]
)


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
    "SMIFileFrameParserDisk",
    "SMIFileFrameParserMemory",
    "SMIFileParserDisk",
    "SMIFileParserMemory",
    "SMIFileFrameDisk",
    "SMIFileFrameMemory",
]
