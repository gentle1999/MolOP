"""
Author: TMJ
Date: 2024-10-19 09:57:26
LastEditors: TMJ
LastEditTime: 2025-12-16 01:19:32
Description: 请填写简介
"""

import importlib.metadata

from molop.config import molopconfig
from molop.io import AutoParser
from molop.io.coords_parsers import (
    GJFFileParserDisk,
    GJFFileParserMemory,
    SDFFileParserDisk,
    SDFFileParserMemory,
    XYZFileParserDisk,
    XYZFileParserMemory,
)
from molop.io.QM_parsers import G16LogFileParserDisk, G16LogFileParserMemory


try:
    __version__ = importlib.metadata.version("molop")
except importlib.metadata.PackageNotFoundError:
    __version__ = "unknown"


__all__ = [
    "__version__",
    "AutoParser",
    "molopconfig",
    "GJFFileParserDisk",
    "GJFFileParserMemory",
    "SDFFileParserDisk",
    "SDFFileParserMemory",
    "XYZFileParserDisk",
    "XYZFileParserMemory",
    "G16LogFileParserDisk",
    "G16LogFileParserMemory",
]
