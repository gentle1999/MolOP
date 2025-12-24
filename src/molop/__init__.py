"""
Author: TMJ
Date: 2024-10-19 09:57:26
LastEditors: TMJ
LastEditTime: 2025-12-22 13:08:29
Description: MolOP is a toolbox for molecule operations and QM information extraction.
"""

import importlib.metadata

from molop.config import molopconfig, moloplogger
from molop.io import AutoParser
from molop.io.base_models.Molecule import Molecule
from molop.io.coords_parsers import (
    GJFFileParserDisk,
    GJFFileParserMemory,
    SDFFileParserDisk,
    SDFFileParserMemory,
    SMIFileParserDisk,
    SMIFileParserMemory,
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
    "moloplogger",
    "GJFFileParserDisk",
    "GJFFileParserMemory",
    "SDFFileParserDisk",
    "SDFFileParserMemory",
    "XYZFileParserDisk",
    "XYZFileParserMemory",
    "G16LogFileParserDisk",
    "G16LogFileParserMemory",
    "SMIFileParserDisk",
    "SMIFileParserMemory",
    "Molecule",
]
