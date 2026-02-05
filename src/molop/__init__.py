"""
Author: TMJ
Date: 2024-10-19 09:57:26
LastEditors: TMJ
LastEditTime: 2026-02-05 16:37:11
Description: MolOP is a toolbox for molecule operations and QM information extraction.
"""

import importlib.metadata

from molop.config import molopconfig, moloplogger
from molop.io import AutoParser


try:
    __version__ = importlib.metadata.version("molop")
except importlib.metadata.PackageNotFoundError:
    __version__ = "unknown"


__all__ = [
    "__version__",
    "AutoParser",
    "molopconfig",
    "moloplogger",
]
