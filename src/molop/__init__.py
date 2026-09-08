"""
Author: TMJ
Date: 2024-10-19 09:57:26
LastEditors: TMJ
LastEditTime: 2026-02-05 16:37:11
Description: MolOP is a toolbox for molecule operations and QM information extraction.
"""

import importlib.metadata

from molop._runtime import configure_native_thread_limits


configure_native_thread_limits()

# Load RDKit's native Chem extension before MolOP imports Open Babel.
from rdkit import Chem as _rdkit_chem  # noqa: E402, F401

from molop.config import molopconfig, moloplogger  # noqa: E402
from molop.io import (  # noqa: E402
    AutoBytesParser,
    AutoFileParser,
    AutoMemoryParser,
    AutoParser,
    AutoParserMemory,
    AutoTextParser,
)
from molop.io.codec_types import ParseOptions  # noqa: E402
from molop.io.parse_outcomes import (  # noqa: E402
    BatchParseResult,
    FileParseOutcome,
    ParseFailure,
)


try:
    __version__ = importlib.metadata.version("molop")
except importlib.metadata.PackageNotFoundError:
    __version__ = "unknown"


__all__ = [
    "__version__",
    "AutoFileParser",
    "AutoBytesParser",
    "AutoMemoryParser",
    "AutoParser",
    "AutoParserMemory",
    "AutoTextParser",
    "BatchParseResult",
    "FileParseOutcome",
    "ParseFailure",
    "ParseOptions",
    "molopconfig",
    "moloplogger",
]
