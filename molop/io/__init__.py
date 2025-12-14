"""
Author: TMJ
Date: 2025-01-15 23:01:22
LastEditors: TMJ
LastEditTime: 2025-12-14 15:38:27
Description: 请填写简介
"""

import glob
import os
from pathlib import Path
from typing import Literal

from molop.io.FileBatchModelDisk import FileBatchModelDisk
from molop.io.FileBatchParserDisk import FileBatchParserDisk
from molop.io.types import (
    G16LogFileParserDisk,
    G16LogFileParserMemory,
    GJFFileParserDisk,
    GJFFileParserMemory,
    SDFFileParserDisk,
    SDFFileParserMemory,
    XYZFileParserDisk,
    XYZFileParserMemory,
)

__all__ = [
    "AutoParser",
    "GJFFileParserMemory",
    "XYZFileParserMemory",
    "SDFFileParserMemory",
    "G16LogFileParserMemory",
    "G16LogFileParserDisk",
    "XYZFileParserDisk",
    "SDFFileParserDisk",
    "GJFFileParserDisk",
]


def split_path_pattern(path_str: str) -> tuple[Path, str]:
    p = Path(path_str)
    parts = p.parts

    split_index = len(parts)
    for i, part in enumerate(parts):
        if glob.has_magic(part):
            split_index = i
            break
    base_path = Path(*parts[:split_index])
    if split_index < len(parts):
        pattern = str(Path(*parts[split_index:]))
    else:
        pattern = ""

    return base_path, pattern


def AutoParser(
    file_path: str,
    *,
    total_charge: int | None = None,
    total_multiplicity: int | None = None,
    n_jobs: int = -1,
    only_extract_structure=False,
    only_last_frame=False,
    release_file_content: bool = True,
    parser_detection: Literal["auto", "gjf", "xyz", "sdf", "g16log"] = "auto",
) -> FileBatchModelDisk:
    """
    The Entrypoint of MolOP

    Parameters:
        file_path (str):
            use regax to match files.
        total_charge (int | None):
            forced charge of the molecule, if not given, will use the charge written in the file or 0.
        total_multiplicity (int | None):
            forced multiplicity of the molecule, if not given, will use the charge written in the file or 1.
        n_jobs (int):
            number of jobs to use, if -1, use cpu with default max number.
        only_extract_structure (bool):
            if True, only extract the structure, else extract the whole file.
        only_last_frame (bool):
            if True, only extract the last frame, else extract all frames.
        release_file_content (bool):
            if True, release the file content after parsing, else keep the file content in memory.
        parser_detection (Literal["auto", "gjf", "xyz", "sdf", "g16log"]):
            if "auto", use the file extension to detect the parser, else use the given parser.

    Returns:
        FileBatchParserDisk

    How to use
    ----------
    ```python
    from molop import AutoParser
    parser = AutoParser("/path/to/file") # return FileBatchParserDisk with singlefile_parser
    parser = AutoParser("/path/to/files/*.log") # return FileBatchParserDisk
    ```
    """
    if os.path.isfile(file_path) and os.path.exists(file_path):
        files = [file_path]
    else:
        base_path, pattern = split_path_pattern(file_path)
        files = list(base_path.glob(pattern))
    return FileBatchParserDisk(n_jobs=n_jobs).parse(
        files,
        total_charge=total_charge,
        total_multiplicity=total_multiplicity,
        only_extract_structure=only_extract_structure,
        only_last_frame=only_last_frame,
        release_file_content=release_file_content,
        parser_detection=parser_detection,
    )
