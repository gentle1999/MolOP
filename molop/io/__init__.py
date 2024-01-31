"""
Author: TMJ
Date: 2023-10-30 15:40:03
LastEditors: TMJ
LastEditTime: 2024-01-30 16:18:40
Description: The IO unit of MolOP
"""
from glob import glob
from typing import Union

from molop.io.file_batch import PARSERTYPES, FileParserBatch, singlefile_parser


def AutoParser(
    file_path: str,
    charge: int = None,
    multiplicity: int = None,
    only_extract_structure=False,
    only_last_frame=False,
) -> Union[FileParserBatch, None]:
    """
    The Entrypoint of MolOP

    Parameters:
        file_path str:
            use regax to match files.
        charge int:
            forced charge of the molecule, if not given, will use the charge written in the file or 0.
        multiplicity int:
            forced multiplicity of the molecule, if not given, will use the charge written in the file or 1.
        only_extract_structure bool:
            if True, only extract the structure, else extract the whole file.
        only_last_frame bool:
            if True, only extract the last frame, else extract all frames.

    Returns:
        FileParserBatch

    How to use
    ----------
    ```python
    from molop import AutoParser
    parser = AutoParser("/path/to/file") # return FileParserBatch with singlefile_parser
    parser = AutoParser("/path/to/files/*.log") # return FileParserBatch
    ```
    """
    files = glob(file_path)
    files.sort()
    if len(files) > 0:
        return FileParserBatch(
            files,
            charge=charge,
            multiplicity=multiplicity,
            only_extract_structure=only_extract_structure,
            only_last_frame=only_last_frame,
        )
    else:
        raise FileNotFoundError("No file found in the path")
