"""
Author: TMJ
Date: 2024-10-19 09:57:26
LastEditors: TMJ
LastEditTime: 2024-10-23 15:01:26
Description: 请填写简介
"""

import os
from glob import glob
from typing import Union

from molop.io.file_batch import FileParserBatch


def AutoParser(
    file_path: str,
    charge: int = 0,
    multiplicity: int = 1,
    n_jobs: int = -1,
    only_extract_structure=False,
    only_last_frame=False,
) -> FileParserBatch:
    """
    The Entrypoint of MolOP

    Parameters:
        file_path str:
            use regax to match files.
        charge int:
            forced charge of the molecule, if not given, will use the charge written in the file or 0.
        multiplicity int:
            forced multiplicity of the molecule, if not given, will use the charge written in the file or 1.
        n_jobs int:
            number of jobs to use, if -1, use all cpu.
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
    if os.path.isfile(file_path) and os.path.exists(file_path):
        files = [file_path]
    else:
        files = glob(file_path)
        files.sort()
    if len(files) > 0:
        batch = FileParserBatch(
            n_jobs=n_jobs,
        )
        batch.add_files(
            files,
            charge=charge,
            multiplicity=multiplicity,
            only_extract_structure=only_extract_structure,
            only_last_frame=only_last_frame,
        )
        return batch
    else:
        raise FileNotFoundError("No file found in the path")
