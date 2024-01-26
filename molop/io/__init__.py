"""
Author: TMJ
Date: 2023-10-30 15:40:03
LastEditors: TMJ
LastEditTime: 2024-01-24 16:29:47
Description: 请填写简介
"""
import os
from typing import Union

from molop.io.file_batch import singlefile_parser, FileParserBatch


def AutoParser(
    file_path,
    charge=None,
    multiplicity=None,
    only_extract_structure=False,
    only_last_frame=False,
):
    if os.path.isfile(file_path):
        return singlefile_parser(
            file_path,
            charge=charge,
            multiplicity=multiplicity,
            only_extract_structure=only_extract_structure,
            only_last_frame=only_last_frame,
        )
    elif os.path.isdir(file_path):
        return FileParserBatch(
            file_path,
            charge=charge,
            multiplicity=multiplicity,
            only_extract_structure=only_extract_structure,
            only_last_frame=only_last_frame,
        )
