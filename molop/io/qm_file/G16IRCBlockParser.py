"""
Author: TMJ
Date: 2024-01-13 12:20:19
LastEditors: TMJ
LastEditTime: 2024-01-13 12:20:29
Description: 请填写简介
"""
from molop.io.qm_file.G16LOGBlockParser import G16LOGBlockParser


class G16IRCBlockParser(G16LOGBlockParser):
    _block_type = "G16 IRC"

    def __init__(
        self,
        block: str,
        charge=0,
        multiplicity=1,
        n_atom=1,
        file_path="",
        version=None,
        parameter_comment=None,
        only_extract_structure=False,
    ):
        super().__init__(
            block,
            charge,
            multiplicity,
            n_atom,
            file_path,
            version,
            parameter_comment,
            only_extract_structure,
        )
