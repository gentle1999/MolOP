"""
Author: TMJ
Date: 2024-02-17 15:17:37
LastEditors: TMJ
LastEditTime: 2024-03-01 19:34:02
Description: 请填写简介
"""
import os
import re
from typing import Literal

from molop.utils import g16fchkpatterns, parameter_comment_parser

from molop.io.bases.file_base import BaseQMFileParser
from molop.io.qm_file.G16FCHKBlockParser import G16FCHKBlockParser
from molop.logger.logger import logger


class G16FCHKParser(BaseQMFileParser):
    _allowed_formats = (".fchk", ".fck", ".fch")

    def __init__(
        self,
        file_path: str,
        charge=None,
        multiplicity=None,
        only_extract_structure=False,
        only_last_frame=False,
    ):
        self._check_formats(file_path)
        super().__init__(file_path, only_extract_structure, only_last_frame)
        self.__force_charge = charge
        self.__force_multiplicity = multiplicity
        self._parse()

    def _parse(self):
        with open(self.file_path, "r") as fr:
            full_text = fr.read()
        fr.close()

        charge = (
            self.__force_charge
            if self.__force_charge
            else int(re.search(g16fchkpatterns["charge"], full_text).group(1))
        )
        multi = (
            self.__force_multiplicity
            if self.__force_multiplicity
            else int(re.search(g16fchkpatterns["multi"], full_text).group(1))
        )
        n_atoms = int(re.search(g16fchkpatterns["n_atoms"], full_text).group(1))
        self._version = re.search(g16fchkpatterns["version"], full_text).group(1)
        self._parameter_comment = re.search(g16fchkpatterns["route"], full_text).group(
            1
        )
        
        (
            _,
            self._route_params,
            self._dieze_tag,
            self._functional,
            self._basis_set,
        ) = parameter_comment_parser("\n"+self._parameter_comment)

        self.append(
            G16FCHKBlockParser(
                block=full_text,
                charge=charge,
                multiplicity=multi,
                n_atom=n_atoms,
                file_path=self._file_path,
                version=self._version,
                parameter_comment=self._parameter_comment,
                only_extract_structure=self._only_extract_structure,
            )
        )


    @property
    def route_params(self) -> dict:
        return self._route_params

    @property
    def dieze_tag(self) -> Literal["#N", "#P", "#T"]:
        return self._dieze_tag

    @property
    def functional(self) -> str:
        return self._functional

    @property
    def basis_set(self) -> str:
        return self._basis_set