"""
Author: TMJ
Date: 2024-01-24 12:33:24
LastEditors: TMJ
LastEditTime: 2024-01-30 17:42:17
Description: 请填写简介
"""
import os
import re

from molop.io.bases.file_base import BaseQMFileParser
from molop.io.qm_file.G16FCHKBlockParser import G16FCHKBlockParser
from molop.logger.logger import logger


class G16FCHKParser(BaseQMFileParser):
    _allowed_formats = (".fchk", ".fck")

    def __init__(
        self,
        file_path: str,
        charge=None,
        multiplicity=None,
        show_progress=False,
        only_extract_structure=False,
        only_last_frame=False,
    ):
        self._check_formats(file_path)
        super().__init__(
            file_path, show_progress, only_extract_structure, only_last_frame
        )
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
            else int(re.findall(r"Charge\s+[A-Z]+\s+([\-0-9]+)", full_text)[0])
        )
        multi = (
            self.__force_multiplicity
            if self.__force_multiplicity
            else int(re.findall(r"Multiplicity\s+[A-Z]+\s+([\-0-9]+)", full_text)[0])
        )
        n_atoms = int(re.findall(r"Number of atoms\s+[A-Z]+\s+([0-9]+)", full_text)[0])
        self._version = re.findall(
            r"Gaussian Version\s+[A-Z\s+\=]+\s+[\-0-9]+\s+([a-zA-Z0-9\-\.]+)", full_text
        )[0]
        self._parameter_comment = (
            "#"
            + re.findall(
                r"\#([a-zA-Z%0-9\.\#\=\s\-\+\(\)\,\"\*\/\\^]+)\nCharge", full_text
            )[0]
        )
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
