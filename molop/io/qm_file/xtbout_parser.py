"""
Author: TMJ
Date: 2024-02-17 15:17:37
LastEditors: TMJ
LastEditTime: 2024-03-02 21:17:14
Description: 请填写简介
"""
import os
import re

from molop.io.bases.file_base import BaseQMFileParser
from molop.io.qm_file.XTBOUTBlockParser import XTBOUTBlockParser
from molop.logger.logger import logger
from molop.utils import xtboutpatterns


class XTBOUTParser(BaseQMFileParser):
    _allowed_formats = (".out",)

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
        """
        Parse the file.
        """
        with open(self.file_path, "r") as fr:
            lines = fr.readlines()
        fr.close()
        full_text = "".join(lines)
        version_match = re.search(xtboutpatterns["version"], full_text)
        if version_match is None:
            version_match = re.search(xtboutpatterns["old version"], full_text)
        if version_match:
            version = version_match.group(1)
        else:
            version = ""
        self._version = version
        charge_match = re.search(xtboutpatterns["charge"], full_text)
        if charge_match is None:
            charge_match = re.search(xtboutpatterns["total charge"], full_text)
        if charge_match:
            charge = int(float(charge_match.group(1)))
        if self.__force_charge is not None:
            charge = self.__force_charge
        multi = 1
        if self.__force_multiplicity is not None:
            multi = self.__force_multiplicity
        parameter_match = re.search(xtboutpatterns["parameter"], full_text)
        if parameter_match:
            self._parameter_comment = parameter_match.group(1)
        else:
            logger.error(
                f"No parameter comment found or illegal characters in {self._file_path}"
            )
            raise ValueError(
                f"No parameter comment found or illegal characters in {self._file_path}"
            )

        self.append(
            XTBOUTBlockParser(
                full_text,
                charge=charge,
                multiplicity=multi,
                version=version,
                file_path=self._file_path,
                parameter_comment=self._parameter_comment,
                only_extract_structure=self._only_extract_structure,
            )
        )
