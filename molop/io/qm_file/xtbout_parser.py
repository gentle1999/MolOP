"""
Author: TMJ
Date: 2024-02-17 15:17:37
LastEditors: TMJ
LastEditTime: 2024-03-25 22:49:52
Description: 请填写简介
"""

import os
import re

from molop.io.bases.file_base import BaseQMFileParser
from molop.io.qm_file.XTBOUTBlockParser import XTBOUTBlockParser
from molop.logger.logger import logger
from molop.utils.xtbpatterns import xtboutpatterns


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
        self._post_parse()

    def _parse(self):
        """
        Parse the file.
        """
        with open(self.file_path, "r") as fr:
            lines = fr.readlines()
        fr.close()
        full_text = "".join(lines)
        version = self._parse_version(full_text)
        charge = self._parse_charge(full_text)
        multi = self._parse_multi(full_text)
        self._parse_parameter(full_text)
        self._parse_functional(full_text)
        self._parse_solvent(full_text)
        self._parse_temperature(full_text)

        self.append(
            XTBOUTBlockParser(
                full_text,
                charge=charge,
                multiplicity=multi,
                version=version,
                basis=self.basis,
                functional=self.functional,
                solvent_model=self.solvent_model,
                solvent=self.solvent,
                temperature=self.temperature,
                file_path=self._file_path,
                parameter_comment=self.parameter_comment,
                only_extract_structure=self._only_extract_structure,
            )
        )

    def _parse_temperature(self, full_text: str):
        temperature_match = re.search(xtboutpatterns["temperature"], full_text)
        if temperature_match:
            self.temperature = float(temperature_match.group(1))

    def _parse_solvent(self, full_text: str):
        solvent_match = re.search(xtboutpatterns["solvent model"], full_text)
        if solvent_match:
            self.solvent_model = solvent_match.group(1)
            self.solvent = re.search(xtboutpatterns["solvent"], full_text).group(1)

    def _parse_functional(self, full_text: str):
        functional_match = re.search(xtboutpatterns["method"], full_text)
        if functional_match:
            self.functional = functional_match.group(1)
            self.basis = self.functional

    def _parse_parameter(self, full_text: str):
        parameter_match = re.search(xtboutpatterns["parameter"], full_text)
        if parameter_match:
            self.parameter_comment = parameter_match.group(1)
        else:
            logger.error(
                f"No parameter comment found or illegal characters in {self._file_path}"
            )
            raise ValueError(
                f"No parameter comment found or illegal characters in {self._file_path}"
            )

    def _parse_multi(self, full_text: str):
        multi_match = re.search(xtboutpatterns["multiplicity"], full_text)
        if multi_match:
            multi = int(float(multi_match.group(2))) + 1
        else:
            multi = 1
        if self.__force_multiplicity is not None:
            multi = self.__force_multiplicity
        return multi

    def _parse_charge(self, full_text: str):
        charge_match = re.search(xtboutpatterns["charge"], full_text)
        if charge_match is None:
            charge_match = re.search(xtboutpatterns["total charge"], full_text)
        if charge_match:
            charge = int(float(charge_match.group(1)))
        else:
            charge = 0
        if self.__force_charge is not None:
            charge = self.__force_charge
        return charge

    def _parse_version(self, full_text: str):
        version_match = re.search(xtboutpatterns["version"], full_text)
        if version_match is None:
            version_match = re.search(xtboutpatterns["old version"], full_text)
        if version_match:
            version = version_match.group(1)
        else:
            version = "unknown"
        self.version = version
        return version
