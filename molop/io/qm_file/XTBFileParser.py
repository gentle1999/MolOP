"""
Author: TMJ
Date: 2024-06-21 11:04:04
LastEditors: TMJ
LastEditTime: 2024-06-25 19:51:15
Description: 请填写简介
"""

import re
from typing import Literal

from packaging.version import Version
from pydantic import Field, computed_field

from molop.io.bases.BaseMolFileParser import BaseQMMolFileParser
from molop.io.qm_file.XTBFrameParser import XTBFrameParser
from molop.logger.logger import moloplogger
from molop.unit import atom_ureg
from molop.utils.xtbpatterns import xtboutpatterns


class XTBFileParser(BaseQMMolFileParser[XTBFrameParser]):
    _allowed_formats = (".out",)
    qm_software: str = Field(default="xTB")
    method: str = Field(default="SEMI-EMPIRICAL")

    def _parse(self):
        """
        Parse the file.
        """
        with open(self.file_path, "r") as fr:
            lines = fr.readlines()
        fr.close()
        full_text = "".join(lines)
        qm_software_version = self._parse_version(full_text)
        if self.qm_software_version is None:
            version = Version("0.0.0")
        else:
            version = Version(self.qm_software_version.split()[1])
        if version <= Version("6.2.2"):
            raise NotImplementedError("xtb <= 6.2.2 not supported.")
        charge = self._parse_charge(full_text)
        multi = self._parse_multi(full_text)
        self._parse_parameter(full_text)
        self._parse_functional(full_text)
        self._parse_solvent(full_text)
        self._parse_temperature(full_text)

        self.__append__(
            XTBFrameParser(
                frame_content=full_text,
                charge=charge,
                multiplicity=multi,
                qm_software_version=qm_software_version,
                method=self.method,
                basis=self.basis,
                functional=self.functional,
                solvent_model=self.solvent_model,
                solvent=self.solvent,
                temperature=self.temperature,
                file_path=self.file_path,
                keywords=self.keywords,
                only_extract_structure=self.only_extract_structure,
            )
        )
        self._parse_time()

    def _parse_temperature(self, full_text: str):
        temperature_match = re.search(xtboutpatterns["temperature"], full_text)
        if temperature_match:
            self.temperature = float(temperature_match.group(1)) * atom_ureg.K

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
            self.keywords = parameter_match.group(1)
        else:
            raise ValueError(
                f"No parameter comment found or illegal characters in {self.file_path}"
            )

    def _parse_multi(self, full_text: str):
        multi_match = re.search(xtboutpatterns["multiplicity"], full_text)
        if multi_match:
            multi = round(float(multi_match.group(2))) + 1
        else:
            multi = 1
        return multi

    def _parse_charge(self, full_text: str):
        charge_match = re.search(xtboutpatterns["charge"], full_text)
        if charge_match is None:
            charge_match = re.search(xtboutpatterns["total charge"], full_text)
        if charge_match:
            charge = round(float(charge_match.group(1)))
        else:
            charge = 0
        return charge

    def _parse_version(self, full_text: str):
        version_match = re.search(xtboutpatterns["version"], full_text)
        if version_match is None:
            version_match = re.search(xtboutpatterns["old version"], full_text)
        if version_match:
            version = version_match.group(1)
        else:
            version = "unknown"
        self.qm_software_version = version
        return version

    def _parse_time(self):
        self.running_time = sum(frame.running_time for frame in self.frames)

    @computed_field
    @property
    def task_type(self) -> Literal["sp", "opt", "freq"]:
        if "hess" in self.keywords.lower():
            return "freq"
        if "opt" in self.keywords.lower():
            return "opt"
        return "sp"
