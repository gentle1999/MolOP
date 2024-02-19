"""
Author: TMJ
Date: 2024-02-17 15:17:37
LastEditors: TMJ
LastEditTime: 2024-02-17 20:03:06
Description: 请填写简介
"""
from typing import Literal
import re
import time
from molop.io.bases.file_base import BaseQMFileParser
from molop.io.qm_file.G16LOGBlockParser import (
    G16LOGBlockParser,
    parameter_comment_parser,
)
from molop.logger.logger import logger


class G16LOGParser(BaseQMFileParser):
    _allowed_formats = (".log", ".g16", ".gal")

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
        self._parse(charge, multiplicity)

    def _parse(self, force_charge, force_multiplicity):
        """
        Parse the file.
        """
        with open(self.file_path, "r") as fr:
            full_text = fr.read()
        fr.close()
        self._parse_charge_multi(full_text, force_charge, force_multiplicity)
        self._parse_input_parameter(full_text)
        n_atom = int(re.findall(r"NAtoms=\s*(\d+)", full_text)[0])
        blocks = []
        matches = re.search(r"Input orientation:", full_text)
        if matches:
            blocks = re.split(r"Input orientation:", full_text)[1:]
            prefix = "Input orientation:"
        else:
            blocks = re.split(r"Standard orientation:", full_text)[1:]
            prefix = "Standard orientation:"
        if self._only_last_frame:
            blocks = blocks[-1:]
        for block in blocks:
            self.append(
                G16LOGBlockParser(
                    prefix + block,
                    charge=self._charge,
                    multiplicity=self._multiplicity,
                    n_atom=n_atom,
                    file_path=self._file_path,
                    version=self.version,
                    parameter_comment=self._parameter_comment,
                    only_extract_structure=self._only_extract_structure,
                ),
            )

    def _parse_input_parameter(self, full_text):
        pattern = re.compile(
            r"""\s+\*+
\s+(Gaussian\s+\d+\:\s+[A-Za-z0-9-.]+\s+\d+-[A-Za-z]{3}-\d{4})
\s+\d+-[A-Za-z]{3}-\d{4}\s+
\s+\*+
([a-zA-Z%0-9.=\s\_\\\/\*\+\-\"]+)
\s+-+
([a-zA-Z%0-9.\=\s\-\+\#\(\),\*\/\\^\n]+)
\s+-+"""
        )
        matches = pattern.search(full_text)
        if matches:
            version, para_1, para_2 = matches.groups()
        else:
            logger.error(f"No version found in {self._file_path}")
            raise ValueError(f"No version found in {self._file_path}")
        self._version = version
        self._parameter_comment = "\n".join(
            (para_1.replace("\n", " "), para_2.replace("\n", " "))
        )

        (
            self._link0,
            self._route_params,
            self._dieze_tag,
            self._functional,
            self._basis_set,
        ) = parameter_comment_parser(self._parameter_comment)

    def _parse_charge_multi(
        self, full_text, force_charge=None, force_multiplicity=None
    ):
        charge, multi = map(
            int,
            re.search(
                r"Charge\s*=\s*([\-\+\d]+)\s+Multiplicity\s*=\s*(\d+)", full_text
            ).groups(),
        )
        if force_charge is not None:
            charge = force_charge
        if force_multiplicity is not None:
            multi = force_multiplicity
        self._charge, self._multiplicity = charge, multi

    @property
    def link0(self) -> dict:
        return self._link0

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
