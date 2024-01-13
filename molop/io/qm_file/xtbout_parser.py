"""
Author: TMJ
Date: 2024-01-12 09:24:49
LastEditors: TMJ
LastEditTime: 2024-01-12 14:03:48
Description: 请填写简介
"""

import os
import re

from molop.io.bases.file_base import BaseQMFileParser
from molop.io.qm_file.XTBOUTBlockParser import XTBOUTBlockParser
from molop.logger.logger import logger


class XTBOUTParser(BaseQMFileParser):
    _allowed_formats = (".out",)

    def __init__(
        self,
        file_path: str,
        charge=None,
        multiplicity=None,
        show_progress=False,
        only_extract_structure=False,
    ):
        self._check_formats(file_path)
        super().__init__(file_path, show_progress, only_extract_structure)
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
        try:
            version = re.findall(
                r"xtb (version \d+\.\d+\.\d+\s+\([0-9a-z]+\) compiled by ['0-9a-zA-Z\@\_\-]+ on \d{4}-\d{2}-\d{2})",
                full_text,
            )[0]
        except:
            try:
                version = re.findall(
                    r"(Version\s+[0-9\.]+\s+[0-9a-zA-Z]+\s+[\(\)a-zA-Z0-9]+)",
                    full_text,
                )[0]
            except:
                version = ""
        self._version = version
        try:
            charge = int(re.findall(r"charge\s+\:\s+([\-\+0-9]+)", full_text)[0])
        except:
            charge = int(
                float(re.findall(r"total charge\s+([\-\+0-9]+)", full_text)[0])
            )
        if self.__force_charge is not None:
            charge = self.__force_charge
        multi = 1
        if self.__force_multiplicity is not None:
            multi = self.__force_multiplicity
        try:
            self._parameter_comment = " ".join(
                re.findall(
                    "program call\s+\:\s+([0-9a-zA-Z\\/_\-\s\.]+)\n",
                    full_text,
                )[0].split()[2:]
            )
        except:
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
                parameter_comment=self._parameter_comment,
                only_extract_structure=self._only_extract_structure,
            )
        )
