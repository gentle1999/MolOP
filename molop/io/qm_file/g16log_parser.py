'''
Author: TMJ
Date: 2024-01-09 20:19:06
LastEditors: TMJ
LastEditTime: 2024-01-16 10:39:02
Description: 请填写简介
'''
import os
import re

from molop.io.bases.file_base import BaseQMFileParser
from molop.io.qm_file.G16LOGBlockParser import G16LOGBlockParser
from molop.logger.logger import logger


class G16LOGParser(BaseQMFileParser):
    _allowed_formats = (".log",)

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
        """
        Parse the file.
        """
        with open(self.file_path, "r") as fr:
            lines = fr.readlines()
        fr.close()
        full_text = "".join(lines)
        charge, multi = map(
            int,
            re.findall(
                r"Charge\s*=\s*([\-\+\d]+)\s+Multiplicity\s*=\s*(\d+)", full_text
            )[0],
        )
        if self.__force_charge is not None:
            charge = self.__force_charge
        if self.__force_multiplicity is not None:
            multi = self.__force_multiplicity
        pattern = r"""\s+\*+
\s+(Gaussian\s+\d+\:\s+[A-Za-z0-9-.]+\s+\d+-[A-Za-z]{3}-\d{4})
\s+\d+-[A-Za-z]{3}-\d{4}\s+
\s+\*+
([a-zA-Z%0-9.=\s\_\\\/\*\+\-]+)
\s+-+
([a-zA-Z%0-9.\=\s\-\+\#\(\),\*\/\\^\n]+)
\s+-+"""
        try:
            version, para_1, para_2 = re.findall(pattern, full_text)[0]
        except:
            logger.error(f"No version found in {self._file_path}")
            raise ValueError(f"No version found in {self._file_path}")
        self._version = version
        self._parameter_comment = "\n".join((para_1, para_2))
        n_atom = int(re.findall(r"NAtoms=\s*(\d+)", full_text)[0])
        block_starts = [
            idx for idx, line in enumerate(lines) if "Input orientation:" in line
        ] + [len(lines)]
        if self._only_last_frame:
            block_starts = block_starts[-2:]
        for idx, start in enumerate(block_starts[:-1]):
            self.append(
                G16LOGBlockParser(
                    "".join(lines[start : block_starts[idx + 1]]),
                    charge=charge,
                    multiplicity=multi,
                    n_atom=n_atom,
                    version=version,
                    parameter_comment=self._parameter_comment,
                    only_extract_structure=self._only_extract_structure,
                ),
            )
