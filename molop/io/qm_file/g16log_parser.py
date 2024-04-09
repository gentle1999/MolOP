"""
Author: TMJ
Date: 2024-02-17 15:17:37
LastEditors: TMJ
LastEditTime: 2024-02-17 20:03:06
Description: 请填写简介
"""
import re
from typing import Literal

from molop.io.bases.file_base import BaseQMFileParser
from molop.io.qm_file.G16LOGBlockParser import G16LOGBlockParser
from molop.logger.logger import logger
from molop.utils.g16patterns import (
    g16logpatterns,
    get_solvent,
    get_solvent_model,
    link0_parser,
    parameter_comment_parser,
)


class G16LOGParser(BaseQMFileParser):
    _allowed_formats = (".log", ".g16", ".gal", ".out", ".irc")

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
        self._post_parse()

    def _parse(self, force_charge, force_multiplicity):
        """
        Parse the file.
        """
        with open(self.file_path, "r") as fr:
            full_text = fr.read()
        fr.close()
        self._parse_charge_multi(full_text, force_charge, force_multiplicity)
        self._parse_input_parameter(full_text)
        self._parse_functional_basis(full_text)
        temperature_match = g16logpatterns["Temperature"].search(full_text)
        if temperature_match:
            self.temperature = float(temperature_match.group(1))
        n_atom_match = re.search(g16logpatterns["n atoms"], full_text)
        if n_atom_match:
            n_atom = int(n_atom_match.group(1))
        else:
            raise RuntimeError("No atoms found in the file.")
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
                    parameter_comment=self.parameter_comment,
                    functional=self.functional,
                    basis=self.basis,
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
        self.version = version
        self.parameter_comment = "\n".join(
            (para_1.replace("\n", " "), para_2.replace("\n ", ""))
        )
        link, route = self.parameter_comment.split("#")
        link = link.replace("\n", " ")
        route = "#" + route.replace("\n", " ")
        self._link0 = link0_parser(link)

        (
            self._route_params,
            self._dieze_tag,
        ) = parameter_comment_parser(route)
        self.solvent_model = get_solvent_model(self.route_params)
        self.solvent = get_solvent(self.route_params)

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

    def _parse_functional_basis(self, full_text: str):
        functional_match = g16logpatterns["functional"].search(full_text)
        if functional_match:
            self.functional = functional_match.group(1).lower()
        else:
            self.functional = "unknown"
        basis_match = g16logpatterns["basis"].search(full_text)
        if basis_match:
            self.basis = basis_match.group(1)
        else:
            self.basis = "unknown"
        if g16logpatterns["Pseudopotential"].search(full_text):
            self.basis = "pseudopotential"

    @property
    def link0(self) -> dict:
        return self._link0

    @property
    def route_params(self) -> dict:
        return self._route_params

    @property
    def dieze_tag(self) -> Literal["#N", "#P", "#T"]:
        return self._dieze_tag
