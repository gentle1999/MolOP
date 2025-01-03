"""
Author: TMJ
Date: 2024-06-21 11:03:06
LastEditors: TMJ
LastEditTime: 2024-06-24 22:10:46
Description: 请填写简介
"""

import re
from typing import Literal

from pydantic import Field, computed_field

from molop.io.bases.BaseMolFileParser import BaseQMMolFileParser
from molop.io.qm_file.G16LogFrameParser import G16LogFrameParser
from molop.logger.logger import moloplogger
from molop.unit import atom_ureg
from molop.utils.g16patterns import (
    g16logpatterns,
    get_solvent,
    get_solvent_model,
    link0_parser,
    parameter_comment_parser,
    semi_empirical_methods,
)


class G16LogFileParser(BaseQMMolFileParser[G16LogFrameParser]):
    _allowed_formats = (".log", ".g16", ".gal", ".out", ".irc")
    qm_software: str = Field(default="Gaussian")

    def _parse(self):
        """
        Parse the file.
        """
        with open(self.file_path, "r") as fr:
            full_text = fr.read()
        fr.close()
        self._parse_charge_multi(full_text)
        self._parse_input_parameter(full_text)
        self._parse_tail(full_text)
        self._add_em()
        temperature_match = g16logpatterns["Temperature"].search(full_text)
        if temperature_match:
            self.temperature = float(temperature_match.group(1))
        n_atom = self._get_n_atoms(full_text)
        frames = []
        matches = re.search(r"Input orientation:", full_text)
        if matches:
            frames = re.split(r"Input orientation:", full_text)[1:]
            prefix = "Input orientation:"
        else:
            frames = re.split(r"Standard orientation:", full_text)[1:]
            prefix = "Standard orientation:"
        if self.only_last_frame:
            frames = frames[-1:]
        for frame in frames:
            self.__append__(
                G16LogFrameParser(
                    frame_content=prefix + frame,
                    charge=self.charge,
                    multiplicity=self.multiplicity,
                    n_atom=n_atom,
                    file_path=self.file_path,
                    qm_software_version=self.qm_software_version,
                    keywords=self.keywords,
                    method=self.method,
                    functional=self.functional,
                    basis=self.basis,
                    only_extract_structure=self.only_extract_structure,
                ),
            )
        self._get_running_time()

    def _get_n_atoms(self, full_text: str):
        n_atom_match = re.search(g16logpatterns["n atoms"], full_text)
        if n_atom_match:
            n_atom = int(n_atom_match.group(1))
        else:
            if re.search(r"Symbolic Z-matrix:\n", full_text):
                n_atom = -1
                for line in full_text.splitlines():
                    if line != " \n":
                        n_atom += 1
                    else:
                        break
            else:
                raise RuntimeError("No atoms found in the file.")
        return n_atom

    def _parse_input_parameter(self, full_text: str):
        version_match = g16logpatterns["version"].search(full_text)
        if version_match:
            version = version_match.group()
        else:
            moloplogger.error(f"No version found in {self.file_path}")
            raise ValueError(f"No version found in {self.file_path}")
        self.qm_software_version = version
        link = ""
        route = ""
        link_start = False
        route_start = False
        for row in full_text[version_match.end() :].splitlines(True):
            if "----" in row and route_start:
                route_start = False
                break
            if route_start:
                route += row
            if "----" in row and not route_start:
                link_start = False
                route_start = True
            if link_start:
                link += row
            if "****" in row:
                link_start = True
        self.keywords = "\n\n".join((link.replace("\n", " "), route.replace("\n ", "")))
        link = link.replace("\n", " ")
        route = route.replace("\n ", "")
        self.__link0 = link0_parser(link)

        (
            self.__route_params,
            self.__dieze_tag,
        ) = parameter_comment_parser(route)
        self.solvent_model = get_solvent_model(self.route_params)
        self.solvent = get_solvent(self.route_params)
        self._parse_functional_basis(full_text)

    def _parse_charge_multi(self, full_text: str):
        charge, multi = map(
            int,
            re.search(
                r"Charge\s*=\s*([\-\+\d]+)\s+Multiplicity\s*=\s*(\d+)", full_text
            ).groups(),
        )
        self.charge, self.multiplicity = charge, multi

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
        if self.functional.endswith("hf"):
            self.method = "HF"

    def _parse_tail(self, full_text: str):
        functional_ref = g16logpatterns["functional"].search(full_text)
        start = g16logpatterns["tail_start"].search(full_text)
        end = g16logpatterns["tail_end"].search(full_text)
        if start and end:
            tail = full_text[start.end() : end.start()]
            tail = tail.replace("\n ", "").replace("|", "\\")

            groups = tail.split("\\")
            self.functional = groups[2].lower()
            self.basis = groups[3].lower()
            if self.basis == "genecp":
                self.basis = "pseudopotential"
            if self.functional[1:] in semi_empirical_methods:
                self.method = "SEMI-EMPIRICAL"
            elif functional_ref and functional_ref.group(1).lower().endswith("hf"):
                if self.functional.endswith("hf"):
                    self.method = "HF"
                else:
                    self.method = "POST-HF"
            else:
                self.method = "DFT"

    def _add_em(self):
        if "em" in self.route_params.keys():
            self.functional = f"{self.functional}-{self.route_params['em']}"
        if "empiricaldispersion" in self.route_params.keys():
            self.functional = (
                f"{self.functional}-{self.route_params['empiricaldispersion']}"
            )

    def _get_running_time(self):
        self.running_time = sum(frame.running_time for frame in self.frames)

    @property
    def link0(self) -> dict:
        return self.__link0

    @property
    def route_params(self) -> dict:
        return self.__route_params

    @property
    def dieze_tag(self) -> Literal["#N", "#P", "#T"]:
        return self.__dieze_tag

    @computed_field
    @property
    def task_type(self) -> Literal["sp", "opt", "freq"]:
        if self.route_params:
            if "opt" in self.route_params:
                return "opt"
            if "freq" in self.route_params:
                return "freq"
        return "sp"
