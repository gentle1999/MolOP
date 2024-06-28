import re
from typing import Literal

from pydantic import Field

from molop.io.bases.BaseMolFileParser import BaseQMMolFileParser
from molop.io.qm_file.G16FchkFrameParser import G16FchkFrameParser
from molop.utils.g16patterns import (
    g16fchkpatterns,
    get_solvent,
    get_solvent_model,
    parameter_comment_parser,
)


class G16FchkFileParser(BaseQMMolFileParser[G16FchkFrameParser]):
    _allowed_formats = (".fchk", ".fck", ".fch")
    qm_software: str = Field(default="Gaussian")

    def _parse(self):
        with open(self.file_path, "r") as fr:
            full_text = fr.read()
        fr.close()

        self.charge = int(re.search(g16fchkpatterns["charge"], full_text).group(1))
        self.multiplicity = int(re.search(g16fchkpatterns["multi"], full_text).group(1))
        n_atoms = int(re.search(g16fchkpatterns["n_atoms"], full_text).group(1))
        self.qm_software_version = re.search(
            g16fchkpatterns["version"], full_text
        ).group(1)
        self.keywords = (
            re.search(g16fchkpatterns["route"], full_text).group(1).replace("\n", "")
        )
        self._parse_functional_basis(full_text)
        (
            self.__route_params,
            self.__dieze_tag,
        ) = parameter_comment_parser(self.keywords)
        self.solvent_model = get_solvent_model(self.route_params)
        self.solvent = get_solvent(self.route_params)
        self.__append__(
            G16FchkFrameParser(
                frame_content=full_text,
                charge=self.charge,
                multiplicity=self.multiplicity,
                n_atom=n_atoms,
                file_path=self.file_path,
                version=self.qm_software_version,
                keywords=self.keywords,
                only_extract_structure=self.only_extract_structure,
            )
        )

    def _parse_functional_basis(self, full_text: str):
        for idx, line in enumerate(full_text.splitlines()):
            if idx == 1:
                self.functional = line.split()[1].lower()
                self.basis = line.split()[2]
                break

    @property
    def route_params(self) -> dict:
        return self.__route_params

    @property
    def dieze_tag(self) -> Literal["#N", "#P", "#T"]:
        return self.__dieze_tag
