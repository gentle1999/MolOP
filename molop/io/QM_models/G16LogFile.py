"""
Author: TMJ
Date: 2025-07-31 20:27:42
LastEditors: TMJ
LastEditTime: 2025-08-20 10:26:43
Description: 请填写简介
"""

from typing import Any, Optional

import numpy as np
from pydantic import Field, model_validator
from typing_extensions import Self

from molop.io.base_models.Bases import BaseDataClassWithUnit
from molop.io.base_models.ChemFile import CalcFileMixin
from molop.io.base_models.Mixins import DiskStorageMixin, MemoryStorageMixin
from molop.io.patterns.G16Patterns import (
    SEMI_EMPIRICAL_METHODS,
    parameter_comment_parser,
)
from molop.io.QM_models.G16LogFileFrame import (
    G16LogFileFrameDisk,
    G16LogFileFrameMemory,
)


class G16LogFileMixin(BaseDataClassWithUnit):
    qm_software: str = Field(default="Gaussian")
    options: str = Field(default="", description="options comment")
    title_card: str = Field(default="", description="title card")
    job_type: str = Field(default="", description="Job type")
    keywords: str = Field(
        default="",
        description="Keywords for the QM parameters",
    )
    method: str = Field(
        default="",
        description="QM method used to perform the calculation. e.g. DFT or GFN2-xTB",
    )
    basis: str = Field(
        default="",
        description="Basis set used in the QM calculation, only for DFT calculations",
    )
    functional: str = Field(
        default="",
        description="Functional used in the QM calculation, only for DFT calculations",
    )

    standard_orientation_transformation_matrix: Optional[np.ndarray] = Field(
        default=None,
        description="Transformation matrix to standard orientation, unit is `angstrom`",
        title="Transformation matrix to standard orientation",
    )

    @model_validator(mode="after")
    def _validate_metadata(self) -> Self:
        if self.basis.lower() == "genecp":
            self.basis = "pseudopotential"
        if any(semi in self.keywords.lower() for semi in SEMI_EMPIRICAL_METHODS):
            self.method = "SEMI-EMPIRICAL"
        else:
            if self.functional.lower().endswith("hf"):
                self.method = "HF"
            if self.functional.lower().endswith("fc"):
                self.method = "FC"
            else:
                self.method = "DFT"
        if "em" in self.route_params.keys():
            self.functional = f"{self.functional}-{self.route_params['em'].upper()}"
        if "empiricaldispersion" in self.route_params.keys():
            self.functional = (
                f"{self.functional}-{self.route_params['empiricaldispersion'].upper()}"
            )
        return self

    @property
    def route_params(self) -> dict[str, Any]:
        return parameter_comment_parser(self.keywords)[0]

    @property
    def dieze_tag(self) -> str | None:
        return parameter_comment_parser(self.keywords)[1]


class G16LogFileMemory(
    MemoryStorageMixin, CalcFileMixin[G16LogFileFrameMemory], G16LogFileMixin
): ...


class G16LogFileDisk(
    DiskStorageMixin, CalcFileMixin[G16LogFileFrameDisk], G16LogFileMixin
):
    _allowed_formats_ = (".log", ".g16", ".gal", ".out", ".irc")
