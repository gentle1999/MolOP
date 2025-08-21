"""
Author: TMJ
Date: 2025-07-31 20:27:55
LastEditors: TMJ
LastEditTime: 2025-08-02 22:13:41
Description: 请填写简介
"""

from typing import Any, Optional

import numpy as np
from pint.facets.numpy.quantity import NumpyQuantity
from pydantic import Field, model_validator
from typing_extensions import Self

from molop.io.base_models.ChemFileFrame import BaseCalcFrame
from molop.io.base_models.Mixins import DiskStorageMixin, MemoryStorageMixin
from molop.io.base_models.Molecule import Molecule
from molop.io.patterns.G16Patterns import (
    SEMI_EMPIRICAL_METHODS,
    parameter_comment_parser,
)
from molop.unit import atom_ureg
from molop.utils.functions import invert_transform_coords


class G16LogFileFrameMixin(Molecule):
    qm_software: str = Field(default="Gaussian")
    options: str = Field(default="", description="options comment")
    title_card: str = Field(default="", description="title card")
    job_type: str = Field(default="", description="Job type")
    keywords: str = Field(
        default="",
        description="Keywords for the QM parameters",
    )

    standard_orientation_transformation_matrix: Optional[np.ndarray] = Field(
        default=None,
        description="Transformation matrix to standard orientation, unit is `angstrom`",
        title="Transformation matrix to standard orientation",
    )
    standard_orientation_coords: Optional[NumpyQuantity] = Field(
        default=None,
        description="Atom coordinates with standard orientation, unit is `angstrom`",
        title="Atom coordinates with standard orientation",
    )

    def _add_default_units(self) -> None:
        super()._add_default_units()
        self._default_units.update({"standard_orientation_coords": atom_ureg.angstrom})

    @model_validator(mode="after")
    def _post_processing(self) -> Self:
        if len(self.coords) != len(self.atoms):
            if (
                self.standard_orientation_coords is not None
                and self.standard_orientation_transformation_matrix is not None
            ):
                self.coords = (
                    invert_transform_coords(
                        self.standard_orientation_coords.m,
                        self.standard_orientation_transformation_matrix,
                    )
                    * self.standard_orientation_coords.u
                )
            else:
                raise ValueError(
                    "The number of atoms and coordinates do not match, and the standard orientation is not provided."
                )
        if self.basis.lower() == "genecp":
            self.basis = "pseudopotential"
        for semi in SEMI_EMPIRICAL_METHODS:
            if semi in self.functional.lower():
                self.method = "SEMI-EMPIRICAL"
                break
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


class G16LogFileFrameMemory(
    MemoryStorageMixin, BaseCalcFrame["G16LogFileFrameMemory"], G16LogFileFrameMixin
): ...


class G16LogFileFrameDisk(
    DiskStorageMixin, BaseCalcFrame["G16LogFileFrameDisk"], G16LogFileFrameMixin
):
    _allowed_formats_ = (".log", ".g16", ".gal", ".out", ".irc")
