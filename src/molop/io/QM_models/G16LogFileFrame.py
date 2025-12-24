"""
Author: TMJ
Date: 2025-07-31 20:27:55
LastEditors: TMJ
LastEditTime: 2025-12-14 21:35:37
Description: 请填写简介
"""

from typing import TYPE_CHECKING, Any, Protocol

import numpy as np
from pint.facets.numpy.quantity import NumpyQuantity
from pydantic import Field, model_validator
from typing_extensions import Self

from molop.io.base_models.ChemFileFrame import BaseCalcFrame
from molop.io.base_models.DataClasses import Vibrations
from molop.io.base_models.Mixins import DiskStorageMixin, MemoryStorageMixin
from molop.io.base_models.Molecule import Molecule
from molop.io.patterns.G16Patterns import (
    SEMI_EMPIRICAL_METHODS,
    route_section_parser,
)
from molop.unit import atom_ureg
from molop.utils.functions import find_rigid_transform, invert_transform_coords, transform_coords


class G16LogFileFrameProtocol(Protocol):
    vibrations: Vibrations | None
    functional: str
    basis_set: str
    method: str

    def log_with_file_info(self, content: str, level: str = "info"): ...


if TYPE_CHECKING:

    class _G16LogFileFrameProtocol(G16LogFileFrameProtocol, Molecule): ...
else:

    class _G16LogFileFrameProtocol(Molecule): ...


class G16LogFileFrameMixin(_G16LogFileFrameProtocol):
    qm_software: str = Field(default="Gaussian")
    options: str = Field(default="", description="options comment")
    title_card: str = Field(default="", description="title card")
    job_type: str = Field(default="", description="Job type")
    keywords: str = Field(
        default="",
        description="Keywords for the QM parameters",
    )

    standard_orientation_transformation_matrix: np.ndarray | None = Field(
        default=None,
        description="Transformation matrix to standard orientation, unit is `angstrom`",
        title="Transformation matrix to standard orientation",
    )
    standard_coords: NumpyQuantity | None = Field(
        default=None,
        description="Atom coordinates with standard orientation, unit is `angstrom`",
        title="Atom coordinates with standard orientation",
    )

    def _add_default_units(self) -> None:
        super()._add_default_units()
        self._default_units.update({"standard_orientation_coords": atom_ureg.angstrom})

    @model_validator(mode="after")
    def _post_processing(self) -> Self:
        if self.standard_coords is not None and (len(self.coords) == len(self.standard_coords)):
            self.standard_orientation_transformation_matrix = find_rigid_transform(
                self.coords.m, self.standard_coords.m
            )

        if len(self.coords) != len(self.atoms):  # no input orientation found
            if self.standard_coords is not None:
                if self.standard_orientation_transformation_matrix is not None:
                    self.coords = (
                        invert_transform_coords(
                            self.standard_coords.m,
                            self.standard_orientation_transformation_matrix,
                        )
                        * self.standard_coords.u
                    )
                    # self.log_with_file_info(
                    #     "To get the correcct input orientation, add `Geom=PrintInputOrient` in the keywords.",
                    #     level="warning",
                    # )
                elif self.standard_orientation_transformation_matrix is None:
                    self.coords = self.standard_coords
            else:  # no standard coords found
                raise ValueError(
                    "The number of atoms and coordinates do not match, "
                    "and the standard orientation is not provided."
                )
        if (
            self.standard_coords is None
            and self.standard_orientation_transformation_matrix is not None
        ):
            self.standard_coords = (
                transform_coords(self.coords.m, self.standard_orientation_transformation_matrix)
                * self.coords.u
            )
        if (
            self.standard_orientation_transformation_matrix is not None
            and self.vibrations is not None
        ):
            self.vibrations.transform_orientation(
                self.standard_orientation_transformation_matrix, inverse=True
            )

        if self.basis_set.lower() == "genecp":
            self.basis_set = "pseudopotential"
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
        if "em" in self.route_params:
            self.functional = f"{self.functional}-{self.route_params['em'].upper()}"
        if "empiricaldispersion" in self.route_params:
            self.functional = (
                f"{self.functional}-{self.route_params['empiricaldispersion'].upper()}"
            )
        return self

    @property
    def route_params(self) -> dict[str, Any]:
        return route_section_parser(self.keywords)[0]

    @property
    def dieze_tag(self) -> str | None:
        return route_section_parser(self.keywords)[1]


class G16LogFileFrameMemory(
    MemoryStorageMixin, G16LogFileFrameMixin, BaseCalcFrame["G16LogFileFrameMemory"]
): ...


class G16LogFileFrameDisk(
    DiskStorageMixin,
    G16LogFileFrameMixin,
    BaseCalcFrame["G16LogFileFrameDisk"],
): ...
