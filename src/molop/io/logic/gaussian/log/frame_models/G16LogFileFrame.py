"""
Author: TMJ
Date: 2025-07-31 20:27:55
LastEditors: TMJ
LastEditTime: 2026-02-04 15:47:33
Description: 请填写简介
"""

from typing import ClassVar, cast

from pint._typing import UnitLike
from pydantic import Field, PrivateAttr, model_validator
from typing_extensions import Self

from molop.io.base_models.ChemFileFrame import BaseCalcFrame, _HasVibrations
from molop.io.base_models.Mixins import DiskStorageMixin, MemoryStorageMixin
from molop.io.logic.gaussian.input.GaussianRoute import (
    GaussianRouteSemanticFieldsMixin,
)
from molop.io.logic.gaussian.log.frame_models.G16Components import (
    G16ComponentTree,
    G16ComponentTreeBuilder,
    _rawify_payload,
)
from molop.unit import atom_ureg
from molop.utils.functions import find_rigid_transform, invert_transform_coords, transform_coords
from molop.utils.types import Array4x4, PintArrayNx3


class G16LogFileFrameMixin(GaussianRouteSemanticFieldsMixin):
    default_units: ClassVar[dict[str, UnitLike]] = {
        "standard_orientation_coords": atom_ureg.angstrom
    }
    standard_orientation_transformation_matrix: Array4x4 | None = Field(
        default=None,
        description="Transformation matrix to standard orientation, unit is `angstrom`",
        title="Transformation matrix to standard orientation",
    )
    standard_coords: PintArrayNx3 | None = Field(
        default=None,
        description="Atom coordinates with standard orientation, unit is `angstrom`",
        title="Atom coordinates with standard orientation",
    )
    _component_tree: G16ComponentTree | None = PrivateAttr(default=None)

    @staticmethod
    def _prepare_component_tree(tree: G16ComponentTree) -> G16ComponentTree:
        G16ComponentTreeBuilder.expand_synthetic_children(tree)
        if not tree.payloads_rawified:
            for node in tree.iter_nodes():
                if node.component is not None:
                    node.component.payload = _rawify_payload(node.component.payload)
            tree.payloads_rawified = True
        return tree

    def _get_private_component_tree(self) -> G16ComponentTree | None:
        private_attrs = getattr(self, "__pydantic_private__", None)
        if isinstance(private_attrs, dict):
            tree = private_attrs.get("_component_tree")
            return tree
        tree = getattr(self, "_component_tree", None)
        return tree

    def _set_private_component_tree(self, tree: G16ComponentTree | None) -> None:
        private_attrs = getattr(self, "__pydantic_private__", None)
        if isinstance(private_attrs, dict):
            private_attrs["_component_tree"] = tree
        else:
            object.__setattr__(self, "_component_tree", tree)

    @model_validator(mode="after")
    def _post_processing(self) -> Self:
        typed_self = cast(_HasVibrations, self)
        if self.standard_coords is not None and (
            len(typed_self.coords) == len(self.standard_coords)
        ):
            self.standard_orientation_transformation_matrix = find_rigid_transform(
                typed_self.coords.m, self.standard_coords.m
            )

        if len(typed_self.coords) != len(typed_self.atoms):  # no input orientation found
            if self.standard_coords is not None:
                if self.standard_orientation_transformation_matrix is not None:
                    typed_self.coords = (
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
                    typed_self.coords = self.standard_coords
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
                transform_coords(
                    typed_self.coords.m, self.standard_orientation_transformation_matrix
                )
                * typed_self.coords.u
            )
        if (
            self.standard_orientation_transformation_matrix is not None
            and typed_self.vibrations is not None
        ):
            typed_self.vibrations.transform_orientation(
                self.standard_orientation_transformation_matrix, inverse=True
            )

        return self

    @property
    def component_tree(self) -> G16ComponentTree:
        tree = self._get_private_component_tree()
        if tree is None:
            tree = G16ComponentTreeBuilder.from_frame_data(self)
        tree = self._prepare_component_tree(tree)
        self._set_private_component_tree(tree)
        return tree

    def render_fakeg(self, **kwargs) -> str:
        return self.component_tree.render_fakeg(**kwargs)


class G16LogFileFrameMemory(
    MemoryStorageMixin, G16LogFileFrameMixin, BaseCalcFrame["G16LogFileFrameMemory"]
): ...


class G16LogFileFrameDisk(
    DiskStorageMixin, G16LogFileFrameMixin, BaseCalcFrame["G16LogFileFrameDisk"]
): ...


G16LogFileFrameMemory.model_rebuild()
G16LogFileFrameDisk.model_rebuild()
