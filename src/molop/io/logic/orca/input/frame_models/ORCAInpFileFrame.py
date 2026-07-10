from __future__ import annotations

from typing import cast

from pydantic import Field, model_validator
from typing_extensions import Self

from molop.io.base_models.ChemFileFrame import BaseQMInputFrame
from molop.io.base_models.DataClasses import ExplicitSolventRequest
from molop.io.base_models.Mixins import DiskStorageMixin, MemoryStorageMixin
from molop.io.logic.orca.common import (
    ORCABlock,
    ORCACommentLine,
    ORCACommonQMFieldsMixin,
    ORCAExcitedStateSemantic,
    ORCAExplicitSolventSemantic,
    ORCAGeometry,
    ORCAKeywordLine,
    ORCAMultiReferenceSemantic,
    ORCAOutputPrintSetting,
    project_orca_geometry_to_qm_frame,
)


class ORCAInpFileFrameMixin(ORCACommonQMFieldsMixin):
    comment_lines: list[ORCACommentLine] = Field(default_factory=list)
    keyword_lines: list[ORCAKeywordLine] = Field(default_factory=list)
    blocks: list[ORCABlock] = Field(default_factory=list)
    geometry: ORCAGeometry | None = Field(default=None)
    excited_state_semantic: ORCAExcitedStateSemantic = Field(
        default_factory=ORCAExcitedStateSemantic,
        description="Structured excited-state task semantics",
    )
    multi_reference_semantic: ORCAMultiReferenceSemantic = Field(
        default_factory=ORCAMultiReferenceSemantic,
        description="Structured multi-reference task semantics",
    )
    explicit_solvent_semantic: ORCAExplicitSolventSemantic = Field(
        default_factory=ORCAExplicitSolventSemantic,
        description="Structured explicit-solvent semantics",
    )
    explicit_solvent_requests: list[ExplicitSolventRequest] = Field(default_factory=list)
    trailing_lines: list[str] = Field(default_factory=list)
    has_mixed_basis: bool = Field(
        default=False,
        description="Whether atom-level basis overrides were found",
    )
    output_print_settings: list[ORCAOutputPrintSetting] = Field(default_factory=list)

    def _render(self, **kwargs) -> str:
        _ = kwargs
        raise NotImplementedError(
            f"{self.__class__.__name__} does not support ORCA input rendering yet."
        )

    @model_validator(mode="after")
    def set_orca_properties(self) -> Self:
        typed_self = cast(BaseQMInputFrame, self)
        self.has_mixed_basis = project_orca_geometry_to_qm_frame(
            typed_self,
            self.geometry,
            default_version="Any",
        )
        return self


class ORCAInpFileFrameMemory(
    MemoryStorageMixin, ORCAInpFileFrameMixin, BaseQMInputFrame["ORCAInpFileFrameMemory"]
): ...


class ORCAInpFileFrameDisk(
    DiskStorageMixin, ORCAInpFileFrameMixin, BaseQMInputFrame["ORCAInpFileFrameDisk"]
): ...
