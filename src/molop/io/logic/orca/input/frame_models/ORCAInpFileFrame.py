from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any, cast

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
from molop.io.logic.orca.input._orca_inp_renderer import (
    adapt_orca_writer_payload,
    render_orca_input_frame,
)


class ORCAInpFileFrameMixin(ORCACommonQMFieldsMixin):
    @classmethod
    def adapt_writer_payload(cls, data: dict[str, Any]) -> dict[str, Any]:
        return adapt_orca_writer_payload(data)

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

    def _render(
        self,
        keywords: str | Sequence[str] | None = None,
        nprocs: int | None = None,
        maxcore: int | None = None,
        blocks: str | Mapping[str, Any] | Sequence[ORCABlock] | None = None,
        charge: int | None = None,
        multiplicity: int | None = None,
        coordinate_decimal_places: int = 10,
        **kwargs: Any,
    ) -> str:
        """Render a canonical ORCA input frame.

        Explicit render arguments replace matching structured values for this
        render only. Unknown ORCA blocks can be supplied as raw block text or a
        mapping whose keys are block names.
        """

        _ = kwargs
        typed_self = cast(BaseQMInputFrame, self)
        return render_orca_input_frame(
            fallback_keywords=typed_self.keywords,
            fallback_keyword_lines=self.keyword_lines,
            fallback_blocks=self.blocks,
            fallback_comments=self.comment_lines,
            fallback_trailing_lines=self.trailing_lines,
            geometry=self.geometry,
            keywords=keywords,
            nprocs=nprocs,
            maxcore=maxcore,
            blocks=blocks,
            charge=charge,
            multiplicity=multiplicity,
            coordinate_decimal_places=coordinate_decimal_places,
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
