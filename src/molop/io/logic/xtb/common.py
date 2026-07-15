from __future__ import annotations

from typing import Any, ClassVar

from pint._typing import UnitLike
from pint.facets.plain import PlainQuantity
from pydantic import Field, model_validator
from typing_extensions import Self

from molop.unit import atom_ureg


def populate_common_xtb_qm_fields(target: Any) -> None:
    target.qm_software = "xTB"
    target.backfill_common_qm_containers_from_legacy()
    target.project_common_qm_fields()


class XTBOutputQMFieldsMixin:
    """Fields shared by xTB output file and frame models."""

    default_units: ClassVar[dict[str, UnitLike]] = {
        "gradient_norm": atom_ureg.Unit("hartree / bohr"),
        "gradient_norm_threshold": atom_ureg.Unit("hartree / bohr"),
    }

    qm_software: str = Field(default="xTB")
    input_file_name: str = Field(
        default="",
        description="Coordinate file name printed in the xTB calculation setup",
    )
    gradient_norm: PlainQuantity | None = Field(
        default=None,
        description="Final source-reported xTB gradient norm",
        exclude_if=lambda value: value is None,
    )
    gradient_norm_threshold: PlainQuantity | None = Field(
        default=None,
        description="Source-reported xTB gradient convergence threshold",
        exclude_if=lambda value: value is None,
    )

    @model_validator(mode="after")
    def _normalize_common_xtb_fields(self) -> Self:
        populate_common_xtb_qm_fields(self)
        return self


__all__ = ["XTBOutputQMFieldsMixin", "populate_common_xtb_qm_fields"]
