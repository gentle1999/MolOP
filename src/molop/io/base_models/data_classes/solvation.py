"""Implicit-solvation result data classes."""

from __future__ import annotations

from pydantic import Field

from molop.io.base_models.Bases import BaseDataClassWithUnit
from molop.io.base_models.summary import SummaryDict, summary_dict_from_fields


class ImplicitSolvation(BaseDataClassWithUnit):
    solvent: str | None = Field(
        default=None,
        description="Solvent used in the QM calculation",
        exclude_if=lambda x: x is None,
    )
    solvent_model: str | None = Field(
        default=None,
        description="Solvent model used in the QM calculation",
        exclude_if=lambda x: x is None,
    )
    atomic_radii: str | None = Field(
        default=None,
        description="Atomic radii used in the QM calculation",
        exclude_if=lambda x: x is None,
    )
    solvent_epsilon: float | None = Field(
        default=None,
        description="Solvent dielectric constant used in the QM calculation",
        exclude_if=lambda x: x is None,
    )
    solvent_epsilon_infinite: float | None = Field(
        default=None,
        description="Solvent epsilon infinite used in the QM calculation",
        exclude_if=lambda x: x is None,
    )

    def to_summary_dict(self, **kwargs) -> SummaryDict:
        return summary_dict_from_fields(self, "ImplicitSolvation", **kwargs)


__all__ = ["ImplicitSolvation"]
