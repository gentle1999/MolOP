"""Calculation status and geometry-convergence data classes."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any, ClassVar

import numpy as np
import pandas as pd
from pint._typing import UnitLike
from pint.facets.plain import PlainQuantity
from pydantic import Field, computed_field, model_validator
from typing_extensions import Self

from molop.io.base_models.Bases import BaseDataClassWithUnit
from molop.io.base_models.summary import (
    SummaryDict,
    summary_column,
    summary_dict_from_fields,
    summary_item,
)
from molop.unit import atom_ureg


class GeometryOptimizationStatus(BaseDataClassWithUnit):
    """
    Geometry optimization status.
    """

    optimization_metrics: ClassVar[tuple[str, ...]] = (
        "energy_change",
        "rms_force",
        "max_force",
        "rms_displacement",
        "max_displacement",
    )
    default_units: ClassVar[dict[str, UnitLike]] = {
        "energy_change": atom_ureg.hartree,
        "energy_change_threshold": atom_ureg.hartree,
        "rms_force": atom_ureg.Unit("hartree/bohr"),
        "rms_force_threshold": atom_ureg.Unit("hartree/bohr"),
        "max_force": atom_ureg.Unit("hartree/bohr"),
        "max_force_threshold": atom_ureg.Unit("hartree/bohr"),
        "rms_displacement": atom_ureg.bohr,
        "rms_displacement_threshold": atom_ureg.bohr,
        "max_displacement": atom_ureg.bohr,
        "max_displacement_threshold": atom_ureg.bohr,
    }
    set_default_units: ClassVar[bool] = True

    source_converged: dict[str, bool | None] | None = Field(
        default=None,
        description="Optional source-reported convergence decisions keyed by metric",
        exclude_if=lambda value: value is None,
    )
    source_labels: dict[str, str] | None = Field(
        default=None,
        description="Optional source labels keyed by normalized optimization metric",
        exclude_if=lambda value: value is None,
    )

    geometry_optimized: bool | None = Field(
        default=None,
        description="Whether the geometry has been optimized",
        exclude_if=lambda x: x is None,
    )
    convergence_multiplier: float = Field(
        default=2.0,
        ge=1.0,
        description="Tolerance multiplier used to judge acceptable geometry optimization convergence.",
    )
    energy_change_threshold: PlainQuantity | None = Field(
        default=None,
        description="Energy change threshold",
        exclude_if=lambda x: x is None,
    )
    rms_force_threshold: PlainQuantity | None = Field(
        default=None,
        description="RMS force threshold in internal some programs use gradient, which has the same absolute value",
        exclude_if=lambda x: x is None,
    )
    max_force_threshold: PlainQuantity | None = Field(
        default=None,
        description="Maximum force threshold in internal some programs use gradient, which has the same absolute value",
        exclude_if=lambda x: x is None,
    )
    rms_displacement_threshold: PlainQuantity | None = Field(
        default=None,
        description="RMS displacement threshold in internal",
        exclude_if=lambda x: x is None,
    )
    max_displacement_threshold: PlainQuantity | None = Field(
        default=None,
        description="Maximum displacement threshold in internal",
        exclude_if=lambda x: x is None,
    )
    energy_change: PlainQuantity | None = Field(
        default=None,
        description="Energy change",
        exclude_if=lambda x: x is None,
    )
    rms_force: PlainQuantity | None = Field(
        default=None,
        description="RMS force some programs use gradient, which has the same absolute value",
        exclude_if=lambda x: x is None,
    )
    max_force: PlainQuantity | None = Field(
        default=None,
        description="Maximum force some programs use gradient, which has the same absolute value",
        exclude_if=lambda x: x is None,
    )
    rms_displacement: PlainQuantity | None = Field(
        default=None,
        description="RMS displacement",
        exclude_if=lambda x: x is None,
    )
    max_displacement: PlainQuantity | None = Field(
        default=None,
        description="Maximum displacement",
        exclude_if=lambda x: x is None,
    )

    @model_validator(mode="before")
    @classmethod
    def attach_default_metric_units(cls, data: Any) -> Any:
        if not isinstance(data, Mapping):
            return data
        normalized = dict(data)
        for field, unit in cls.default_units.items():
            value = normalized.get(field)
            if value is not None and not hasattr(value, "units"):
                normalized[field] = value * unit
        return normalized

    def _metric_converged(self, metric: str) -> bool | None:
        value = getattr(self, metric)
        threshold = getattr(self, f"{metric}_threshold")
        if value is None or threshold is None:
            return None
        return bool(abs(value) <= threshold * self.convergence_multiplier)

    @computed_field()  # type: ignore[prop-decorator]
    @property
    def energy_change_converged(self) -> bool | None:
        """
        Whether the energy change has converged.
        """
        return self._metric_converged("energy_change")

    @computed_field()  # type: ignore[prop-decorator]
    @property
    def rms_force_converged(self) -> bool | None:
        """
        Whether the RMS force has converged.
        """
        return self._metric_converged("rms_force")

    @computed_field()  # type: ignore[prop-decorator]
    @property
    def max_force_converged(self) -> bool | None:
        """
        Whether the maximum force has converged.
        """
        return self._metric_converged("max_force")

    @computed_field()  # type: ignore[prop-decorator]
    @property
    def rms_displacement_converged(self) -> bool | None:
        """
        Whether the RMS displacement has converged.
        """
        return self._metric_converged("rms_displacement")

    @computed_field()  # type: ignore[prop-decorator]
    @property
    def max_displacement_converged(self) -> bool | None:
        """
        Whether the maximum displacement has converged.
        """
        return self._metric_converged("max_displacement")

    def not_converged_num(self) -> int:
        """
        Return the number of not converged properties.
        """
        return sum(
            1
            for metric in (
                self.energy_change_converged,
                self.rms_force_converged,
                self.max_force_converged,
                self.rms_displacement_converged,
                self.max_displacement_converged,
            )
            if metric is False
        )

    def __vector__(self) -> np.ndarray[Any, np.dtype[np.floating[Any]]]:
        return np.array(
            [
                abs(getattr(self, metric).magnitude)
                if getattr(self, metric) is not None
                else np.nan
                for metric in self.optimization_metrics
            ],
            dtype=float,
        )

    def to_df(self) -> pd.DataFrame:
        metrics_list = list(self.optimization_metrics)
        df = pd.DataFrame(
            index=pd.Index(metrics_list),
            columns=pd.Index(["value", "threshold", "converged"]),
        )
        df["value"] = [getattr(self, metric) for metric in self.optimization_metrics]
        df["threshold"] = [
            getattr(self, f"{metric}_threshold") for metric in self.optimization_metrics
        ]
        df["converged"] = [
            getattr(self, f"{metric}_converged") for metric in self.optimization_metrics
        ]
        return df

    def __le__(self, other: GeometryOptimizationStatus):
        # TODO: more accurate comparison
        if not isinstance(other, GeometryOptimizationStatus):
            raise NotImplementedError
        if self.not_converged_num() > other.not_converged_num():
            return False
        self_vector = self.__vector__()
        other_vector = other.__vector__()
        shared = np.isfinite(self_vector) & np.isfinite(other_vector)
        if not shared.any():
            return True
        return bool(np.all(self_vector[shared] <= other_vector[shared]))

    def __gt__(self, other: GeometryOptimizationStatus):
        return not self.__le__(other)

    @model_validator(mode="after")
    def __check_geometry_optimized__(self) -> Self:
        status = (
            self.energy_change_converged,
            self.rms_force_converged,
            self.max_force_converged,
            self.rms_displacement_converged,
            self.max_displacement_converged,
        )
        if self.geometry_optimized is True:
            return self
        known_statuses = [metric_status for metric_status in status if metric_status is not None]
        if known_statuses:
            self.geometry_optimized = all(known_statuses)
        return self

    def to_summary_dict(self, **kwargs) -> SummaryDict:
        summary = {
            summary_column(
                "GeometryOptimizationStatus", "geometry_optimized"
            ): self.geometry_optimized,
            summary_column(
                "GeometryOptimizationStatus", "convergence_multiplier"
            ): self.convergence_multiplier,
        }
        for metric in self.optimization_metrics:
            if (value := getattr(self, metric)) is not None:
                item = summary_item("GeometryOptimizationStatus", metric, value)
                if item is not None:
                    column, summary_value = item
                    summary[column] = summary_value
            if (threshold := getattr(self, f"{metric}_threshold")) is not None:
                item = summary_item(
                    "GeometryOptimizationStatus",
                    f"{metric}_threshold",
                    threshold,
                )
                if item is not None:
                    column, summary_value = item
                    summary[column] = summary_value
        return summary


class Status(BaseDataClassWithUnit):
    scf_converged: bool | None = Field(
        default=None,
        description="Whether the SCF has converged, or None when the output has no SCF evidence",
    )
    normal_terminated: bool | None = Field(
        default=None,
        description=(
            "Whether the calculation segment terminated normally, or None when no "
            "segment-termination evidence is attached to this object"
        ),
    )

    def to_summary_dict(self, **kwargs) -> SummaryDict:
        return summary_dict_from_fields(self, "Status", **kwargs)


__all__ = ["GeometryOptimizationStatus", "Status"]
