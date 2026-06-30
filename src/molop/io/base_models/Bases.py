"""
Author: TMJ
Date: 2025-07-26 19:08:33
LastEditors: TMJ
LastEditTime: 2026-06-30 16:05:44
Description: 请填写简介
"""

from collections.abc import Iterator, Mapping, Sequence
from typing import Any, ClassVar, TypeAlias, cast, overload

import numpy as np
import pandas as pd
from pint._typing import UnitLike
from pint.facets.numpy.quantity import NumpyQuantity
from pint.facets.plain import PlainQuantity
from pydantic import BaseModel, ConfigDict, Field, model_validator
from typing_extensions import Self

from molop.config import molopconfig, moloplogger
from molop.unit import unit_transform


PropertyScalarValue: TypeAlias = (
    str | int | float | bool | PlainQuantity | NumpyQuantity | None
)
PropertyColumnValue: TypeAlias = (
    list[str | int | float | bool | None] | np.ndarray | PlainQuantity | NumpyQuantity
)


class BaseDataClassWithUnit(BaseModel):
    model_config = ConfigDict(arbitrary_types_allowed=True)
    default_units: ClassVar[dict[str, UnitLike]] = {}
    set_default_units: ClassVar[bool] = False

    @staticmethod
    def __unitless_dump__(item: Any, **kwargs) -> Any:
        if hasattr(item, "m") and hasattr(item, "units"):
            if isinstance(item.m, np.ndarray):
                return item.m.tolist()
            else:
                return item.m
        if isinstance(item, np.ndarray):
            return item.tolist()
        if isinstance(item, list):
            return [BaseDataClassWithUnit.__unitless_dump__(i, **kwargs) for i in item]
        if isinstance(item, tuple):
            return tuple(BaseDataClassWithUnit.__unitless_dump__(i, **kwargs) for i in item)
        if isinstance(item, Mapping):
            return {
                k: BaseDataClassWithUnit.__unitless_dump__(v, **kwargs) for k, v in item.items()
            }
        if isinstance(item, BaseDataClassWithUnit):
            return item.to_unitless_dump(**kwargs)
        return item

    def to_unitless_dump(self, **kwargs) -> dict[str, Any]:
        """
        Parameters:
            kwargs:
                follow the same format as
                [`model_dump`](https://docs.pydantic.dev/latest/concepts/serialization/) method.
        """
        return {
            k: BaseDataClassWithUnit.__unitless_dump__(getattr(self, k), **kwargs)
            for k, v in self.model_dump(**kwargs).items()
        }

    @model_validator(mode="after")
    def __unit_transform__(self) -> Self:
        if self.set_default_units or molopconfig.force_unit_transform:
            self._transform_units(self.default_units)
            if moloplogger.isEnabledFor(10):
                moloplogger.debug(f"Data class {self.__class__.__name__} parsed.")
        return self

    def _transform_units(self, unit_dict: Mapping[str, UnitLike]) -> None:
        for key, unit in unit_dict.items():
            if hasattr(self, key):
                setattr(self, key, unit_transform(getattr(self, key), unit))

    def to_summary_series(self, **kwargs) -> pd.Series:
        return pd.Series(self.to_summary_dict(**kwargs))

    def to_summary_dict(self, **kwargs) -> dict[tuple[str, str], Any]:
        return {}

    @classmethod
    def __pydantic_init_subclass__(cls, **kwargs):
        super().__pydantic_init_subclass__(**kwargs)
        final_default_units: dict[str, UnitLike] = {}
        for base in reversed(cls.__mro__):
            if hasattr(base, "default_units"):
                base_default_units = getattr(base, "default_units", {})
                final_default_units.update(base_default_units)
        cls.default_units = final_default_units


class PropertyPoint(BaseDataClassWithUnit):
    label: str | None = Field(default=None, description="Point label or identifier")
    value: PlainQuantity | NumpyQuantity | float | None = Field(
        default=None, description="Point value, stored as a pint quantity or dimensionless float"
    )
    uncertainty: PlainQuantity | NumpyQuantity | float | None = Field(
        default=None,
        description="Optional uncertainty, stored as a pint quantity or dimensionless float",
    )
    metadata: dict[str, Any] = Field(default_factory=dict, description="Extra point metadata")


class PropertySeries(BaseDataClassWithUnit, Sequence[PropertyPoint]):
    points: list[PropertyPoint] = Field(default_factory=list, description="Ordered property points")
    axis_label: str | None = Field(default=None, description="Independent-axis label")
    value_label: str | None = Field(default=None, description="Dependent-axis label")
    metadata: dict[str, Any] = Field(default_factory=dict, description="Series-level metadata")

    def __iter__(self) -> Iterator[PropertyPoint]:  # type: ignore[override]
        return iter(self.points)

    def __len__(self) -> int:
        return len(self.points)

    @overload
    def __getitem__(self, index: int) -> PropertyPoint: ...

    @overload
    def __getitem__(self, index: slice) -> list[PropertyPoint]: ...

    def __getitem__(self, index: int | slice) -> PropertyPoint | list[PropertyPoint]:
        return self.points[index]


class TensorProperty(BaseDataClassWithUnit):
    label: str | None = Field(default=None, description="Property label")
    tensor: NumpyQuantity | np.ndarray | None = Field(default=None, description="Tensor payload")
    frame: str | None = Field(default=None, description="Reference frame")
    metadata: dict[str, Any] = Field(default_factory=dict, description="Extra tensor metadata")


class PropertyTable(BaseDataClassWithUnit):
    columns: dict[str, PropertyColumnValue] = Field(
        default_factory=dict,
        description="Column-oriented property table. Quantity columns keep pint units.",
    )
    row_labels: list[str] = Field(default_factory=list, description="Optional row labels")
    metadata: dict[str, Any] = Field(default_factory=dict, description="Table-level metadata")

    @staticmethod
    def _column_length(column: PropertyColumnValue) -> int:
        if isinstance(column, list | np.ndarray):
            return len(column)
        magnitude = column.magnitude
        if isinstance(magnitude, np.ndarray):
            return len(magnitude)
        return 1

    @staticmethod
    def _column_item(column: PropertyColumnValue, row_index: int) -> PropertyScalarValue:
        if isinstance(column, list):
            return column[row_index]
        if isinstance(column, np.ndarray):
            item = column[row_index]
            if isinstance(item, np.generic):
                return cast(PropertyScalarValue, item.item())
            return cast(PropertyScalarValue, item)
        return cast(PropertyScalarValue, cast(NumpyQuantity, column)[row_index])

    @model_validator(mode="after")
    def validate_property_table(self) -> Self:
        if not self.columns:
            return self
        column_lengths = {name: self._column_length(column) for name, column in self.columns.items()}
        expected_length = next(iter(column_lengths.values()))
        mismatched = {
            name: length for name, length in column_lengths.items() if length != expected_length
        }
        if mismatched:
            raise ValueError(f"PropertyTable columns have inconsistent lengths: {mismatched}")
        if self.row_labels and len(self.row_labels) != expected_length:
            raise ValueError(
                f"PropertyTable row_labels length {len(self.row_labels)} does not match "
                f"column length {expected_length}"
            )
        return self

    def __len__(self) -> int:
        if not self.columns:
            return 0
        return self._column_length(next(iter(self.columns.values())))

    def __getitem__(self, column_name: str) -> PropertyColumnValue:
        return self.columns[column_name]

    @property
    def column_names(self) -> list[str]:
        return list(self.columns)

    def row(self, row_index: int) -> dict[str, PropertyScalarValue]:
        if row_index < 0:
            row_index += len(self)
        if row_index < 0 or row_index >= len(self):
            raise IndexError(f"PropertyTable row index {row_index} out of range")

        row_data: dict[str, PropertyScalarValue] = {}
        for column_name, column in self.columns.items():
            row_data[column_name] = self._column_item(column, row_index)
        return row_data


class PropertyTransition(BaseDataClassWithUnit):
    label: str | None = Field(default=None, description="Transition label")
    initial_state: int | None = Field(default=None, description="Initial state index")
    final_state: int | None = Field(default=None, description="Final state index")
    energy: PlainQuantity | None = Field(default=None, description="Transition energy")
    wavelength: PlainQuantity | None = Field(default=None, description="Transition wavelength")
    oscillator_strength: float | None = Field(default=None, description="Oscillator strength")
    rotatory_strength: float | None = Field(default=None, description="Rotatory strength")
    transition_dipole: NumpyQuantity | None = Field(default=None, description="Transition dipole")
    properties: dict[str, PropertyScalarValue] = Field(
        default_factory=dict, description="Extra transition properties"
    )


class SpectralBand(BaseDataClassWithUnit):
    label: str | None = Field(default=None, description="Band label")
    center: PlainQuantity | None = Field(default=None, description="Band center")
    width: PlainQuantity | None = Field(default=None, description="Band width")
    intensity: PlainQuantity | float | None = Field(default=None, description="Band intensity")
    assignment: str | None = Field(default=None, description="Band assignment")
    metadata: dict[str, Any] = Field(default_factory=dict, description="Extra band metadata")


class Spectrum(BaseDataClassWithUnit):
    label: str | None = Field(default=None, description="Spectrum label")
    x_label: str | None = Field(default=None, description="X-axis label")
    y_label: str | None = Field(default=None, description="Y-axis label")
    bands: list[SpectralBand] = Field(default_factory=list, description="Spectrum bands")
    transitions: list[PropertyTransition] = Field(
        default_factory=list, description="Resolved transitions"
    )
    metadata: dict[str, Any] = Field(default_factory=dict, description="Spectrum metadata")


class PropertyBundle(BaseDataClassWithUnit):
    scalar_properties: dict[str, PropertyScalarValue] = Field(
        default_factory=dict, description="Scalar results"
    )
    vector_properties: dict[str, NumpyQuantity] = Field(
        default_factory=dict, description="Vector results"
    )
    tensor_properties: dict[str, TensorProperty] = Field(
        default_factory=dict, description="Tensor-valued results"
    )
    tables: dict[str, PropertyTable] = Field(default_factory=dict, description="Tabular results")
    spectra: dict[str, Spectrum] = Field(default_factory=dict, description="Spectral results")
    series: dict[str, PropertySeries] = Field(default_factory=dict, description="Generic series")
    transitions: dict[str, list[PropertyTransition]] = Field(
        default_factory=dict, description="Transition collections"
    )
    metadata: dict[str, Any] = Field(default_factory=dict, description="Bundle metadata")
