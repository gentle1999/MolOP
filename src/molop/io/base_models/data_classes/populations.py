"""Atomic population and spin data classes."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any, Literal, cast

import numpy as np
from pydantic import ConfigDict, Field, model_validator
from typing_extensions import Self

from molop.io.base_models.Bases import (
    BaseDataClassWithUnit,
    PropertyBundle,
    PropertyColumnValue,
    PropertyScalarValue,
    PropertyTable,
)
from molop.io.base_models.summary import SummaryDict, summary_dict_from_fields


class AtomicPopulationSeries(BaseDataClassWithUnit):
    """One extensible atom-aligned population series."""

    model_config = ConfigDict(arbitrary_types_allowed=True, extra="forbid")

    scheme: str = Field(min_length=1, description="Population analysis scheme")
    quantity: str = Field(
        min_length=1,
        description="Physical quantity, for example charge or spin_density",
    )
    values: list[float] = Field(min_length=1, description="Values in source atom order")
    spin_channel: Literal["alpha", "beta", "total"] | None = Field(
        default=None,
        description="Spin channel when the source resolves one",
    )
    source_label: str | None = Field(
        default=None,
        description="Exact or normalized source table/record label",
    )
    metadata: dict[str, Any] = Field(
        default_factory=dict,
        description="Format-specific population metadata",
    )


class ChargeSpinPopulations(BaseDataClassWithUnit):
    model_config = ConfigDict(arbitrary_types_allowed=True, extra="forbid")

    populations: dict[str, AtomicPopulationSeries] = Field(
        default_factory=dict,
        description=(
            "Extensible atom-aligned population series keyed by a stable snake_case identifier"
        ),
        exclude_if=lambda x: len(x) == 0,
    )

    def population_items(self) -> list[tuple[str, AtomicPopulationSeries]]:
        """Return all population series in insertion order."""

        return list(self.populations.items())

    def __getitem__(self, name: str) -> AtomicPopulationSeries:
        return self.populations[name]

    def __len__(self) -> int:
        return len(self.populations)

    def get_population(self, name: str) -> AtomicPopulationSeries | None:
        """Return one population series by stable key."""

        return self.populations.get(name)

    @property
    def population_names(self) -> list[str]:
        return list(self.populations)

    def to_population_table(self, atom_symbols: Sequence[str] | None = None) -> PropertyTable:
        population_items = self.population_items()
        population_lengths = [len(series.values) for _name, series in population_items]
        if population_lengths:
            num_atoms = population_lengths[0]
        elif atom_symbols is not None:
            num_atoms = len(atom_symbols)
        else:
            num_atoms = 0

        if atom_symbols is not None and len(atom_symbols) != num_atoms:
            raise ValueError(
                f"atom_symbols length {len(atom_symbols)} does not match population length {num_atoms}"
            )

        columns: dict[str, PropertyColumnValue] = {}
        if num_atoms:
            columns["atom_index"] = np.arange(num_atoms)
        if atom_symbols is not None:
            columns["atom_symbol"] = list(atom_symbols)
        for population_name, series in population_items:
            columns[population_name] = cast(PropertyColumnValue, series.values)

        return PropertyTable(
            columns=columns,
            row_labels=[str(i) for i in range(num_atoms)],
            metadata={
                "source": "ChargeSpinPopulations",
                "index_base": 0,
                "population_names": self.population_names,
            },
        )

    def to_property_bundle(self, atom_symbols: Sequence[str] | None = None) -> PropertyBundle:
        table = self.to_population_table(atom_symbols=atom_symbols)
        return PropertyBundle(
            tables={"atomic_populations": table} if len(table) else {},
            metadata={"source": "ChargeSpinPopulations"},
        )

    @model_validator(mode="after")
    def validate_charge_spin_populations(self) -> Self:
        population_items = self.population_items()
        if not population_items:
            return self
        population_lengths = {name: len(series.values) for name, series in population_items}
        expected_length = next(iter(population_lengths.values()))
        mismatched = {
            name: length for name, length in population_lengths.items() if length != expected_length
        }
        if mismatched:
            raise ValueError(
                "All populations must have the same length; "
                f"expected {expected_length}, got {mismatched}"
            )
        return self

    def to_summary_dict(self, **kwargs) -> SummaryDict:
        return {}


class TotalSpin(BaseDataClassWithUnit):
    spin_square: float | None = Field(
        default=None,
        description="Spin square of the molecule",
        exclude_if=lambda x: x is None,
    )
    spin_quantum_number: float | None = Field(
        default=None,
        description="Spin quantum number of the molecule",
        exclude_if=lambda x: x is None,
    )

    def to_summary_dict(self, **kwargs) -> SummaryDict:
        return summary_dict_from_fields(self, "TotalSpin", **kwargs)

    def to_property_bundle(self) -> PropertyBundle:
        scalar_properties: dict[str, PropertyScalarValue] = {
            key: value
            for key, value in {
                "spin_square": self.spin_square,
                "spin_quantum_number": self.spin_quantum_number,
            }.items()
            if value is not None
        }
        return PropertyBundle(
            scalar_properties=scalar_properties,
            metadata={"source": "TotalSpin"},
        )


__all__ = ["AtomicPopulationSeries", "ChargeSpinPopulations", "TotalSpin"]
