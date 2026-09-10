"""Vibrational analysis data classes."""

from __future__ import annotations

from collections.abc import Iterator, Sequence
from typing import Any, ClassVar, Literal, cast, overload

import numpy as np
from pint._typing import UnitLike
from pint.facets.plain import PlainQuantity
from pydantic import Field, computed_field, model_validator
from typing_extensions import Self

from molop.io.base_models.Bases import BaseDataClassWithUnit, SpectralBand, Spectrum
from molop.io.base_models.summary import (
    SummaryDict,
    summary_column,
    summary_dict_from_fields,
)
from molop.unit import atom_ureg
from molop.utils.functions import invert_transform_displacements, transform_displacements
from molop.utils.types import PintArrayN, PintArrayNx3


def _is_quantity(value: Any) -> bool:
    return isinstance(value, PlainQuantity)


class Vibration(BaseDataClassWithUnit):
    default_units: ClassVar[dict[str, UnitLike]] = {
        "frequency": atom_ureg.cm_1,
        "reduced_mass": atom_ureg.amu,
        "force_constant": atom_ureg.Unit("mdyne/angstrom"),
        "IR_intensity": atom_ureg.Unit("km/mol"),
        "vibration_mode": atom_ureg.angstrom,
    }
    set_default_units: ClassVar[bool] = True

    frequency: PlainQuantity | None = Field(
        default=None,
        description="Frequency of each mode, unit is `cm^-1`",
        exclude_if=lambda x: x is None,
    )
    reduced_mass: PlainQuantity | None = Field(
        default=None,
        description="Reduced mass of each mode, unit is `amu`",
        exclude_if=lambda x: x is None,
    )
    force_constant: PlainQuantity | None = Field(
        default=None,
        description="Force constant of each mode, unit is `mdyne/angstrom`",
        exclude_if=lambda x: x is None,
    )
    IR_intensity: PlainQuantity | None = Field(
        default=None,
        description="IR intensity of each mode, unit is `km/mol`",
        exclude_if=lambda x: x is None,
    )
    vibration_mode: PintArrayNx3 = Field(
        default=np.zeros((0, 3)) * atom_ureg.angstrom,
        description="Vibration mode of each mode, unit is `angstrom`",
        exclude_if=lambda x: len(x) == 0,
    )

    @computed_field()  # type: ignore[prop-decorator]
    @property
    def is_imaginary(self) -> bool:
        return bool(self.frequency is not None and cast(Any, self.frequency) < 0)

    def to_summary_dict(self, **kwargs) -> SummaryDict:
        return summary_dict_from_fields(self, "Vibration", **kwargs)

    def transform_orientation(
        self, transformation_matrix: np.ndarray, inverse: bool = False
    ) -> None:
        if inverse:
            mode = cast(Any, self.vibration_mode)
            self.vibration_mode = (
                invert_transform_displacements(mode.m, transformation_matrix) * mode.u
            )
        else:
            mode = cast(Any, self.vibration_mode)
            self.vibration_mode = transform_displacements(mode.m, transformation_matrix) * mode.u


class Vibrations(BaseDataClassWithUnit, Sequence[Vibration]):
    default_units: ClassVar[dict[str, UnitLike]] = {
        "frequency": atom_ureg.cm_1,
        "reduced_mass": atom_ureg.amu,
        "force_constant": atom_ureg.Unit("mdyne/angstrom"),
        "IR_intensity": atom_ureg.Unit("km/mol"),
        "vibration_mode": atom_ureg.angstrom,
    }
    set_default_units: ClassVar[bool] = True

    frequencies: PintArrayN = Field(
        default=np.array([]) * atom_ureg.cm_1,
        description="Frequency of each mode, unit is `cm^-1`",
        exclude_if=lambda x: len(x) == 0,
    )
    reduced_masses: PintArrayN = Field(
        default=np.array([]) * atom_ureg.amu,
        description="Reduced mass of each mode, unit is `amu`",
        exclude_if=lambda x: len(x) == 0,
    )
    force_constants: PintArrayN = Field(
        default=np.array([]) * atom_ureg.mdyne / atom_ureg.angstrom,
        description="Force constant of each mode, unit is `mdyne/angstrom`",
        exclude_if=lambda x: len(x) == 0,
    )
    IR_intensities: PintArrayN = Field(
        default=np.array([]) * atom_ureg.km / atom_ureg.mol,
        description="IR intensity of each mode, unit is `km/mol`",
        exclude_if=lambda x: len(x) == 0,
    )
    vibration_modes: list[PintArrayNx3] = Field(
        default_factory=list,
        description="Vibration mode of each mode, unit is `angstrom`",
        exclude_if=lambda x: len(x) == 0,
    )
    mode_indices: list[int] = Field(
        default_factory=list,
        description="Source frequency-mode index for each stored mode",
        exclude_if=lambda x: len(x) == 0,
    )
    axis_order: tuple[Literal["mode"], Literal["atom"], Literal["cartesian"]] | None = Field(
        default=None,
        description="Axis order of the conceptual stacked vibration_modes array",
        exclude_if=lambda x: x is None,
    )
    atom_order: Literal["source"] | None = Field(
        default=None,
        description="Atom ordering used by vibration_modes",
        exclude_if=lambda x: x is None,
    )
    normalization: Literal["unknown", "source_program"] | None = Field(
        default=None,
        description="Normalization convention of vibration_modes",
        exclude_if=lambda x: x is None,
    )
    mass_weighting: Literal["unknown", "mass_weighted", "not_mass_weighted"] | None = Field(
        default=None,
        description="Mass-weighting convention of vibration_modes",
        exclude_if=lambda x: x is None,
    )

    @model_validator(mode="after")
    def validate_mode_metadata(self) -> Self:
        mode_count = len(self.frequencies)
        if not self.mode_indices and mode_count:
            self.mode_indices = list(range(mode_count))
        if self.mode_indices and len(self.mode_indices) != mode_count:
            raise ValueError("mode_indices must match frequencies")
        if len(set(self.mode_indices)) != len(self.mode_indices) or any(
            index < 0 for index in self.mode_indices
        ):
            raise ValueError("mode_indices must be unique non-negative indices")
        if self.vibration_modes:
            if len(self.vibration_modes) != mode_count:
                raise ValueError("vibration_modes must match frequencies")
            mode_shapes = {tuple(mode.shape) for mode in self.vibration_modes}
            if len(mode_shapes) != 1 or next(iter(mode_shapes))[-1:] != (3,):
                raise ValueError("vibration_modes must share an (atom, 3) shape")
            self.axis_order = self.axis_order or ("mode", "atom", "cartesian")
            self.atom_order = self.atom_order or "source"
            self.normalization = self.normalization or "unknown"
            self.mass_weighting = self.mass_weighting or "unknown"
        return self

    def __iter__(self) -> Iterator[Vibration]:  # type: ignore[override]
        for i in range(len(self)):
            yield self[i]

    @overload
    def __getitem__(self, frame: int) -> Vibration: ...

    @overload
    def __getitem__(self, frame: slice) -> list[Vibration]: ...

    @overload
    def __getitem__(self, frame: Sequence) -> list[Vibration]: ...

    def __getitem__(self, frame: int | slice | Sequence) -> Vibration | list[Vibration]:
        item_dict = {}
        if isinstance(frame, int):
            if len(self.frequencies) > frame:
                item_dict["frequency"] = self.frequencies[frame]
            if len(self.reduced_masses) > frame:
                item_dict["reduced_mass"] = self.reduced_masses[frame]
            if len(self.force_constants) > frame:
                item_dict["force_constant"] = self.force_constants[frame]
            if len(self.IR_intensities) > frame:
                item_dict["IR_intensity"] = self.IR_intensities[frame]
            if len(self.vibration_modes) > frame:
                item_dict["vibration_mode"] = self.vibration_modes[frame]
            return Vibration.model_validate(item_dict)
        if isinstance(frame, slice):
            return [self[idx] for idx in range(*frame.indices(len(self.frequencies)))]
        return [self[idx] for idx in frame]

    def __len__(self) -> int:
        return len(self.frequencies)

    @property
    def num_imaginary(self) -> int:
        return len(self.imaginary_idxs)

    @property
    def imaginary_idxs(self) -> list[int]:
        return [i for i, freq in enumerate(self) if freq.is_imaginary]

    @property
    def imaginary_vibrations(self) -> Vibrations:
        imaginary_idxs = self.imaginary_idxs
        return self.model_validate(
            {
                "frequencies": (
                    cast(Any, self.frequencies)[imaginary_idxs]
                    if _is_quantity(self.frequencies)
                    else None
                ),
                "reduced_masses": (
                    cast(Any, self.reduced_masses)[imaginary_idxs]
                    if _is_quantity(self.reduced_masses)
                    else None
                ),
                "force_constants": (
                    cast(Any, self.force_constants)[imaginary_idxs]
                    if _is_quantity(self.force_constants)
                    else None
                ),
                "IR_intensities": (
                    cast(Any, self.IR_intensities)[imaginary_idxs]
                    if _is_quantity(self.IR_intensities)
                    else None
                ),
                "vibration_modes": (
                    [self.vibration_modes[i] for i in imaginary_idxs]
                    if len(self.vibration_modes)
                    else []
                ),
                "mode_indices": [self.mode_indices[i] for i in imaginary_idxs],
                "axis_order": self.axis_order,
                "atom_order": self.atom_order,
                "normalization": self.normalization,
                "mass_weighting": self.mass_weighting,
            }
        )

    def transform_orientation(
        self, transformation_matrix: np.ndarray, inverse: bool = False
    ) -> None:
        if inverse:
            self.vibration_modes = [
                invert_transform_displacements(cast(Any, mode).m, transformation_matrix)
                * cast(Any, mode).u
                for mode in self.vibration_modes
            ]
        else:
            self.vibration_modes = [
                transform_displacements(cast(Any, mode).m, transformation_matrix)
                * cast(Any, mode).u
                for mode in self.vibration_modes
            ]

    def to_ir_spectrum(self, label: str = "IR") -> Spectrum:
        bands: list[SpectralBand] = []
        for mode_idx, frequency in enumerate(self.frequencies):
            intensity = (
                cast(Any, self.IR_intensities)[mode_idx]
                if len(self.IR_intensities) > mode_idx
                else None
            )
            bands.append(
                SpectralBand(
                    label=f"mode {mode_idx + 1}",
                    center=frequency,
                    intensity=intensity,
                    metadata={"mode_index": self.mode_indices[mode_idx]},
                )
            )
        return Spectrum(
            label=label,
            x_label="wavenumber",
            y_label="IR intensity",
            bands=bands,
            metadata={"source": "Vibrations"},
        )

    def to_summary_dict(self, **kwargs) -> SummaryDict:
        return {
            summary_column("Vibration", "num_imaginary"): self.num_imaginary,
            summary_column("Vibration", "num_vibrations"): len(self),
        }


__all__ = ["Vibration", "Vibrations"]
