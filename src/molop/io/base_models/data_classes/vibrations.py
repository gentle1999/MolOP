"""Vibrational analysis data classes."""

from __future__ import annotations

from collections.abc import Iterator, Sequence
from typing import Any, ClassVar, Literal, cast, overload

import numpy as np
from pint._typing import UnitLike
from pint.facets.plain import PlainQuantity
from pydantic import Field, computed_field, field_validator, model_validator
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


class Vibration(BaseDataClassWithUnit):
    default_units: ClassVar[dict[str, UnitLike]] = {
        "frequency": atom_ureg.cm_1,
        "reduced_mass": atom_ureg.amu,
        "force_constant": atom_ureg.Unit("mdyne/angstrom"),
        "IR_intensity": atom_ureg.Unit("km/mol"),
        "vibration_mode": atom_ureg.angstrom,
        "raman_activity": atom_ureg.Unit("angstrom**4/amu"),
        "raman_depolarization_ratio": atom_ureg.dimensionless,
        "raman_depolarization_ratio_polarized": atom_ureg.dimensionless,
        "raman_depolarization_ratio_unpolarized": atom_ureg.dimensionless,
    }
    set_default_units: ClassVar[bool] = True

    mode_index: int | None = Field(
        default=None, ge=0, description="Stable 0-based internal mode index"
    )
    source_mode_index: int | None = Field(
        default=None, ge=0, description="Original printed mode number, when available"
    )
    atom_order: Literal["source"] | None = None
    normalization: Literal["unknown", "source_program"] | None = None
    mass_weighting: Literal["unknown", "mass_weighted", "not_mass_weighted"] | None = None
    source_program: str | None = Field(default=None, description="Program producing these modes")
    frequency_method: str | None = None
    raman_activity_method: str | None = None
    orientation: Literal["input", "standard", "source", "unknown"] = "unknown"

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
    raman_activity: PlainQuantity | None = Field(
        default=None,
        description="Raman scattering activity, unit is `angstrom^4/amu`",
        exclude_if=lambda x: x is None,
    )
    raman_depolarization_ratio: PlainQuantity | None = Field(
        default=None,
        description="Dimensionless Raman depolarization ratio; source polarization unspecified",
        exclude_if=lambda x: x is None,
    )
    raman_depolarization_ratio_polarized: PlainQuantity | None = Field(
        default=None,
        description="Raman depolarization ratio for plane-polarized incident light (P)",
        exclude_if=lambda x: x is None,
    )
    raman_depolarization_ratio_unpolarized: PlainQuantity | None = Field(
        default=None,
        description="Raman depolarization ratio for unpolarized incident light (U)",
        exclude_if=lambda x: x is None,
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
        "frequencies": atom_ureg.cm_1,
        "reduced_masses": atom_ureg.amu,
        "force_constants": atom_ureg.Unit("mdyne/angstrom"),
        "IR_intensities": atom_ureg.Unit("km/mol"),
        "raman_activities": atom_ureg.Unit("angstrom**4/amu"),
        "raman_depolarization_ratios": atom_ureg.dimensionless,
        "raman_depolarization_ratios_polarized": atom_ureg.dimensionless,
        "raman_depolarization_ratios_unpolarized": atom_ureg.dimensionless,
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
    raman_activities: PintArrayN = Field(
        default=np.array([]) * atom_ureg.angstrom**4 / atom_ureg.amu,
        description="Raman scattering activity of each mode, unit is `angstrom^4/amu`",
        exclude_if=lambda x: len(x) == 0,
    )
    raman_depolarization_ratios: PintArrayN = Field(
        default=np.array([]) * atom_ureg.dimensionless,
        description="Raman depolarization ratios whose source polarization is unspecified",
        exclude_if=lambda x: len(x) == 0,
    )
    raman_depolarization_ratios_polarized: PintArrayN = Field(
        default=np.array([]) * atom_ureg.dimensionless,
        description="Raman depolarization ratios for plane-polarized incident light (P)",
        exclude_if=lambda x: len(x) == 0,
    )
    raman_depolarization_ratios_unpolarized: PintArrayN = Field(
        default=np.array([]) * atom_ureg.dimensionless,
        description="Raman depolarization ratios for unpolarized incident light (U)",
        exclude_if=lambda x: len(x) == 0,
    )
    mode_indices: list[int] = Field(
        default_factory=list,
        description="Stable 0-based internal mode index for each stored mode",
        exclude_if=lambda x: len(x) == 0,
    )
    source_mode_indices: list[int | None] = Field(
        default_factory=list, description="Original printed mode numbers; None for unnumbered modes"
    )
    source_program: str | None = Field(default=None, description="Program producing these modes")
    frequency_method: str | None = Field(
        default=None, description="Printed method used for frequencies"
    )
    raman_activity_method: str | None = Field(
        default=None, description="Printed method used for Raman activities"
    )
    orientation: Literal["input", "standard", "source", "unknown"] = "unknown"

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

    @field_validator("vibration_modes")
    @classmethod
    def normalize_displacement_units(cls, modes: list[PintArrayNx3]) -> list[PintArrayNx3]:
        return [cast(PintArrayNx3, mode.to(atom_ureg.angstrom)) for mode in modes]

    @model_validator(mode="after")
    def validate_mode_metadata(self) -> Self:
        mode_count = len(self.frequencies)
        for name in (
            "reduced_masses",
            "force_constants",
            "IR_intensities",
            "raman_activities",
            "raman_depolarization_ratios",
            "raman_depolarization_ratios_polarized",
            "raman_depolarization_ratios_unpolarized",
        ):
            if len(getattr(self, name)) not in (0, mode_count):
                raise ValueError(f"{name} must match frequencies")
        if not self.mode_indices and mode_count:
            self.mode_indices = list(range(mode_count))
        if not self.source_mode_indices:
            self.source_mode_indices = [None] * mode_count
        if len(self.source_mode_indices) != mode_count:
            raise ValueError("source_mode_indices must match frequencies")
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
        item_dict: dict[str, Any] = {}
        if isinstance(frame, int):
            frame = range(len(self))[frame]
            item_dict.update(
                mode_index=self.mode_indices[frame],
                source_mode_index=self.source_mode_indices[frame],
                source_program=self.source_program,
                frequency_method=self.frequency_method,
                raman_activity_method=self.raman_activity_method,
                orientation=self.orientation,
                atom_order=self.atom_order,
                normalization=self.normalization,
                mass_weighting=self.mass_weighting,
            )
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
            if len(self.raman_activities) > frame:
                item_dict["raman_activity"] = self.raman_activities[frame]
            if len(self.raman_depolarization_ratios) > frame:
                item_dict["raman_depolarization_ratio"] = self.raman_depolarization_ratios[frame]
            if len(self.raman_depolarization_ratios_polarized) > frame:
                item_dict["raman_depolarization_ratio_polarized"] = (
                    self.raman_depolarization_ratios_polarized[frame]
                )
            if len(self.raman_depolarization_ratios_unpolarized) > frame:
                item_dict["raman_depolarization_ratio_unpolarized"] = (
                    self.raman_depolarization_ratios_unpolarized[frame]
                )
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
        payload = self.model_dump(exclude_computed_fields=True)
        for field in (
            "frequencies",
            "reduced_masses",
            "force_constants",
            "IR_intensities",
            "raman_activities",
            "raman_depolarization_ratios",
            "raman_depolarization_ratios_polarized",
            "raman_depolarization_ratios_unpolarized",
        ):
            values = getattr(self, field)
            if len(values):
                payload[field] = values[imaginary_idxs]
        payload["mode_indices"] = [self.mode_indices[i] for i in imaginary_idxs]
        payload["source_mode_indices"] = [self.source_mode_indices[i] for i in imaginary_idxs]
        payload["vibration_modes"] = (
            [self.vibration_modes[i] for i in imaginary_idxs] if self.vibration_modes else []
        )
        return self.model_validate(payload)

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

    def _mode_metadata(self, index: int) -> dict[str, Any]:
        return {
            "mode_index": self.mode_indices[index],
            "source_mode_index": self.source_mode_indices[index],
            "source_program": self.source_program,
            "frequency_method": self.frequency_method,
            "raman_activity_method": self.raman_activity_method,
            "orientation": self.orientation,
            "normalization": self.normalization,
            "mass_weighting": self.mass_weighting,
            "atom_order": self.atom_order,
        }

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
                    label=f"mode {self.mode_indices[mode_idx]}",
                    center=frequency,
                    intensity=intensity,
                    metadata=self._mode_metadata(mode_idx),
                )
            )
        return Spectrum(
            label=label,
            x_label="wavenumber",
            y_label="IR intensity",
            bands=bands,
            metadata={"source": "Vibrations"},
        )

    def to_raman_spectrum(self, label: str = "Raman") -> Spectrum:
        """Expose Raman activities without assuming an experimental setup."""
        bands: list[SpectralBand] = []
        for index, frequency in enumerate(self.frequencies):
            metadata: dict[str, Any] = self._mode_metadata(index)
            if len(self.raman_depolarization_ratios):
                metadata["depolarization_ratio"] = float(
                    self.raman_depolarization_ratios[index].magnitude
                )
                metadata["depolarization_convention"] = "unspecified"
            for polarization in ("polarized", "unpolarized"):
                values = getattr(self, f"raman_depolarization_ratios_{polarization}")
                if len(values):
                    metadata[f"depolarization_ratio_{polarization}"] = float(
                        values[index].magnitude
                    )
            bands.append(
                SpectralBand(
                    label=f"mode {self.mode_indices[index]}",
                    center=frequency,
                    intensity=self.raman_activities[index] if len(self.raman_activities) else None,
                    metadata=metadata,
                )
            )
        return Spectrum(
            label=label,
            x_label="wavenumber",
            y_label="Raman activity",
            bands=bands,
            metadata={"source": "Vibrations"},
        )

    def to_summary_dict(self, **kwargs) -> SummaryDict:
        return {
            summary_column("Vibration", "num_imaginary"): self.num_imaginary,
            summary_column("Vibration", "num_vibrations"): len(self),
        }


__all__ = ["Vibration", "Vibrations"]
