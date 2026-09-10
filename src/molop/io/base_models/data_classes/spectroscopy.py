"""NMR, shielding, and solvation data classes."""

from __future__ import annotations

from typing import ClassVar, Literal

import numpy as np
from pint._typing import UnitLike
from pint.facets.plain import PlainQuantity
from pydantic import Field, model_validator
from typing_extensions import Self

from molop.io.base_models.Bases import BaseDataClassWithUnit
from molop.io.base_models.summary import SummaryDict, summary_dict_from_fields
from molop.unit import atom_ureg
from molop.utils.types import PintArray3, PintArray3x3, PintSquareMatrix


class ShieldingTensor(BaseDataClassWithUnit):
    """Nuclear magnetic shielding tensor for one source-order atom."""

    default_units: ClassVar[dict[str, UnitLike]] = {
        "shielding_tensor": atom_ureg.ppm,
        "isotropic": atom_ureg.ppm,
        "anisotropy": atom_ureg.ppm,
        "principal_values": atom_ureg.ppm,
    }
    set_default_units: ClassVar[bool] = True

    atom_index: int = Field(ge=0, description="Zero-based atom index in source atom order")
    atom_symbol: str = Field(min_length=1, description="Atom element symbol")
    shielding_tensor: PintArray3x3 = Field(
        description="Full 3 x 3 magnetic shielding tensor, unit is `ppm`",
    )
    isotropic: PlainQuantity | None = Field(
        default=None,
        description="Observed isotropic shielding, or trace(tensor) / 3 when absent",
    )
    anisotropy: PlainQuantity | None = Field(
        default=None,
        description="Observed shielding anisotropy, or the Haeberlen value when absent",
    )
    principal_values: PintArray3 | None = Field(
        default=None,
        description="Observed or derived three principal shielding values",
    )
    anisotropy_convention: str | None = Field(
        default=None,
        description="Convention used for the anisotropy value",
    )
    orientation: Literal["input", "standard", "source", "unknown"] = Field(
        default="unknown",
        description="Cartesian orientation used by the tensor components",
    )

    @model_validator(mode="after")
    def validate_shielding_tensor(self) -> Self:
        tensor_magnitude = np.asarray(self.shielding_tensor.magnitude, dtype=float)
        if not np.isfinite(tensor_magnitude).all():
            raise ValueError("shielding_tensor must contain only finite values")

        if self.isotropic is None:
            self.isotropic = float(np.trace(tensor_magnitude) / 3.0) * self.shielding_tensor.units

        if self.principal_values is None:
            symmetric_tensor = (tensor_magnitude + tensor_magnitude.T) / 2.0
            eigenvalues = np.linalg.eigvalsh(symmetric_tensor)
            self.principal_values = (
                np.sort(np.asarray(eigenvalues, dtype=float)) * self.shielding_tensor.units
            )

        principal_values_quantity = self.principal_values
        assert principal_values_quantity is not None
        if self.anisotropy is None:
            principal_values = np.asarray(principal_values_quantity.magnitude, dtype=float)
            isotropic = float(np.mean(principal_values))
            ordered = sorted(
                principal_values,
                key=lambda value: abs(float(value) - isotropic),
                reverse=True,
            )
            self.anisotropy = (
                float(ordered[0] - (ordered[1] + ordered[2]) / 2.0)
                * principal_values_quantity.units
            )
            self.anisotropy_convention = self.anisotropy_convention or "Haeberlen"
        return self

    def to_summary_dict(self, **kwargs) -> SummaryDict:
        return summary_dict_from_fields(self, "ShieldingTensor", **kwargs)


class NMR(BaseDataClassWithUnit):
    default_units: ClassVar[dict[str, UnitLike]] = {
        "spin_spin_coupling_k": atom_ureg.Hz,
        "spin_spin_coupling_j": atom_ureg.Hz,
    }
    set_default_units: ClassVar[bool] = True

    gauge: str | None = Field(
        default=None,
        description="Magnetic gauge/origin scheme, such as GIAO or CSGT",
    )
    shielding_tensors: list[ShieldingTensor] = Field(
        default_factory=list, description="NMR shielding tensors"
    )
    coupling_atom_indices: list[int] = Field(
        default_factory=list,
        description="Source-order atom indices corresponding to coupling matrix axes",
    )
    spin_spin_coupling_k: PintSquareMatrix | None = Field(
        default=None,
        description="Reduced isotropic spin-spin coupling matrix, unit is `Hz`",
    )
    spin_spin_coupling_j: PintSquareMatrix | None = Field(
        default=None,
        description="Isotropic spin-spin coupling matrix, unit is `Hz`",
    )
    spin_spin_coupling_k_components: dict[Literal["FC", "SD", "PSO", "DSO"], PintSquareMatrix] = (
        Field(
            default_factory=dict,
            description="Reduced spin-spin coupling contribution matrices in Hz",
        )
    )
    spin_spin_coupling_j_components: dict[Literal["FC", "SD", "PSO", "DSO"], PintSquareMatrix] = (
        Field(
            default_factory=dict,
            description="Spin-spin coupling contribution matrices in Hz",
        )
    )

    @model_validator(mode="after")
    def validate_nmr(self) -> Self:
        shielding_indices = [item.atom_index for item in self.shielding_tensors]
        if len(shielding_indices) != len(set(shielding_indices)):
            raise ValueError("shielding_tensors contain duplicate atom indices")

        matrix_sizes: set[int] = set()
        for field_name in ("spin_spin_coupling_k", "spin_spin_coupling_j"):
            value = getattr(self, field_name)
            if value is None:
                continue
            matrix_sizes.add(int(value.shape[0]))
        for field_name in (
            "spin_spin_coupling_k_components",
            "spin_spin_coupling_j_components",
        ):
            components = getattr(self, field_name)
            for component_name, value in components.items():
                normalized = value.to(atom_ureg.Hz)
                components[component_name] = normalized
                matrix_sizes.add(int(normalized.shape[0]))
        if len(matrix_sizes) > 1:
            raise ValueError(
                "spin-spin coupling total and component matrices must have matching shapes"
            )
        if self.coupling_atom_indices:
            if len(self.coupling_atom_indices) != len(set(self.coupling_atom_indices)):
                raise ValueError("coupling_atom_indices must be unique")
            if matrix_sizes and len(self.coupling_atom_indices) != next(iter(matrix_sizes)):
                raise ValueError("coupling_atom_indices length must match coupling matrix size")
        return self

    def to_summary_dict(self, **kwargs) -> SummaryDict:
        return summary_dict_from_fields(self, "NMR", **kwargs)


__all__ = ["NMR", "ShieldingTensor"]
