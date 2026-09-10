"""Single-point molecular property data classes."""

from __future__ import annotations

from collections.abc import Sequence
from typing import ClassVar, cast

import numpy as np
from pint._typing import UnitLike
from pint.facets.numpy.quantity import NumpyQuantity
from pint.facets.plain import PlainQuantity
from pydantic import Field

from molop.io.base_models.Bases import (
    BaseDataClassWithUnit,
    PropertyBundle,
    PropertyColumnValue,
    PropertyScalarValue,
    PropertyTable,
    TensorProperty,
)
from molop.io.base_models.summary import SummaryDict, summary_dict_from_fields
from molop.unit import atom_ureg
from molop.utils.types import (
    PintArray3,
    PintArray6,
    PintArray6Or3x3,
    PintArray10,
    PintArray15,
    SquareArray,
)


class Polarizability(BaseDataClassWithUnit):
    default_units: ClassVar[dict[str, UnitLike]] = {
        "electronic_spatial_extent": atom_ureg.bohr**2,
        "isotropic_polarizability": atom_ureg.bohr**3,
        "anisotropic_polarizability": atom_ureg.bohr**3,
        "polarizability_tensor": atom_ureg.bohr**3,
        "electric_dipole_moment": atom_ureg.debye,
        "dipole": atom_ureg.debye,
        "quadrupole": atom_ureg.debye * atom_ureg.angstrom,
        "octapole": atom_ureg.debye * atom_ureg.angstrom**2,
        "hexadecapole": atom_ureg.debye * atom_ureg.angstrom**3,
    }
    set_default_units: ClassVar[bool] = True

    # polarizability
    electronic_spatial_extent: PlainQuantity | None = Field(
        default=None,
        description="Electronic spatial extent, unit is bohr^2",
        exclude_if=lambda x: x is None,
    )
    isotropic_polarizability: PlainQuantity | None = Field(
        default=None,
        description="Isotropic polarizability, unit is bohr^3",
        exclude_if=lambda x: x is None,
    )
    anisotropic_polarizability: PlainQuantity | None = Field(
        default=None,
        description="Anisotropic polarizability, unit is bohr^3",
        exclude_if=lambda x: x is None,
    )
    polarizability_tensor: PintArray6Or3x3 | None = Field(
        default=None,
        description="Polarizability tensor",
        exclude_if=lambda x: (x is None) or (len(x) == 0),
    )
    electric_dipole_moment: PintArray3 | None = Field(
        default=None,
        description="Electric dipole moment, unit is `debye`",
        exclude_if=lambda x: (x is None) or (len(x) == 0),
    )
    dipole: PintArray3 | None = Field(
        default=None,
        description="Dipole moment, unit is `debye`",
        exclude_if=lambda x: (x is None) or (len(x) == 0),
    )
    quadrupole: PintArray6 | None = Field(
        default=None,
        description="Quadrupole moment, unit is `debye*angstrom`",
        exclude_if=lambda x: (x is None) or (len(x) == 0),
    )
    traceless_quadrupole: PintArray6 | None = Field(
        default=None,
        description="Traceless quadrupole moment, unit is `debye*angstrom`",
        exclude_if=lambda x: (x is None) or (len(x) == 0),
    )
    octapole: PintArray10 | None = Field(
        default=None,
        description="Octapole moment, unit is `debye*angstrom**2`",
        exclude_if=lambda x: (x is None) or (len(x) == 0),
    )
    hexadecapole: PintArray15 | None = Field(
        default=None,
        description="Hexadecapole moment, unit is `debye*angstrom**3`",
        exclude_if=lambda x: (x is None) or (len(x) == 0),
    )

    def to_summary_dict(self, **kwargs) -> SummaryDict:
        return summary_dict_from_fields(self, "Polarizability", **kwargs)

    @staticmethod
    def _has_quantity_payload(value: PlainQuantity | NumpyQuantity | None) -> bool:
        if value is None:
            return False
        magnitude = value.magnitude
        if not isinstance(magnitude, np.ndarray):
            return True
        return len(magnitude) > 0

    def to_property_bundle(self) -> PropertyBundle:
        scalar_properties: dict[str, PropertyScalarValue] = {
            key: value
            for key, value in {
                "electronic_spatial_extent": self.electronic_spatial_extent,
                "isotropic_polarizability": self.isotropic_polarizability,
                "anisotropic_polarizability": self.anisotropic_polarizability,
            }.items()
            if value is not None
        }
        vector_properties: dict[str, NumpyQuantity] = {}
        for key in ("electric_dipole_moment", "dipole"):
            value = getattr(self, key)
            if self._has_quantity_payload(value):
                vector_properties[key] = cast(NumpyQuantity, value)

        tensor_properties: dict[str, TensorProperty] = {}
        for key in (
            "polarizability_tensor",
            "quadrupole",
            "traceless_quadrupole",
            "octapole",
            "hexadecapole",
        ):
            value = getattr(self, key)
            if self._has_quantity_payload(value):
                tensor_properties[key] = TensorProperty(
                    label=key,
                    tensor=cast(NumpyQuantity, value),
                    frame="cartesian",
                )

        return PropertyBundle(
            scalar_properties=scalar_properties,
            vector_properties=vector_properties,
            tensor_properties=tensor_properties,
            metadata={"source": "Polarizability"},
        )


class BondOrders(BaseDataClassWithUnit):
    wiberg_bond_order: SquareArray = Field(
        default=np.zeros((0, 0)),
        description="Wiberg bond order",
        exclude_if=lambda x: x.shape == (0, 0),
    )
    mo_bond_order: SquareArray = Field(
        default=np.zeros((0, 0)),
        description="MO bond order, ∑[i∈A]∑[j∈B]P(i,j)",
        exclude_if=lambda x: x.shape == (0, 0),
    )
    mayer_bond_order: SquareArray = Field(
        default=np.zeros((0, 0)),
        description="MAYER POPULATION ANALYSIS bond order, ∑[i∈A]∑[j∈B]P(i,j)",
        exclude_if=lambda x: x.shape == (0, 0),
    )
    atom_atom_overlap_bond_order: SquareArray = Field(
        default=np.zeros((0, 0)),
        description="Atom-atom overlap bond order",
        exclude_if=lambda x: x.shape == (0, 0),
    )
    nbo_bond_order: SquareArray = Field(
        default=np.zeros((0, 0)),
        description="NBO bond order",
        exclude_if=lambda x: x.shape == (0, 0),
    )
    nbo_bond_order_for_alpha_spin: SquareArray = Field(
        default=np.zeros((0, 0)),
        description="NBO bond order for alpha spin",
        exclude_if=lambda x: x.shape == (0, 0),
    )
    nbo_bond_order_for_beta_spin: SquareArray = Field(
        default=np.zeros((0, 0)),
        description="NBO bond order for beta spin",
        exclude_if=lambda x: x.shape == (0, 0),
    )

    def to_summary_dict(self, **kwargs) -> SummaryDict:
        return {}

    @staticmethod
    def _has_matrix_payload(value: np.ndarray) -> bool:
        return value.shape != (0, 0)

    def to_property_bundle(self) -> PropertyBundle:
        tensor_properties: dict[str, TensorProperty] = {
            name: TensorProperty(label=name, tensor=value, frame="atom-pair")
            for name, value in {
                "wiberg_bond_order": self.wiberg_bond_order,
                "mo_bond_order": self.mo_bond_order,
                "mayer_bond_order": self.mayer_bond_order,
                "atom_atom_overlap_bond_order": self.atom_atom_overlap_bond_order,
                "nbo_bond_order": self.nbo_bond_order,
                "nbo_bond_order_for_alpha_spin": self.nbo_bond_order_for_alpha_spin,
                "nbo_bond_order_for_beta_spin": self.nbo_bond_order_for_beta_spin,
            }.items()
            if self._has_matrix_payload(value)
        }
        return PropertyBundle(
            tensor_properties=tensor_properties,
            metadata={"source": "BondOrders"},
        )


class Dispersions(BaseDataClassWithUnit):
    default_units: ClassVar[dict[str, UnitLike]] = {
        "C6AA": atom_ureg.bohr**6,
        "C8AA": atom_ureg.bohr**8,
    }
    set_default_units: ClassVar[bool] = True

    C6AA: PlainQuantity | None = Field(
        default=None,
        description="Mol. C6AA dispersion, unit is `bohr^6`",
        exclude_if=lambda x: x is None,
    )
    C8AA: PlainQuantity | None = Field(
        default=None,
        description="Mol. C8AA dispersion, unit is `bohr^8`",
        exclude_if=lambda x: x is None,
    )

    def to_summary_dict(self, **kwargs) -> SummaryDict:
        return summary_dict_from_fields(self, "Dispersion", **kwargs)

    def to_property_bundle(self) -> PropertyBundle:
        scalar_properties: dict[str, PropertyScalarValue] = {
            key: value
            for key, value in {
                "C6AA": self.C6AA,
                "C8AA": self.C8AA,
            }.items()
            if value is not None
        }
        return PropertyBundle(
            scalar_properties=scalar_properties,
            metadata={"source": "Dispersions"},
        )


class SinglePointProperties(BaseDataClassWithUnit):
    """
    Single point properties.
    """

    default_units: ClassVar[dict[str, UnitLike]] = {
        "vip": atom_ureg.Unit("eV / particle"),
        "vea": atom_ureg.Unit("eV / particle"),
        "gei": atom_ureg.Unit("eV / particle"),
    }
    set_default_units: ClassVar[bool] = True

    vip: PlainQuantity | None = Field(
        default=None,
        description="Vertical ionization potential, unit is `eV/particle`",
        exclude_if=lambda x: x is None,
    )
    vea: PlainQuantity | None = Field(
        default=None,
        description="Vertical electron affinity, unit is `eV/particle`",
        exclude_if=lambda x: x is None,
    )
    gei: PlainQuantity | None = Field(
        default=None,
        description="Global Electrophilicity Index, unit is `eV/particle`",
        exclude_if=lambda x: x is None,
    )
    fukui_positive: list[float] = Field(
        default_factory=list, description="Fukui Index f(+)", exclude_if=lambda x: len(x) == 0
    )
    fukui_negative: list[float] = Field(
        default_factory=list, description="Fukui Index f(-)", exclude_if=lambda x: len(x) == 0
    )
    fukui_zero: list[float] = Field(
        default_factory=list, description="Fukui Index f(0)", exclude_if=lambda x: len(x) == 0
    )
    fod: list[float] = Field(
        default_factory=list, description="fractional occupation density population"
    )

    @property
    def atomic_property_names(self) -> list[str]:
        return [
            name
            for name in ("fukui_positive", "fukui_negative", "fukui_zero", "fod")
            if len(getattr(self, name)) > 0
        ]

    def to_atomic_property_table(self, atom_symbols: Sequence[str] | None = None) -> PropertyTable:
        property_lengths = [len(getattr(self, name)) for name in self.atomic_property_names]
        if property_lengths:
            num_atoms = property_lengths[0]
        elif atom_symbols is not None:
            num_atoms = len(atom_symbols)
        else:
            num_atoms = 0

        if atom_symbols is not None and len(atom_symbols) != num_atoms:
            raise ValueError(
                f"atom_symbols length {len(atom_symbols)} does not match atomic property length "
                f"{num_atoms}"
            )

        columns: dict[str, PropertyColumnValue] = {}
        if num_atoms:
            columns["atom_index"] = np.arange(num_atoms)
        if atom_symbols is not None:
            columns["atom_symbol"] = list(atom_symbols)
        for property_name in self.atomic_property_names:
            columns[property_name] = getattr(self, property_name)

        return PropertyTable(
            columns=columns,
            row_labels=[str(i) for i in range(num_atoms)],
            metadata={
                "source": "SinglePointProperties",
                "index_base": 0,
                "atomic_property_names": self.atomic_property_names,
            },
        )

    def to_property_bundle(self, atom_symbols: Sequence[str] | None = None) -> PropertyBundle:
        scalar_properties: dict[str, PropertyScalarValue] = {
            key: value
            for key, value in {
                "vip": self.vip,
                "vea": self.vea,
                "gei": self.gei,
            }.items()
            if value is not None
        }
        table = self.to_atomic_property_table(atom_symbols=atom_symbols)
        return PropertyBundle(
            scalar_properties=scalar_properties,
            tables={"atomic_properties": table} if len(table) else {},
            metadata={"source": "SinglePointProperties"},
        )

    def to_summary_dict(self, **kwargs) -> SummaryDict:
        return {}


__all__ = ["BondOrders", "Dispersions", "Polarizability", "SinglePointProperties"]
