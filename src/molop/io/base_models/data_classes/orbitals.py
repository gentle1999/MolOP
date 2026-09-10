"""Molecular and natural orbital data classes."""

from __future__ import annotations

from collections.abc import Iterator, Sequence
from typing import Any, ClassVar, cast, overload

import numpy as np
from pint._typing import UnitLike
from pint.facets.plain import PlainQuantity
from pydantic import Field, computed_field, model_validator
from typing_extensions import Self

from molop.io.base_models.Bases import (
    BaseDataClassWithUnit,
    PropertyBundle,
    PropertyScalarValue,
    PropertyTable,
)
from molop.io.base_models.summary import SummaryDict, summary_dict_from_fields, summary_item
from molop.unit import atom_ureg
from molop.utils.types import ArrayN, PintArrayN


class MoleculeOrbital(BaseDataClassWithUnit):
    default_units: ClassVar[dict[str, UnitLike]] = {
        "alpha_energy": atom_ureg.hartree,
        "beta_energy": atom_ureg.hartree,
    }
    set_default_units: ClassVar[bool] = True

    alpha_energy: PlainQuantity | None = Field(
        default=None, description="alpha orbital energy, unit is `hartree`"
    )
    beta_energy: PlainQuantity | None = Field(
        default=None, description="beta orbital energy, unit is `hartree`"
    )
    alpha_occupancy: float | bool | None = Field(
        default=None, description="alpha orbital occupancy"
    )
    alpha_symmetry: str | None = Field(default=None, description="alpha orbital symmetry")
    beta_occupancy: float | bool | None = Field(default=None, description="beta orbital occupancy")
    beta_symmetry: str | None = Field(default=None, description="beta orbital symmetry")
    coefficient: ArrayN | None = Field(default=None, description="coefficient of the orbital")

    def to_summary_dict(self, **kwargs) -> SummaryDict:
        return summary_dict_from_fields(self, "Orbital", **kwargs)


class MolecularOrbitals(BaseDataClassWithUnit, Sequence[MoleculeOrbital]):
    default_units: ClassVar[dict[str, UnitLike]] = {
        "alpha_energies": atom_ureg.hartree,
        "beta_energies": atom_ureg.hartree,
    }
    set_default_units: ClassVar[bool] = True

    # orbital energies
    electronic_state: str | None = Field(
        default=None, description="electronic state of the molecule"
    )
    alpha_energies: PintArrayN = Field(
        default=np.array([]) * atom_ureg.hartree,
        description="alpha orbital energies, unit is `hartree`",
    )
    beta_energies: PintArrayN = Field(
        default=np.array([]) * atom_ureg.hartree,
        description="beta orbital energies, unit is `hartree`",
    )
    alpha_occupancies: list[float | bool | None] = Field(
        default_factory=list, description="alpha orbital occupancies"
    )
    beta_occupancies: list[float | bool | None] = Field(
        default_factory=list, description="beta orbital occupancies"
    )
    alpha_symmetries: list[str | None] = Field(
        default_factory=list, description="alpha orbital symmetries"
    )
    beta_symmetries: list[str | None] = Field(
        default_factory=list, description="beta orbital symmetries"
    )
    coefficients: list[ArrayN | None] = Field(
        default_factory=list, description="coefficients of the orbitals"
    )

    @computed_field(description="HOMO orbital idx")  # type: ignore[prop-decorator]
    @property
    def HOMO_id(self) -> int | None:
        return self._frontier_occupied_id(self.alpha_occupancies)

    @computed_field(description="LUMO orbital idx")  # type: ignore[prop-decorator]
    @property
    def LUMO_id(self) -> int | None:
        if self.HOMO_id is None:
            return None
        return self.HOMO_id + 1

    @computed_field(description="beta HOMO orbital idx")  # type: ignore[prop-decorator]
    @property
    def beta_HOMO_id(self) -> int | None:
        return self._frontier_occupied_id(self.beta_occupancies)

    @computed_field(description="beta LUMO orbital idx")  # type: ignore[prop-decorator]
    @property
    def beta_LUMO_id(self) -> int | None:
        if self.beta_HOMO_id is None:
            return None
        return self.beta_HOMO_id + 1

    @computed_field(description="SOMO orbital idx")  # type: ignore[prop-decorator]
    @property
    def SOMO_ids(self) -> list[int]:
        if len(self.beta_occupancies) == 0:
            return []
        return [
            i
            for i, (alpha_occ, beta_occ) in enumerate(
                zip(self.alpha_occupancies, self.beta_occupancies, strict=True)
            )
            if self._occupancy_value(alpha_occ) > 0.0 and self._occupancy_value(beta_occ) <= 0.0
        ]

    @computed_field(description="NHOMO orbital idx")  # type: ignore[prop-decorator]
    @property
    def NHOMO_id(self) -> int | None:
        if self.HOMO_id is None:
            return None
        if self.HOMO_id == 0:
            return None
        return self.HOMO_id - 1

    @computed_field(description="SLUMO orbital idx")  # type: ignore[prop-decorator]
    @property
    def SLUMO_id(self) -> int | None:
        if self.LUMO_id is None:
            return None
        if self.LUMO_id == len(self.alpha_occupancies):
            return None
        return self.LUMO_id + 1

    @computed_field(description="HOMO energy")  # type: ignore[prop-decorator]
    @property
    def HOMO_energy(self) -> PlainQuantity | None:
        if self.HOMO_id is None:
            return None
        if len(self.alpha_energies) <= self.HOMO_id:
            return None
        return self.alpha_energies[self.HOMO_id]

    @computed_field(description="LUMO energy")  # type: ignore[prop-decorator]
    @property
    def LUMO_energy(self) -> PlainQuantity | None:
        if self.LUMO_id is None:
            return None
        if len(self.alpha_energies) <= self.LUMO_id:
            return None
        return self.alpha_energies[self.LUMO_id]

    @computed_field(description="NHOMO energy")  # type: ignore[prop-decorator]
    @property
    def NHOMO_energy(self) -> PlainQuantity | None:
        if self.NHOMO_id is None:
            return None
        if len(self.alpha_energies) <= self.NHOMO_id:
            return None
        return self.alpha_energies[self.NHOMO_id]

    @computed_field(description="SLUMO energy")  # type: ignore[prop-decorator]
    @property
    def SLUMO_energy(self) -> PlainQuantity | None:
        if self.SLUMO_id is None:
            return None
        if len(self.alpha_energies) <= self.SLUMO_id:
            return None
        return self.alpha_energies[self.SLUMO_id]

    @computed_field(description="HOMO-LUMO gap")  # type: ignore[prop-decorator]
    @property
    def HOMO_LUMO_gap(self) -> PlainQuantity | None:
        if self.HOMO_energy is None or self.LUMO_energy is None:
            return None
        return cast(Any, self.LUMO_energy) - cast(Any, self.HOMO_energy)

    def __iter__(self) -> Iterator[MoleculeOrbital]:  # type: ignore[override]
        for i in range(len(self)):
            yield self[i]

    @computed_field(description="beta HOMO energy")  # type: ignore[prop-decorator]
    @property
    def beta_HOMO_energy(self) -> PlainQuantity | None:
        if self.beta_HOMO_id is None:
            return None
        if len(self.beta_energies) <= self.beta_HOMO_id:
            return None
        return self.beta_energies[self.beta_HOMO_id]

    @computed_field(description="beta LUMO energy")  # type: ignore[prop-decorator]
    @property
    def beta_LUMO_energy(self) -> PlainQuantity | None:
        if self.beta_LUMO_id is None:
            return None
        if len(self.beta_energies) <= self.beta_LUMO_id:
            return None
        return self.beta_energies[self.beta_LUMO_id]

    @computed_field(description="beta HOMO-LUMO gap")  # type: ignore[prop-decorator]
    @property
    def beta_HOMO_LUMO_gap(self) -> PlainQuantity | None:
        if self.beta_HOMO_energy is None or self.beta_LUMO_energy is None:
            return None
        return cast(Any, self.beta_LUMO_energy) - cast(Any, self.beta_HOMO_energy)

    @computed_field(description="beta NHOMO orbital idx")  # type: ignore[prop-decorator]
    @property
    def beta_NHOMO_id(self) -> int | None:
        if self.beta_HOMO_id is None:
            return None
        if self.beta_HOMO_id == 0:
            return None
        return self.beta_HOMO_id - 1

    @computed_field(description="beta NHOMO energy")  # type: ignore[prop-decorator]
    @property
    def beta_NHOMO_energy(self) -> PlainQuantity | None:
        if self.beta_NHOMO_id is None:
            return None
        if len(self.beta_energies) <= self.beta_NHOMO_id:
            return None
        return self.beta_energies[self.beta_NHOMO_id]

    @computed_field(description="beta SLUMO orbital idx")  # type: ignore[prop-decorator]
    @property
    def beta_SLUMO_id(self) -> int | None:
        if self.beta_LUMO_id is None:
            return None
        if self.beta_LUMO_id == len(self.beta_occupancies):
            return None
        return self.beta_LUMO_id + 1

    @computed_field(description="beta SLUMO energy")  # type: ignore[prop-decorator]
    @property
    def beta_SLUMO_energy(self) -> PlainQuantity | None:
        if self.beta_SLUMO_id is None:
            return None
        if len(self.beta_energies) <= self.beta_SLUMO_id:
            return None
        return self.beta_energies[self.beta_SLUMO_id]

    @overload
    def __getitem__(self, orbitalIDX: int) -> MoleculeOrbital: ...
    @overload
    def __getitem__(self, orbitalIDX: slice) -> list[MoleculeOrbital]: ...
    @overload
    def __getitem__(self, orbitalIDX: Sequence) -> list[MoleculeOrbital]: ...
    def __getitem__(
        self, orbitalIDX: int | slice | Sequence
    ) -> MoleculeOrbital | list[MoleculeOrbital]:
        def get_item(seq: Any, idx: int):
            if len(seq) > abs(idx):
                return seq[idx]
            else:
                return None

        if isinstance(orbitalIDX, int):
            assert max(
                (
                    len(self.alpha_energies),
                    len(self.beta_energies),
                    len(self.alpha_occupancies),
                    len(self.beta_occupancies),
                    len(self.alpha_symmetries),
                    len(self.beta_symmetries),
                    len(self.coefficients),
                )
            ) > abs(orbitalIDX), f"orbital index {orbitalIDX} out of range"
            return MoleculeOrbital.model_validate(
                {
                    "alpha_energy": get_item(cast(Sequence[Any], self.alpha_energies), orbitalIDX),
                    "beta_energy": get_item(cast(Sequence[Any], self.beta_energies), orbitalIDX),
                    "alpha_occupancy": get_item(self.alpha_occupancies, orbitalIDX),
                    "beta_occupancy": get_item(self.beta_occupancies, orbitalIDX),
                    "alpha_symmetry": get_item(self.alpha_symmetries, orbitalIDX),
                    "beta_symmetry": get_item(self.beta_symmetries, orbitalIDX),
                    "coefficient": get_item(self.coefficients, orbitalIDX),
                }
            )
        if isinstance(orbitalIDX, slice):
            return [
                self[orbital_id]
                for orbital_id in range(*orbitalIDX.indices(len(self.alpha_energies)))
            ]
        else:
            return [self[orbital_id] for orbital_id in orbitalIDX]

    def __len__(self) -> int:
        return len(self.alpha_energies)

    @staticmethod
    def _occupancy_value(occupancy: float | bool | None) -> float:
        if occupancy is None:
            return 0.0
        if isinstance(occupancy, bool):
            return 1.0 if occupancy else 0.0
        return float(occupancy)

    @classmethod
    def _frontier_occupied_id(cls, occupancies: Sequence[float | bool | None]) -> int | None:
        last_occupied: int | None = None
        for idx, occupancy in enumerate(occupancies):
            value = cls._occupancy_value(occupancy)
            if value > 0.0:
                last_occupied = idx
                continue
            if last_occupied is not None:
                return last_occupied
        return last_occupied

    @property
    def HOMO(self) -> MoleculeOrbital | None:
        if self.HOMO_id is None:
            return None
        return self[self.HOMO_id]

    @property
    def LUMO(self) -> MoleculeOrbital | None:
        if self.LUMO_id is None:
            return None
        return self[self.LUMO_id]

    @property
    def beta_HOMO(self) -> MoleculeOrbital | None:
        if self.beta_HOMO_id is None:
            return None
        return self[self.beta_HOMO_id]

    @property
    def beta_LUMO(self) -> MoleculeOrbital | None:
        if self.beta_LUMO_id is None:
            return None
        return self[self.beta_LUMO_id]

    @property
    def SOMOs(self) -> list[MoleculeOrbital]:
        return self[self.SOMO_ids]

    def to_orbital_table(self, spin: str = "alpha") -> PropertyTable:
        spin_key = spin.lower()
        if spin_key not in {"alpha", "beta"}:
            raise ValueError(f"Unsupported spin channel {spin!r}; expected 'alpha' or 'beta'")

        energies = self.alpha_energies if spin_key == "alpha" else self.beta_energies
        occupancies = self.alpha_occupancies if spin_key == "alpha" else self.beta_occupancies
        symmetries = self.alpha_symmetries if spin_key == "alpha" else self.beta_symmetries
        num_orbitals = len(energies)

        columns = {
            "orbital_index": np.arange(num_orbitals),
            "energy": energies,
        }
        if occupancies:
            columns["occupancy"] = occupancies
        if symmetries:
            columns["symmetry"] = symmetries

        return PropertyTable(
            columns=columns,
            row_labels=[str(i) for i in range(num_orbitals)],
            metadata={
                "source": "MolecularOrbitals",
                "spin": spin_key,
                "index_base": 0,
                "electronic_state": self.electronic_state,
                "has_coefficients": bool(self.coefficients),
            },
        )

    def to_property_bundle(self) -> PropertyBundle:
        scalar_properties: dict[str, PropertyScalarValue] = {
            key: value
            for key, value in {
                "HOMO_id": self.HOMO_id,
                "LUMO_id": self.LUMO_id,
                "HOMO_energy": self.HOMO_energy,
                "LUMO_energy": self.LUMO_energy,
                "HOMO_LUMO_gap": self.HOMO_LUMO_gap,
                "beta_HOMO_id": self.beta_HOMO_id,
                "beta_LUMO_id": self.beta_LUMO_id,
                "beta_HOMO_energy": self.beta_HOMO_energy,
                "beta_LUMO_energy": self.beta_LUMO_energy,
                "beta_HOMO_LUMO_gap": self.beta_HOMO_LUMO_gap,
            }.items()
            if value is not None
        }
        tables: dict[str, PropertyTable] = {}
        if len(self.alpha_energies):
            tables["alpha_orbitals"] = self.to_orbital_table("alpha")
        if len(self.beta_energies):
            tables["beta_orbitals"] = self.to_orbital_table("beta")
        return PropertyBundle(
            scalar_properties=scalar_properties,
            tables=tables,
            metadata={"source": "MolecularOrbitals"},
        )

    @model_validator(mode="after")
    def validate_molecular_orbitals(self) -> Self:
        assert len(self.alpha_energies) == len(self.alpha_occupancies), (
            "alpha orbital energies and occupancies must have the same length"
        )
        assert len(self.beta_energies) == len(self.beta_occupancies), (
            "beta orbital energies and occupancies must have the same length"
        )
        if self.alpha_symmetries:
            assert len(self.alpha_symmetries) == len(self.alpha_energies), (
                "alpha orbital symmetries and energies must have the same length"
            )
        if self.beta_symmetries:
            assert len(self.beta_symmetries) == len(self.beta_energies), (
                "beta orbital symmetries and energies must have the same length"
            )
        return self

    def to_summary_dict(self, **kwargs) -> SummaryDict:
        summary: SummaryDict = {}
        for field, value in {
            "electronic_state": self.electronic_state,
            "HOMO_energy": self.HOMO_energy,
            "LUMO_energy": self.LUMO_energy,
            "HOMO-LUMO_gap": self.HOMO_LUMO_gap,
        }.items():
            item = summary_item("Orbitals", field, value)
            if item is not None:
                column, magnitude = item
                summary[column] = magnitude
        return summary


class NaturalAtomicOrbital(BaseDataClassWithUnit):
    default_units: ClassVar[dict[str, UnitLike]] = {"energy": atom_ureg.Unit("hartree")}
    set_default_units: ClassVar[bool] = True

    orbital_index: int = Field(
        default=0,
        description="index of the orbital",
    )
    element: str = Field(
        default="",
        description="The atomic elements to which the natural bond orbitals belong",
        exclude_if=lambda x: x == "",
    )
    atom_index: int = Field(
        default=0,
        description="The index of the atom to which the natural bond orbital belongs",
        exclude_if=lambda x: x == 0,
    )
    angular_momentum: str = Field(
        default="",
        description="The angular momentum of the natural bond orbital",
        exclude_if=lambda x: x == "",
    )
    ao_type: str = Field(
        default="",
        description="The type of the atomic orbital to which the natural bond orbital belongs",
        exclude_if=lambda x: x == "",
    )
    occupancy: float = Field(
        default=0.0,
        description="The occupancy of the natural bond orbital",
        exclude_if=lambda x: x == 0.0,
    )
    energy: PlainQuantity | None = Field(
        default=None,
        description="The energy of the natural bond orbital",
        exclude_if=lambda x: x is None,
    )

    def to_summary_dict(self, **kwargs) -> SummaryDict:
        return summary_dict_from_fields(self, "NaturalAtomicOrbital", **kwargs)


class NaturalAtomicOrbitals(BaseDataClassWithUnit):
    orbitals: list[NaturalAtomicOrbital] = Field(
        default_factory=list,
        description="The list of natural atomic orbitals",
        exclude_if=lambda x: len(x) == 0,
    )
    alpha_spin_orbitals: list[NaturalAtomicOrbital] = Field(
        default_factory=list,
        description="The list of alpha spin natural atomic orbitals",
        exclude_if=lambda x: len(x) == 0,
    )
    beta_spin_orbitals: list[NaturalAtomicOrbital] = Field(
        default_factory=list,
        description="The list of beta spin natural atomic orbitals",
        exclude_if=lambda x: len(x) == 0,
    )

    def to_summary_dict(self, **kwargs) -> SummaryDict:
        return {}


# TODO: ready for NaturalBondOrbitals
class NaturalBondOrbital(BaseDataClassWithUnit):
    default_units: ClassVar[dict[str, UnitLike]] = {"energy": atom_ureg.Unit("hartree")}
    set_default_units: ClassVar[bool] = True

    orbital_index: int = Field(
        default=0,
        description="index of the orbital",
    )
    orbital_type: str = Field(
        default="",
        description="The type of the natural bond orbital",
        exclude_if=lambda x: x == "",
    )
    sub_index: int = Field(
        default=0,
        description="The sub-index of the natural bond orbital",
        exclude_if=lambda x: x == 0,
    )
    bonding_atoms_dict: dict[str, int] = Field(
        default_factory=dict,
        description="The dictionary of bonding atoms, key is the atom element, value is the atom index (1-based)",
        exclude_if=lambda x: len(x) == 0,
    )
    occupancy: float = Field(
        default=0.0,
        description="The occupancy of the natural bond orbital",
        exclude_if=lambda x: x == 0.0,
    )
    energy: PlainQuantity | None = Field(
        default=None,
        description="The energy of the natural bond orbital",
        exclude_if=lambda x: x is None,
    )
    principal_elocalizations_dict: dict[str, float] = Field(
        default_factory=dict,
        description="The dictionary of principal delocalizations, key is the delocalization type "
        "(geminal,vicinal,remote), value is the bond orbital index",
        exclude_if=lambda x: len(x) == 0,
    )

    def to_summary_dict(self, **kwargs) -> SummaryDict:
        return summary_dict_from_fields(self, "NaturalBondOrbital", **kwargs)


class NaturalBondOrbitals(BaseDataClassWithUnit):
    orbitals: list[NaturalBondOrbital] = Field(
        default_factory=list,
        description="The list of natural bond orbitals",
        exclude_if=lambda x: len(x) == 0,
    )
    alpha_spin_orbitals: list[NaturalBondOrbital] = Field(
        default_factory=list,
        description="The list of alpha spin natural bond orbitals",
        exclude_if=lambda x: len(x) == 0,
    )
    beta_spin_orbitals: list[NaturalBondOrbital] = Field(
        default_factory=list,
        description="The list of beta spin natural bond orbitals",
        exclude_if=lambda x: len(x) == 0,
    )

    def to_summary_dict(self, **kwargs) -> SummaryDict:
        return {}


__all__ = [
    "MolecularOrbitals",
    "MoleculeOrbital",
    "NaturalAtomicOrbital",
    "NaturalAtomicOrbitals",
    "NaturalBondOrbital",
    "NaturalBondOrbitals",
]
