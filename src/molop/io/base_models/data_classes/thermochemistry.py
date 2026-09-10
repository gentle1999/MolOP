"""Energy and thermochemistry data classes."""

from __future__ import annotations

from typing import ClassVar, Literal

from pint._typing import UnitLike
from pint.facets.plain import PlainQuantity
from pydantic import ConfigDict, Field, computed_field, model_validator
from typing_extensions import Self

from molop.io.base_models.Bases import BaseDataClassWithUnit, PropertyBundle, PropertyScalarValue
from molop.io.base_models.summary import SummaryDict, summary_dict_from_fields
from molop.unit import atom_ureg
from molop.utils.types import PintArrayN


class EnergyObservation(BaseDataClassWithUnit):
    """One source-labeled energy observation retained by the core energy model."""

    default_units: ClassVar[dict[str, UnitLike]] = {"value": atom_ureg.hartree}
    set_default_units: ClassVar[bool] = True

    method: str = Field(min_length=1, description="Energy method or reference label")
    quantity_semantics: Literal[
        "total_energy",
        "correlation_correction",
        "component",
    ] = Field(description="Physical meaning of the observed energy")
    value: PlainQuantity = Field(description="Observed energy, normalized to hartree")
    source_label: str = Field(min_length=1, description="Label printed by the source program")


class Energies(BaseDataClassWithUnit):
    model_config = ConfigDict(arbitrary_types_allowed=True, extra="forbid")
    default_units: ClassVar[dict[str, UnitLike]] = {
        "electronic_energy": atom_ureg.hartree,
        "reference_energy": atom_ureg.hartree,
        "mp2_energy": atom_ureg.hartree,
        "mp3_energy": atom_ureg.hartree,
        "mp4_energy": atom_ureg.hartree,
        "mp5_energy": atom_ureg.hartree,
        "ccsd_energy": atom_ureg.hartree,
        "ccsd_t_energy": atom_ureg.hartree,
    }
    set_default_units: ClassVar[bool] = True

    observations: list[EnergyObservation] = Field(
        default_factory=list,
        description="Optional source-labeled energy observations captured during parsing",
        exclude_if=lambda observations: not observations,
    )

    # energies
    electronic_energy: PlainQuantity | None = Field(
        default=None,
        description="Electronic energy of the molecule, unit is `hartree`",
    )
    reference_energy: PlainQuantity | None = Field(
        default=None,
        description="Reference electronic energy such as an HF or Kohn-Sham energy, unit is `hartree`",
    )
    mp2_energy: PlainQuantity | None = Field(
        default=None, description="MP2 energy of the molecule, unit is `hartree`"
    )
    mp3_energy: PlainQuantity | None = Field(
        default=None, description="MP3 energy of the molecule, unit is `hartree`"
    )
    mp4_energy: PlainQuantity | None = Field(
        default=None, description="MP4 energy of the molecule, unit is `hartree`"
    )
    mp5_energy: PlainQuantity | None = Field(
        default=None, description="MP5 energy of the molecule, unit is `hartree`"
    )
    ccsd_energy: PlainQuantity | None = Field(
        default=None, description="CCSD energy of the molecule, unit is `hartree`"
    )
    ccsd_t_energy: PlainQuantity | None = Field(
        default=None, description="CCSD(T) energy of the molecule, unit is `hartree`"
    )

    @property
    def energy(self) -> dict[str, PlainQuantity]:
        energy_fields = (
            "electronic_energy",
            "ccsd_t_energy",
            "ccsd_energy",
            "mp5_energy",
            "mp4_energy",
            "mp3_energy",
            "mp2_energy",
            "reference_energy",
        )
        return {
            energy_field: getattr(self, energy_field)
            for energy_field in energy_fields
            if getattr(self, energy_field) is not None
        }

    @computed_field(description="Total energy, unit is `hartree`")  # type: ignore[prop-decorator]
    @property
    def total_energy(self) -> PlainQuantity | None:
        keys = list(self.energy.keys())
        if len(keys) > 0:
            return self.energy[keys[0]].to(atom_ureg.hartree)
        else:
            return None

    def to_summary_dict(self, **kwargs) -> SummaryDict:
        return summary_dict_from_fields(self, "Energy", **kwargs)

    def to_property_bundle(self) -> PropertyBundle:
        scalar_properties: dict[str, PropertyScalarValue] = {
            key: value
            for key, value in {
                **self.energy,
                "total_energy": self.total_energy,
            }.items()
            if value is not None
        }
        return PropertyBundle(
            scalar_properties=scalar_properties,
            metadata={"source": "Energies"},
        )


class ThermalInformations(BaseDataClassWithUnit):
    """
    Thermal energy data

    ref: https://www.cup.uni-muenchen.de/ch/compchem/vib/thermo1.html
    Theorietically, the energy follow the relationship below:
    U_0 = E_tot + ZPVE
    U_T(?K) = E_tot + TCE
    H_T(?K) = E_tot + TCH
    G_T(?K) = E_tot + TCG
    G_T(?K) = H_T(?K) - T * S(?K)
    """

    default_units: ClassVar[dict[str, UnitLike]] = {
        "ZPVE": atom_ureg.kcal / atom_ureg.mol,
        "U_0": atom_ureg.kcal / atom_ureg.mol,
        "TCE": atom_ureg.kcal / atom_ureg.mol,
        "TCH": atom_ureg.kcal / atom_ureg.mol,
        "TCG": atom_ureg.kcal / atom_ureg.mol,
        "U_T": atom_ureg.kcal / atom_ureg.mol,
        "H_T": atom_ureg.kcal / atom_ureg.mol,
        "G_T": atom_ureg.kcal / atom_ureg.mol,
        "S": atom_ureg.calorie / atom_ureg.mol / atom_ureg.kelvin,
        "C_V": atom_ureg.calorie / atom_ureg.mol / atom_ureg.kelvin,
        "molecular_mass": atom_ureg.amu,
        "moments_of_inertia": atom_ureg.amu * atom_ureg.bohr**2,
        "rotational_temperatures": atom_ureg.K,
        "rotational_constants": atom_ureg.gigahertz,
        "vibrational_temperatures": atom_ureg.K,
    }
    set_default_units: ClassVar[bool] = True

    ZPVE: PlainQuantity | None = Field(
        default=None,
        description="Zero-point vibrational energy, unit is `kcal/mol`",
        exclude_if=lambda x: x is None,
    )
    TCE: PlainQuantity | None = Field(
        default=None,
        description="thermal correction to the internal energy at ?K, unit is `kcal/mol`",
        exclude_if=lambda x: x is None,
    )
    TCH: PlainQuantity | None = Field(
        default=None,
        description="thermal correction to the enthalpy at ?K, unit is `kcal/mol`",
        exclude_if=lambda x: x is None,
    )
    TCG: PlainQuantity | None = Field(
        default=None,
        description="thermal correction to the Gibbs free energy at ?K, unit is `kcal/mol`",
        exclude_if=lambda x: x is None,
    )
    U_0: PlainQuantity | None = Field(
        default=None,
        description="Zero-point energy, unit is `kcal/mol`",
        exclude_if=lambda x: x is None,
    )
    U_T: PlainQuantity | None = Field(
        default=None,
        description="thermal energy at ?K, unit is `kcal/mol`",
        exclude_if=lambda x: x is None,
    )
    H_T: PlainQuantity | None = Field(
        default=None,
        description="enthalpy at ?K, unit is `kcal/mol`",
        exclude_if=lambda x: x is None,
    )
    G_T: PlainQuantity | None = Field(
        default=None,
        description="Gibbs Free Energy at ?K, unit is `kcal/mol`",
        exclude_if=lambda x: x is None,
    )
    S: PlainQuantity | None = Field(
        default=None,
        description="entropy at ?K, unit is `cal/mol/K`",
        exclude_if=lambda x: x is None,
    )
    C_V: PlainQuantity | None = Field(
        default=None,
        description="heat capacity at constant volume, unit is `cal/mol/K`",
        exclude_if=lambda x: x is None,
    )
    molecular_mass: PlainQuantity | None = Field(
        default=None,
        description="molecular mass from the thermochemistry section, unit is `amu`",
        exclude_if=lambda x: x is None,
    )
    moments_of_inertia: PintArrayN | None = Field(
        default=None,
        description="principal moments of inertia, unit is `amu*bohr**2`",
        exclude_if=lambda x: x is None,
    )
    rotational_symmetry_number: int | None = Field(
        default=None,
        description="rotational symmetry number from the thermochemistry section",
        exclude_if=lambda x: x is None,
    )
    rotational_temperatures: PintArrayN | None = Field(
        default=None,
        description="rotational temperatures, unit is `K`",
        exclude_if=lambda x: x is None,
    )
    rotational_constants: PintArrayN | None = Field(
        default=None,
        description="rotational constants from frequency thermochemistry, unit is `GHz`",
        exclude_if=lambda x: x is None,
    )
    vibrational_temperatures: PintArrayN | None = Field(
        default=None,
        description="vibrational temperatures, unit is `K`",
        exclude_if=lambda x: x is None,
    )
    vibrational_temperature_mode_indices: list[int] | None = Field(
        default=None,
        description="Frequency-mode indices corresponding to vibrational_temperatures",
        exclude_if=lambda x: x is None,
    )

    @model_validator(mode="after")
    def validate_vibrational_temperature_mode_indices(self) -> Self:
        if self.vibrational_temperature_mode_indices is None:
            return self
        if self.vibrational_temperatures is None:
            raise ValueError(
                "vibrational_temperature_mode_indices require vibrational_temperatures"
            )
        if len(self.vibrational_temperature_mode_indices) != len(self.vibrational_temperatures):
            raise ValueError(
                "vibrational_temperature_mode_indices must match vibrational_temperatures"
            )
        if len(set(self.vibrational_temperature_mode_indices)) != len(
            self.vibrational_temperature_mode_indices
        ) or any(index < 0 for index in self.vibrational_temperature_mode_indices):
            raise ValueError(
                "vibrational_temperature_mode_indices must be unique non-negative indices"
            )
        return self

    def to_summary_dict(self, **kwargs) -> SummaryDict:
        return summary_dict_from_fields(self, "Thermal", **kwargs)

    def to_property_bundle(self) -> PropertyBundle:
        scalar_properties: dict[str, PropertyScalarValue] = {
            key: value
            for key, value in {
                "ZPVE": self.ZPVE,
                "TCE": self.TCE,
                "TCH": self.TCH,
                "TCG": self.TCG,
                "U_0": self.U_0,
                "U_T": self.U_T,
                "H_T": self.H_T,
                "G_T": self.G_T,
                "S": self.S,
                "C_V": self.C_V,
                "molecular_mass": self.molecular_mass,
                "rotational_symmetry_number": self.rotational_symmetry_number,
            }.items()
            if value is not None
        }
        return PropertyBundle(
            scalar_properties=scalar_properties,
            metadata={"source": "ThermalInformations"},
        )


__all__ = ["Energies", "EnergyObservation", "ThermalInformations"]
