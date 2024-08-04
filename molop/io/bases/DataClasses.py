"""
Author: TMJ
Date: 2024-06-18 20:42:42
LastEditors: TMJ
LastEditTime: 2024-06-18 20:43:03
Description: 请填写简介
"""

from typing import List, Sequence, Union, overload

import numpy as np
from pint.facets.plain import PlainQuantity
from pydantic import BaseModel, ConfigDict, Field, PrivateAttr, computed_field

from molop.unit import atom_ureg


class BaseDataClassWithUnit(BaseModel):
    model_config = ConfigDict(arbitrary_types_allowed=True)

    # @model_serializer
    # def serialize_model(self):
    #    """
    #    Serialize the model to a dictionary.
    #    If the field is a Quantity, serialize it as a dictionary without units.
    #    """
    #    fields = {
    #        k: (
    #            getattr(self, k).magnitude
    #            if isinstance(getattr(self, k), PlainQuantity)
    #            else getattr(self, k)
    #        )
    #        for k in self.model_computed_fields
    #    }
    #    computed_fields = {
    #        k: (
    #            getattr(self, k).magnitude
    #            if isinstance(getattr(self, k), PlainQuantity)
    #            else getattr(self, k)
    #        )
    #        for k in self.model_fields
    #    }
    #    return {**fields, **computed_fields}

    def to_unitless_dump(self, **kwargs):
        return {
            k: (
                (v.m.tolist() if isinstance(v.m, np.ndarray) else v.m)
                if isinstance(v, PlainQuantity)
                else (
                    getattr(self, k).to_unitless_dump(**kwargs)
                    if hasattr(getattr(self, k), "to_unitless_dump")
                    else (
                        [
                            (
                                v1.to_unitless_dump(**kwargs)
                                if hasattr(v1, "to_unitless_dump")
                                else v1
                            )
                            for v1 in getattr(self, k)
                        ]
                        if isinstance(v, List)
                        else v.tolist() if isinstance(v, np.ndarray) else v
                    )
                )
            )
            for k, v in self.model_dump(**kwargs).items()
        }


class WaveFunction(BaseDataClassWithUnit): ...


class Energies(BaseDataClassWithUnit):
    # energies
    electronic_energy: Union[PlainQuantity, None] = Field(
        default=None,
        description="Electronic energy of the molecule, unit is `hartree/particle`",
    )
    scf_energy: Union[PlainQuantity, None] = Field(
        default=None,
        description="SCF energy of the molecule, unit is `hartree/particle`",
    )
    mp2_energy: Union[PlainQuantity, None] = Field(
        default=None,
        description="MP2 energy of the molecule, unit is `hartree/particle`",
    )
    mp3_energy: Union[PlainQuantity, None] = Field(
        default=None,
        description="MP3 energy of the molecule, unit is `hartree/particle`",
    )
    mp4_energy: Union[PlainQuantity, None] = Field(
        default=None,
        description="MP4 energy of the molecule, unit is `hartree/particle`",
    )
    ccsd_energy: Union[PlainQuantity, None] = Field(
        default=None,
        description="CCSD energy of the molecule, unit is `hartree/particle`",
    )

    @property
    def energy(self):
        return {
            energy_type: getattr(self, f"{energy_type}_energy")
            for energy_type in (
                "ccsd",
                "mp4",
                "mp3",
                "mp2",
                "scf",
                "electronic",
            )  # precision order
            if getattr(self, f"{energy_type}_energy") is not None
        }

    @computed_field(description="Total energy, unit is `hartree/particle`")
    @property
    def total_energy(self) -> PlainQuantity:
        keys = list(self.energy.keys())
        if len(keys) > 0:
            return self.energy[keys[0]]
        else:
            return None


class ThermalEnergies(BaseDataClassWithUnit):
    """
    Thermal energy data

    ref: https://www.cup.uni-muenchen.de/ch/compchem/vib/thermo1.html
    Theorietically, the energy follow the relationship below:
    U_0 = E_tot + ZPVE
    U_T(298) = U_0 + TCE
    H_T(298) = U_0 + TCH
    G_T(298) = U_0 + TCG
    G_T(298) = H_T(298) - T * S(298)
    """

    ZPVE: Union[PlainQuantity, None] = Field(
        default=None,
        description="Zero-point vibrational energy, unit is `kcal/mol`",
    )
    U_0: Union[PlainQuantity, None] = Field(
        default=None,
        description="Zero-point energy, unit is `kcal/mol`",
    )
    TCE: Union[PlainQuantity, None] = Field(
        default=None,
        description="thermal correction to the internal energy at 298.15K, unit is `kcal/mol`",
    )
    TCH: Union[PlainQuantity, None] = Field(
        default=None,
        description="thermal correction to the enthalpy at 298.15K, unit is `kcal/mol`",
    )
    TCG: Union[PlainQuantity, None] = Field(
        default=None,
        description="thermal correction to the Gibbs free energy at 298.15K, unit is `kcal/mol`",
    )
    U_T: Union[PlainQuantity, None] = Field(
        default=None,
        description="thermal energy at 298.15K, unit is `kcal/mol`",
    )
    H_T: Union[PlainQuantity, None] = Field(
        default=None,
        description="enthalpy at 298.15K, unit is `kcal/mol`",
    )
    G_T: Union[PlainQuantity, None] = Field(
        default=None,
        description="Gibbs Free Energy at 298.15K, unit is `kcal/mol`",
    )
    S: Union[PlainQuantity, None] = Field(
        default=None,
        description="entropy at 298.15K, unit is `cal/mol/K`",
    )
    C_V: Union[PlainQuantity, None] = Field(
        default=None,
        description="heat capacity at constant volume, unit is `cal/mol/K`",
    )


class MoleculeOrbital(BaseDataClassWithUnit):
    alpha_energy: PlainQuantity = Field(
        default=None,
        description="alpha orbital energy, unit is `hartree/particle`",
    )
    beta_energy: PlainQuantity = Field(
        default=None,
        description="beta orbital energy, unit is `hartree/particle`",
    )
    alpha_occupancy: bool = Field(
        default=None,
        description="alpha orbital occupancy",
    )
    beta_occupancy: bool = Field(
        default=None,
        description="beta orbital occupancy",
    )


class MolecularOrbitals(BaseDataClassWithUnit):
    __index: int = PrivateAttr(default=0)
    # orbital energies
    alpha_energies: PlainQuantity = Field(
        default=np.array([]) * atom_ureg.hartree / atom_ureg.particle,
        description="alpha orbital energies, unit is `hartree/particle`",
    )
    beta_energies: PlainQuantity = Field(
        default=np.array([]) * atom_ureg.hartree / atom_ureg.particle,
        description="beta orbital energies, unit is `hartree/particle`",
    )
    alpha_occupancies: List[bool] = Field(
        default=[],
        description="alpha orbital occupancies",
    )
    beta_occupancies: List[bool] = Field(
        default=[],
        description="beta orbital occupancies",
    )

    @computed_field(
        description="HOMO orbital idx",
    )
    @property
    def HOMO_id(self) -> int:
        for i, alpha_occ in enumerate(self.alpha_occupancies):
            if alpha_occ == 0:
                return i - 1
        return 0

    @computed_field(
        description="LUMO orbital idx",
    )
    @property
    def LUMO_id(self) -> int:
        return max(0, self.HOMO_id + 1, len(self.alpha_occupancies) - 1)

    @computed_field(
        description="SOMO orbital idx",
    )
    @property
    def SOMO_ids(self) -> List[int]:
        return [
            i
            for i, (alpha_occ, beta_occ) in enumerate(
                zip(self.alpha_occupancies, self.beta_occupancies)
            )
            if alpha_occ and not beta_occ
        ]

    @computed_field(
        description="NHOMO orbital idx",
    )
    @property
    def NHOMO_id(self) -> int:
        return max(0, self.HOMO_id - 1)

    @computed_field(
        description="SLUMO orbital idx",
    )
    @property
    def SLUMO_id(self) -> int:
        return min(len(self.alpha_occupancies) - 1, self.LUMO_id + 1)

    @computed_field(
        description="HOMO energy",
    )
    @property
    def HOMO_energy(self) -> Union[PlainQuantity, None]:
        try:
            return self.alpha_energies[self.HOMO_id]
        except:
            return None

    @computed_field(
        description="LUMO energy",
    )
    @property
    def LUMO_energy(self) -> Union[PlainQuantity, None]:
        try:
            return self.alpha_energies[self.LUMO_id]
        except:
            return None

    @computed_field(
        description="NHOMO energy",
    )
    @property
    def NHOMO_energy(self) -> Union[PlainQuantity, None]:
        try:
            return self.alpha_energies[self.NHOMO_id]
        except:
            return None

    @computed_field(
        description="SLUMO energy",
    )
    @property
    def SLUMO_energy(self) -> Union[PlainQuantity, None]:
        try:
            return self.alpha_energies[self.SLUMO_id]
        except:
            return None

    @computed_field(
        description="HOMO-LUMO gap",
    )
    @property
    def HOMO_LUMO_gap(self) -> Union[PlainQuantity, None]:
        try:
            return self.LUMO_energy - self.HOMO_energy
        except:
            return None

    def __iter__(self):
        self.__index = 0
        return self

    def __next__(
        self,
    ) -> MoleculeOrbital:
        if self.__index >= len(self):
            raise StopIteration
        else:
            self.__index += 1
            return self[self.__index - 1]

    @overload
    def __getitem__(self, orbitalID: int) -> MoleculeOrbital: ...

    @overload
    def __getitem__(self, orbitalID: slice) -> List[MoleculeOrbital]: ...

    @overload
    def __getitem__(self, orbitalID: Union[Sequence]) -> List[MoleculeOrbital]: ...

    def __getitem__(
        self, orbitalID: Union[int, slice]
    ) -> Union[MoleculeOrbital, List[MoleculeOrbital]]:
        if isinstance(orbitalID, int):
            return MoleculeOrbital.model_validate(
                {
                    "alpha_energy": self.alpha_energies[orbitalID],
                    "beta_energy": self.beta_energies[orbitalID],
                    "alpha_occupancy": self.alpha_occupancies[orbitalID],
                    "beta_occupancy": self.beta_occupancies[orbitalID],
                }
            )
        if isinstance(orbitalID, slice):
            return [
                MoleculeOrbital.model_validate(
                    {
                        "alpha_energy": self.alpha_energies[orbital_id],
                        "beta_energy": self.beta_energies[orbital_id],
                        "alpha_occupancy": self.alpha_occupancies[orbital_id],
                        "beta_occupancy": self.beta_occupancies[orbital_id],
                    },
                )
                for orbital_id in range(*orbitalID.indices(len(self.alpha_energies)))
            ]
        else:
            return [
                MoleculeOrbital.model_validate(
                    {
                        "alpha_energy": self.alpha_energies[orbital_id],
                        "beta_energy": self.beta_energies[orbital_id],
                        "alpha_occupancy": self.alpha_occupancies[orbital_id],
                        "beta_occupancy": self.beta_occupancies[orbital_id],
                    },
                )
                for orbital_id in orbitalID
            ]

    def __len__(self) -> int:
        return len(self.alpha_energies)


class Vibration(BaseDataClassWithUnit):
    frequency: Union[PlainQuantity, None] = Field(
        default=None,
        description="Frequency of each mode, unit is `cm^-1`",
    )
    reduced_mass: Union[PlainQuantity, None] = Field(
        default=None,
        description="Reduced mass of each mode, unit is `amu`",
    )
    force_constant: Union[PlainQuantity, None] = Field(
        default=None,
        description="Force constant of each mode, unit is `mdyne/angstrom`",
    )
    IR_intensity: Union[PlainQuantity, None] = Field(
        default=None,
        description="IR intensity of each mode, unit is `kmol/mol`",
    )
    vibration_mode: Union[PlainQuantity, None] = Field(
        default=np.array([[]]) * atom_ureg.angstrom,
        description="Vibration mode of each mode, unit is `angstrom`",
    )

    @computed_field
    @property
    def is_imaginary(self) -> bool:
        return isinstance(self.frequency, PlainQuantity) and self.frequency < 0


class Vibrations(BaseDataClassWithUnit):
    __index: int = PrivateAttr(default=0)
    frequencies: PlainQuantity = Field(
        default=np.array([]) * atom_ureg.cm_1,
        description="Frequency of each mode, unit is `cm^-1`",
    )
    reduced_masses: PlainQuantity = Field(
        default=np.array([]) * atom_ureg.amu,
        description="Reduced mass of each mode, unit is `amu`",
    )
    force_constants: PlainQuantity = Field(
        default=np.array([]) * atom_ureg.mdyne / atom_ureg.angstrom,
        description="Force constant of each mode, unit is `mdyne/angstrom`",
    )
    IR_intensities: PlainQuantity = Field(
        default=np.array([]) * atom_ureg.kmol / atom_ureg.mol,
        description="IR intensity of each mode, unit is `kmol/mol`",
    )
    vibration_modes: List[Union[PlainQuantity, None]] = Field(
        default=[],
        description="Vibration mode of each mode, unit is `angstrom`",
    )

    def __iter__(self):
        self.__index = 0
        return self

    def __next__(
        self,
    ) -> Vibration:
        if self.__index >= len(self):
            raise StopIteration
        else:
            self.__index += 1
            return self[self.__index - 1]

    def __getitem__(self, frameID: int) -> Vibration:
        item_dict = {}
        try:
            item_dict["frequency"] = self.frequencies[frameID]
        except:
            pass
        try:
            item_dict["reduced_mass"] = self.reduced_masses[frameID]
        except:
            pass
        try:
            item_dict["force_constant"] = self.force_constants[frameID]
        except:
            pass
        try:
            item_dict["IR_intensity"] = self.IR_intensities[frameID]
        except:
            pass
        try:
            item_dict["vibration_mode"] = self.vibration_modes[frameID]
        except:
            pass
        return Vibration.model_validate(item_dict)

    def __len__(self) -> int:
        return len(self.frequencies)

    @property
    def imaginary_idxs(self) -> List[int]:
        return [i for i, freq in enumerate(self) if freq.is_imaginary]

    @property
    def imaginary_vibrations(self) -> "Vibrations":
        imaginary_idxs = self.imaginary_idxs
        return self.model_validate(
            {
                "frequencies": (
                    self.frequencies[imaginary_idxs]
                    if isinstance(self.frequencies, PlainQuantity)
                    else None
                ),
                "reduced_masses": (
                    self.reduced_masses[imaginary_idxs]
                    if isinstance(self.reduced_masses, PlainQuantity)
                    else None
                ),
                "force_constants": (
                    self.force_constants[imaginary_idxs]
                    if isinstance(self.force_constants, PlainQuantity)
                    else None
                ),
                "IR_intensities": (
                    self.IR_intensities[imaginary_idxs]
                    if isinstance(self.IR_intensities, PlainQuantity)
                    else None
                ),
                "vibration_modes": (
                    [self.vibration_modes[i] for i in imaginary_idxs]
                    if len(self.vibration_modes)
                    else []
                ),
            }
        )


class ChargeSpinPopulations(BaseDataClassWithUnit):
    # charge and spin populations
    mulliken_charges: List[float] = Field(
        default=[],
        description="Mulliken charges",
    )
    mulliken_spins: List[float] = Field(
        default=[],
        description="Mulliken spin densities",
    )
    apt_charges: List[float] = Field(
        default=[],
        description="Atomic polarizability tensor charges",
    )
    lowdin_charges: List[float] = Field(
        default=[],
        description="Lowdin charges",
    )
    hirshfeld_charges: List[float] = Field(
        default=[],
        description="Hirshfeld charges",
    )
    hirshfeld_spins: List[float] = Field(
        default=[],
        description="Hirshfeld spins",
    )
    hirshfeld_q_cm5: List[float] = Field(
        default=[],
        description="Hirshfeld charges in cm5",
    )
    npa_charges: List[float] = Field(
        default=[],
        description="NPA charges",
    )


class TotalSpin(BaseDataClassWithUnit):
    spin_square: Union[float, None] = Field(
        default=None,
        description="Spin square of the molecule",
    )
    spin_quantum_number: Union[float, None] = Field(
        default=None,
        description="Spin quantum number of the molecule",
    )


class Polarizability(BaseDataClassWithUnit):
    # polarizability
    electronic_spatial_extent: Union[PlainQuantity, None] = Field(
        default=None,
        description="Electronic spatial extent, unit is bohr^2",
    )
    isotropic_polarizability: Union[PlainQuantity, None] = Field(
        default=None,
        description="Isotropic polarizability, unit is bohr^3",
    )
    polarizability: List[float] = Field(
        default=[],
        description="Polarizability tensor",
    )
    dipole: Union[PlainQuantity, None] = Field(
        default=np.array([[]]) * atom_ureg.debye,
        description="Dipole moment, unit is `debye`",
    )
    quadrupole: Union[PlainQuantity, None] = Field(
        default=np.array([[]]) * atom_ureg.debye * atom_ureg.angstrom,
        description="Quadrupole moment, unit is `debye*angstrom`",
    )
    octapole: Union[PlainQuantity, None] = Field(
        default=np.array([[]]) * atom_ureg.debye * atom_ureg.angstrom**2,
        description="Octapole moment, unit is `debye*angstrom**2`",
    )
    hexadecapole: Union[PlainQuantity, None] = Field(
        default=np.array([[]]) * atom_ureg.debye * atom_ureg.angstrom**3,
        description="Hexadecapole moment, unit is `debye*angstrom**3`",
    )


class BondOrders(BaseDataClassWithUnit):
    wiberg_bond_order: np.ndarray = Field(
        default=np.array([[]]),
        description="Wiberg bond order",
    )
    mo_bond_order: np.ndarray = Field(
        default=np.array([[]]),
        description="MO bond order",
    )
    atom_atom_overlap_bond_order: np.ndarray = Field(
        default=np.array([[]]),
        description="Atom-atom overlap bond order",
    )
    nbo_bond_order: np.ndarray = Field(
        default=np.array([[]]),
        description="NBO bond order",
    )


class Dispersions(BaseDataClassWithUnit):
    C6AA: Union[PlainQuantity, None] = Field(
        default=None,
        description="Mol. C6AA dispersion, unit is `bohr^6`",
    )
    C8AA: Union[PlainQuantity, None] = Field(
        default=None,
        description="Mol. C8AA dispersion, unit is `bohr^8`",
    )


class SinglePointProperties(BaseDataClassWithUnit):
    """
    Single point properties.
    """

    vip: Union[PlainQuantity, None] = Field(
        default=None,
        description="Vertical ionization potential, unit is `eV/particle`",
    )
    vea: Union[PlainQuantity, None] = Field(
        default=None,
        description="Vertical electron affinity, unit is `eV/particle`",
    )
    gei: Union[PlainQuantity, None] = Field(
        default=None,
        description="Global Electrophilicity Index, unit is `eV/particle`",
    )
    fukui_positive: List[float] = Field(
        default=[],
        description="Fukui Index f(+)",
    )
    fukui_negative: List[float] = Field(
        default=[],
        description="Fukui Index f(-)",
    )
    fukui_zero: List[float] = Field(
        default=[],
        description="Fukui Index f(0)",
    )
    fod: List[float] = Field(
        default=[],
        description="fractional occupation density population",
    )
