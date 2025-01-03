"""
Author: TMJ
Date: 2024-06-18 20:42:42
LastEditors: TMJ
LastEditTime: 2024-06-18 20:43:03
Description: 请填写简介
"""

from abc import abstractmethod
from typing import Any, Dict, List, Sequence, Union, overload

import numpy as np
import pandas as pd
from pint.facets.plain import PlainQuantity
from pint.facets.plain.quantity import Magnitude
from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    PrivateAttr,
    computed_field,
    model_validator,
)
from typing_extensions import Self

from molop.config import molopconfig
from molop.logger.logger import moloplogger
from molop.unit import atom_ureg, unit_transform


class BaseDataClassWithUnit(BaseModel):
    model_config = ConfigDict(arbitrary_types_allowed=True)

    @overload
    def __unitless_dump__(item: PlainQuantity[np.ndarray]) -> list: ...

    @overload
    def __unitless_dump__(item: PlainQuantity[Magnitude]) -> float: ...

    @overload
    def __unitless_dump__(item: np.ndarray) -> list: ...

    @overload
    def __unitless_dump__(item: float) -> float: ...

    @overload
    def __unitless_dump__(item: int) -> int: ...

    @overload
    def __unitless_dump__(item: str) -> str: ...

    @overload
    def __unitless_dump__(item: list) -> list: ...

    @overload
    def __unitless_dump__(item: tuple) -> tuple: ...

    @overload
    def __unitless_dump__(item: dict) -> dict: ...

    @overload
    def __unitless_dump__(item: "BaseDataClassWithUnit") -> "BaseDataClassWithUnit": ...

    @staticmethod
    def __unitless_dump__(
        item: Union[
            PlainQuantity[np.ndarray],
            PlainQuantity[Magnitude],
            np.ndarray,
            float,
            int,
            str,
            list,
            tuple,
            dict,
            "BaseDataClassWithUnit",
        ]
    ) -> Union[float, int, str, list, tuple, dict]:
        if isinstance(item, PlainQuantity):
            if isinstance(item.m, np.ndarray):
                return item.m.tolist()
            else:
                return item.m
        if isinstance(item, np.ndarray):
            return item.tolist()
        if isinstance(item, list):
            return [BaseDataClassWithUnit.__unitless_dump__(i) for i in item]
        if isinstance(item, tuple):
            return tuple(BaseDataClassWithUnit.__unitless_dump__(i) for i in item)
        if isinstance(item, dict):
            return {
                k: BaseDataClassWithUnit.__unitless_dump__(v) for k, v in item.items()
            }
        if isinstance(item, BaseDataClassWithUnit):
            return item.to_unitless_dump()
        return item

    def to_unitless_dump(
        self, **kwargs
    ) -> Dict[str, Union[float, int, str, list, tuple, dict]]:
        return {
            k: BaseDataClassWithUnit.__unitless_dump__(getattr(self, k))
            for k, v in self.model_dump(**kwargs).items()
        }

    @model_validator(mode="after")
    def __parse_frame__(self) -> Self:
        if molopconfig.force_unit_transform:
            self._set_default_units()
            moloplogger.debug(
                f"Data class {self.__class__.__name__} parsed.\n" f"{self}"
            )
        return self

    @abstractmethod
    def _set_default_units(self):
        raise NotImplementedError

    def _transform_units(self, unit_dict: Dict[str, str]):
        for key, unit in unit_dict.items():
            if hasattr(self, key):
                setattr(self, key, unit_transform(getattr(self, key), unit))


class WaveFunction(BaseDataClassWithUnit): ...


# TODO


class Energies(BaseDataClassWithUnit):
    # energies
    electronic_energy: Union[PlainQuantity, None] = Field(
        default=None,
        description="Electronic energy of the molecule, unit is `hartree`",
    )
    scf_energy: Union[PlainQuantity, None] = Field(
        default=None,
        description="SCF energy of the molecule, unit is `hartree`",
    )
    mp2_energy: Union[PlainQuantity, None] = Field(
        default=None,
        description="MP2 energy of the molecule, unit is `hartree`",
    )
    mp3_energy: Union[PlainQuantity, None] = Field(
        default=None,
        description="MP3 energy of the molecule, unit is `hartree`",
    )
    mp4_energy: Union[PlainQuantity, None] = Field(
        default=None,
        description="MP4 energy of the molecule, unit is `hartree`",
    )
    ccsd_energy: Union[PlainQuantity, None] = Field(
        default=None,
        description="CCSD energy of the molecule, unit is `hartree`",
    )

    @property
    def energy(self) -> Dict[str, PlainQuantity]:
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

    @computed_field(description="Total energy, unit is `hartree`")
    @property
    def total_energy(self) -> Union[PlainQuantity, None]:
        keys = list(self.energy.keys())
        if len(keys) > 0:
            return self.energy[keys[0]].to(atom_ureg.hartree)
        else:
            return None

    def _set_default_units(self):
        self._transform_units(
            {
                "electronic_energy": atom_ureg.hartree,
                "scf_energy": atom_ureg.hartree,
                "mp2_energy": atom_ureg.hartree,
                "mp3_energy": atom_ureg.hartree,
                "mp4_energy": atom_ureg.hartree,
                "ccsd_energy": atom_ureg.hartree,
            }
        )


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
        description="Zero-point vibrational energy, unit is `hartree/particle`",
    )
    U_0: Union[PlainQuantity, None] = Field(
        default=None,
        description="Zero-point energy, unit is `hartree/particle`",
    )
    TCE: Union[PlainQuantity, None] = Field(
        default=None,
        description="thermal correction to the internal energy at 298.15K, unit is `hartree/particle`",
    )
    TCH: Union[PlainQuantity, None] = Field(
        default=None,
        description="thermal correction to the enthalpy at 298.15K, unit is `hartree/particle`",
    )
    TCG: Union[PlainQuantity, None] = Field(
        default=None,
        description="thermal correction to the Gibbs free energy at 298.15K, unit is `hartree/particle`",
    )
    U_T: Union[PlainQuantity, None] = Field(
        default=None,
        description="thermal energy at 298.15K, unit is `hartree/particle`",
    )
    H_T: Union[PlainQuantity, None] = Field(
        default=None,
        description="enthalpy at 298.15K, unit is `hartree/particle`",
    )
    G_T: Union[PlainQuantity, None] = Field(
        default=None,
        description="Gibbs Free Energy at 298.15K, unit is `hartree/particle`",
    )
    S: Union[PlainQuantity, None] = Field(
        default=None,
        description="entropy at 298.15K, unit is `cal/mol/K`",
    )
    C_V: Union[PlainQuantity, None] = Field(
        default=None,
        description="heat capacity at constant volume, unit is `cal/mol/K`",
    )

    def _set_default_units(self):
        self._transform_units(
            {
                "ZPVE": atom_ureg.hartree / atom_ureg.particle,
                "U_0": atom_ureg.hartree / atom_ureg.particle,
                "TCE": atom_ureg.hartree / atom_ureg.particle,
                "TCH": atom_ureg.hartree / atom_ureg.particle,
                "TCG": atom_ureg.hartree / atom_ureg.particle,
                "U_T": atom_ureg.hartree / atom_ureg.particle,
                "H_T": atom_ureg.hartree / atom_ureg.particle,
                "G_T": atom_ureg.hartree / atom_ureg.particle,
                "S": atom_ureg.calorie / atom_ureg.mol / atom_ureg.kelvin,
                "C_V": atom_ureg.calorie / atom_ureg.mol / atom_ureg.kelvin,
            }
        )


class MoleculeOrbital(BaseDataClassWithUnit):
    alpha_energy: PlainQuantity = Field(
        default=None,
        description="alpha orbital energy, unit is `hartree`",
    )
    beta_energy: PlainQuantity = Field(
        default=None,
        description="beta orbital energy, unit is `hartree`",
    )
    alpha_occupancy: bool = Field(
        default=None,
        description="alpha orbital occupancy",
    )
    beta_occupancy: bool = Field(
        default=None,
        description="beta orbital occupancy",
    )

    def _set_default_units(self):
        self._transform_units(
            {
                "alpha_energy": atom_ureg.hartree,
                "beta_energy": atom_ureg.hartree,
            }
        )


class MolecularOrbitals(BaseDataClassWithUnit):
    __index: int = PrivateAttr(default=0)
    # orbital energies
    alpha_energies: PlainQuantity = Field(
        default=np.array([]) * atom_ureg.hartree,
        description="alpha orbital energies, unit is `hartree`",
    )
    beta_energies: PlainQuantity = Field(
        default=np.array([]) * atom_ureg.hartree,
        description="beta orbital energies, unit is `hartree`",
    )
    alpha_occupancies: List[bool] = Field(
        default=[],
        description="alpha orbital occupancies",
    )
    beta_occupancies: List[bool] = Field(
        default=[],
        description="beta orbital occupancies",
    )

    def _set_default_units(self):
        self._transform_units(
            {
                "alpha_energies": atom_ureg.hartree,
                "beta_energies": atom_ureg.hartree,
            }
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
        return max(
            0,
            self.HOMO_id + 1,
            min(len(self.alpha_occupancies) - 1, sum(self.alpha_occupancies)),
        )

    @computed_field(
        description="beta HOMO orbital idx",
    )
    @property
    def beta_HOMO_id(self) -> int:
        for i, beta_occ in enumerate(self.beta_occupancies):
            if beta_occ == 0:
                return i - 1
        return 0

    @computed_field(
        description="beta LUMO orbital idx",
    )
    @property
    def beta_LUMO_id(self) -> int:
        return max(
            0,
            self.beta_HOMO_id + 1,
            min(len(self.beta_occupancies) - 1, sum(self.beta_occupancies)),
        )

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

    @computed_field(
        description="beta HOMO energy",
    )
    @property
    def beta_HOMO_energy(self) -> Union[PlainQuantity, None]:
        try:
            return self.beta_energies[self.beta_HOMO_id]
        except:
            return None

    @computed_field(
        description="beta LUMO energy",
    )
    @property
    def beta_LUMO_energy(self) -> Union[PlainQuantity, None]:
        try:
            return self.beta_energies[self.beta_LUMO_id]
        except:
            return None

    @computed_field(
        description="beta HOMO-LUMO gap",
    )
    @property
    def beta_HOMO_LUMO_gap(self) -> Union[PlainQuantity, None]:
        try:
            return self.beta_LUMO_energy - self.beta_HOMO_energy
        except:
            return None

    @computed_field(
        description="beta NHOMO orbital idx",
    )
    @property
    def beta_NHOMO_id(self) -> int:
        return max(0, self.beta_HOMO_id - 1)

    @computed_field(
        description="beta NHOMO energy",
    )
    @property
    def beta_NHOMO_energy(self) -> Union[PlainQuantity, None]:
        try:
            return self.beta_energies[self.beta_NHOMO_id]
        except:
            return None

    @computed_field(
        description="beta SLUMO orbital idx",
    )
    @property
    def beta_SLUMO_id(self) -> int:
        return min(len(self.beta_occupancies) - 1, self.beta_LUMO_id + 1)

    @computed_field(
        description="beta SLUMO energy",
    )
    @property
    def beta_SLUMO_energy(self) -> Union[PlainQuantity, None]:
        try:
            return self.beta_energies[self.beta_SLUMO_id]
        except:
            return None

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
    def __getitem__(self, orbitalID: Sequence) -> List[MoleculeOrbital]: ...

    def __getitem__(
        self, orbitalID: Union[int, slice, Sequence]
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

    @property
    def HOMO(self) -> MoleculeOrbital:
        if self.HOMO_id is None:
            return None
        return self[self.HOMO_id]

    @property
    def LUMO(self) -> MoleculeOrbital:
        if self.LUMO_id is None:
            return None
        return self[self.LUMO_id]

    @property
    def beta_HOMO(self) -> MoleculeOrbital:
        if self.beta_HOMO_id is None:
            return None
        return self[self.beta_HOMO_id]

    @property
    def beta_LUMO(self) -> MoleculeOrbital:
        if self.beta_LUMO_id is None:
            return None
        return self[self.beta_LUMO_id]
    
    @property
    def SOMOs(self) -> List[MoleculeOrbital]:
        return self[self.SOMO_ids]


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
        description="IR intensity of each mode, unit is `km/mol`",
    )
    vibration_mode: Union[PlainQuantity, None] = Field(
        default=np.array([[]]) * atom_ureg.angstrom,
        description="Vibration mode of each mode, unit is `angstrom`",
    )

    def _set_default_units(self):
        self._transform_units(
            {
                "frequency": atom_ureg.cm_1,
                "reduced_mass": atom_ureg.amu,
                "force_constant": atom_ureg.mdyne / atom_ureg.angstrom,
                "IR_intensity": atom_ureg.km / atom_ureg.mol,
                "vibration_mode": atom_ureg.angstrom,
            }
        )

    @computed_field
    @property
    def is_imaginary(self) -> bool:
        if self.frequency is not None and self.frequency < 0:
            return True
        else:
            return False


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
        default=np.array([]) * atom_ureg.km / atom_ureg.mol,
        description="IR intensity of each mode, unit is `km/mol`",
    )
    vibration_modes: List[Union[PlainQuantity, None]] = Field(
        default=[],
        description="Vibration mode of each mode, unit is `angstrom`",
    )

    def _set_default_units(self):
        self._transform_units(
            {
                "frequencies": atom_ureg.cm_1,
                "reduced_masses": atom_ureg.amu,
                "force_constants": atom_ureg.mdyne / atom_ureg.angstrom,
                "IR_intensities": atom_ureg.km / atom_ureg.mol,
            }
        )
        for i, vibration_mode in enumerate(self.vibration_modes):
            if vibration_mode is not None:
                self.vibration_modes[i] = unit_transform(
                    vibration_mode, atom_ureg.angstrom
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

    def _set_default_units(self): ...


class TotalSpin(BaseDataClassWithUnit):
    spin_square: Union[float, None] = Field(
        default=None,
        description="Spin square of the molecule",
    )
    spin_quantum_number: Union[float, None] = Field(
        default=None,
        description="Spin quantum number of the molecule",
    )

    def _set_default_units(self): ...


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
    anisotropic_polarizability: Union[PlainQuantity, None] = Field(
        default=None,
        description="Anisotropic polarizability, unit is bohr^3",
    )
    polarizability_tensor: Union[PlainQuantity, None] = Field(
        default=np.array([]) * atom_ureg.bohr**3,
        description="Polarizability tensor",
    )
    electric_dipole_moment: Union[PlainQuantity, None] = Field(
        default=np.array([]) * atom_ureg.debye,
        description="Electric dipole moment, unit is `debye`",
    )
    dipole: Union[PlainQuantity, None] = Field(
        default=np.array([]) * atom_ureg.debye,
        description="Dipole moment, unit is `debye`",
    )
    quadrupole: Union[PlainQuantity, None] = Field(
        default=np.array([]) * atom_ureg.debye * atom_ureg.angstrom,
        description="Quadrupole moment, unit is `debye*angstrom`",
    )
    octapole: Union[PlainQuantity, None] = Field(
        default=np.array([]) * atom_ureg.debye * atom_ureg.angstrom**2,
        description="Octapole moment, unit is `debye*angstrom**2`",
    )
    hexadecapole: Union[PlainQuantity, None] = Field(
        default=np.array([]) * atom_ureg.debye * atom_ureg.angstrom**3,
        description="Hexadecapole moment, unit is `debye*angstrom**3`",
    )

    def _set_default_units(self):
        self._transform_units(
            {
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
        )


class BondOrders(BaseDataClassWithUnit):
    wiberg_bond_order: np.ndarray = Field(
        default=np.array([[]]),
        description="Wiberg bond order",
    )
    mo_bond_order: np.ndarray = Field(
        default=np.array([[]]),
        description="MO bond order, ∑[i∈A]∑[j∈B]P(i,j)",
    )
    mayer_bond_order: np.ndarray = Field(
        default=np.array([[]]),
        description="MAYER POPULATION ANALYSIS bond order, ∑[i∈A]∑[j∈B]P(i,j)",
    )
    atom_atom_overlap_bond_order: np.ndarray = Field(
        default=np.array([[]]),
        description="Atom-atom overlap bond order",
    )
    nbo_bond_order: np.ndarray = Field(
        default=np.array([[]]),
        description="NBO bond order",
    )

    def _set_default_units(self): ...


class Dispersions(BaseDataClassWithUnit):
    C6AA: Union[PlainQuantity, None] = Field(
        default=None,
        description="Mol. C6AA dispersion, unit is `bohr^6`",
    )
    C8AA: Union[PlainQuantity, None] = Field(
        default=None,
        description="Mol. C8AA dispersion, unit is `bohr^8`",
    )

    def _set_default_units(self):
        self._transform_units(
            {
                "C6AA": atom_ureg.bohr**6,
                "C8AA": atom_ureg.bohr**8,
            }
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

    def _set_default_units(self):
        self._transform_units(
            {
                "vip": atom_ureg.eV / atom_ureg.particle,
                "vea": atom_ureg.eV / atom_ureg.particle,
                "gei": atom_ureg.eV / atom_ureg.particle,
            }
        )


class GeometryOptimizationStatus(BaseDataClassWithUnit):
    """
    Geometry optimization status.
    """

    geometry_optimized: Union[bool, None] = Field(
        default=None, description="Whether the geometry has been optimized"
    )
    energy_change_threshold: float = Field(
        default=None, description="Energy change threshold"
    )
    rms_force_threshold: float = Field(
        default=None,
        description="RMS force threshold in internal some programs use gradient, which has the same absolute value",
    )
    max_force_threshold: float = Field(
        default=None,
        description="Maximum force threshold in internal some programs use gradient, which has the same absolute value",
    )
    rms_displacement_threshold: float = Field(
        default=None, description="RMS displacement threshold in internal"
    )
    max_displacement_threshold: float = Field(
        default=None, description="Maximum displacement threshold in internal"
    )
    energy_change: float = Field(default=float("inf"), description="Energy change")
    rms_force: float = Field(
        default=float("inf"),
        description="RMS force some programs use gradient, which has the same absolute value",
    )
    max_force: float = Field(
        default=float("inf"),
        description="Maximum force some programs use gradient, which has the same absolute value",
    )
    rms_displacement: float = Field(
        default=float("inf"), description="RMS displacement"
    )
    max_displacement: float = Field(
        default=float("inf"), description="Maximum displacement"
    )

    def _set_default_units(self): ...

    @computed_field()
    @property
    def energy_change_converged(self) -> Union[bool, None]:
        """
        Whether the energy change has converged.
        """
        if self.energy_change_threshold is None:
            return None
        return self.energy_change < self.energy_change_threshold

    @computed_field()
    @property
    def rms_force_converged(self) -> Union[bool, None]:
        """
        Whether the RMS force has converged.
        """
        if self.rms_force_threshold is None:
            return None
        return self.rms_force < self.rms_force_threshold

    @computed_field()
    @property
    def max_force_converged(self) -> Union[bool, None]:
        """
        Whether the maximum force has converged.
        """
        if self.max_force_threshold is None:
            return None
        return self.max_force < self.max_force_threshold

    @computed_field()
    @property
    def rms_displacement_converged(self) -> Union[bool, None]:
        """
        Whether the RMS displacement has converged.
        """
        if self.rms_displacement_threshold is None:
            return None
        return self.rms_displacement < self.rms_displacement_threshold

    @computed_field()
    @property
    def max_displacement_converged(self) -> Union[bool, None]:
        """
        Whether the maximum displacement has converged.
        """
        if self.max_displacement_threshold is None:
            return None
        return self.max_displacement < self.max_displacement_threshold

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

    def __vector__(self):
        return abs(
            np.array(
                [
                    self.energy_change,
                    self.rms_force,
                    self.max_force,
                    self.rms_displacement,
                    self.max_displacement,
                ]
            )
        )

    def to_df(self) -> pd.DataFrame:
        metrics = (
            "energy_change",
            "rms_force",
            "max_force",
            "rms_displacement",
            "max_displacement",
        )
        df = pd.DataFrame(index=metrics, columns=["value", "threshold", "converged"])
        df["value"] = [getattr(self, metric) for metric in metrics]
        df["threshold"] = [getattr(self, f"{metric}_threshold") for metric in metrics]
        df["converged"] = [getattr(self, f"{metric}_converged") for metric in metrics]
        return df

    def __le__(self, other: "GeometryOptimizationStatus"):
        if not isinstance(other, GeometryOptimizationStatus):
            raise NotImplementedError
        if self.not_converged_num() > other.not_converged_num():
            return False
        else:
            return sum(self.__vector__() <= other.__vector__()) >= 4

    def __gt__(self, other: "GeometryOptimizationStatus"):
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
        if False not in status and any(status):
            self.geometry_optimized = True
        return self


class Status(BaseDataClassWithUnit):
    scf_converged: bool = Field(
        default=False, description="Whether the SCF has converged"
    )
    normal_terminated: bool = Field(
        default=False, description="Whether the calculation has terminated normally"
    )

    def _set_default_units(self): ...


# TODO: add NMR dataclass