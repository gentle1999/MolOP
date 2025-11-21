from typing import Any, Dict, List, Optional, Sequence, Union, overload

import numpy as np
import pandas as pd
from pint.facets.numpy.quantity import NumpyQuantity
from pint.facets.plain import PlainQuantity
from pydantic import (
    Field,
    PrivateAttr,
    computed_field,
    model_validator,
)
from typing_extensions import Self

from molop.unit import atom_ureg
from molop.utils.functions import invert_transform_coords, transform_coords

from .Bases import BaseDataClassWithUnit


class Energies(BaseDataClassWithUnit):
    # energies
    electronic_energy: Union[PlainQuantity, None] = Field(
        default=None,
        description="Electronic energy of the molecule, unit is `hartree`",
    )
    scf_energy: Union[PlainQuantity, None] = Field(
        default=None, description="SCF energy of the molecule, unit is `hartree`"
    )
    mp2_energy: Union[PlainQuantity, None] = Field(
        default=None, description="MP2 energy of the molecule, unit is `hartree`"
    )
    mp3_energy: Union[PlainQuantity, None] = Field(
        default=None, description="MP3 energy of the molecule, unit is `hartree`"
    )
    mp4_energy: Union[PlainQuantity, None] = Field(
        default=None, description="MP4 energy of the molecule, unit is `hartree`"
    )
    ccsd_energy: Union[PlainQuantity, None] = Field(
        default=None, description="CCSD energy of the molecule, unit is `hartree`"
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

    def _add_default_units(self) -> None:
        self._default_units.update(
            {
                "electronic_energy": atom_ureg.hartree,
                "scf_energy": atom_ureg.hartree,
                "mp2_energy": atom_ureg.hartree,
                "mp3_energy": atom_ureg.hartree,
                "mp4_energy": atom_ureg.hartree,
                "ccsd_energy": atom_ureg.hartree,
            }
        )

    def to_summary_dict(self, **kwargs) -> Dict[str, Any]:
        return {
            f"{key} ({getattr(self, key).units})"
            if isinstance(getattr(self, key), PlainQuantity)
            else key: getattr(self, key).m
            if isinstance(getattr(self, key), PlainQuantity)
            else getattr(self, key)
            for key in self.model_dump(**kwargs).keys()
            if getattr(self, key) is not None
        }


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

    _set_default_units: bool = PrivateAttr(default=True)

    ZPVE: Optional[PlainQuantity] = Field(
        default=None,
        description="Zero-point vibrational energy, unit is `kcal/mol`",
    )
    TCE: Optional[PlainQuantity] = Field(
        default=None,
        description="thermal correction to the internal energy at ?K, unit is `kcal/mol`",
    )
    TCH: Optional[PlainQuantity] = Field(
        default=None,
        description="thermal correction to the enthalpy at ?K, unit is `kcal/mol`",
    )
    TCG: Optional[PlainQuantity] = Field(
        default=None,
        description="thermal correction to the Gibbs free energy at ?K, unit is `kcal/mol`",
    )
    U_0: Optional[PlainQuantity] = Field(
        default=None,
        description="Zero-point energy, unit is `kcal/mol`",
    )
    U_T: Optional[PlainQuantity] = Field(
        default=None,
        description="thermal energy at ?K, unit is `kcal/mol`",
    )
    H_T: Optional[PlainQuantity] = Field(
        default=None,
        description="enthalpy at ?K, unit is `kcal/mol`",
    )
    G_T: Optional[PlainQuantity] = Field(
        default=None,
        description="Gibbs Free Energy at ?K, unit is `kcal/mol`",
    )
    S: Optional[PlainQuantity] = Field(
        default=None,
        description="entropy at ?K, unit is `cal/mol/K`",
    )
    C_V: Optional[PlainQuantity] = Field(
        default=None,
        description="heat capacity at constant volume, unit is `cal/mol/K`",
    )

    def _add_default_units(self) -> None:
        self._default_units.update(
            {
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
            }
        )

    def to_summary_dict(self, **kwargs) -> Dict[str, Any]:
        return {
            f"{key} ({getattr(self, key).units})"
            if isinstance(getattr(self, key), PlainQuantity)
            else key: getattr(self, key).m
            if isinstance(getattr(self, key), PlainQuantity)
            else getattr(self, key)
            for key in self.model_dump(**kwargs).keys()
            if getattr(self, key) is not None
        }


class MoleculeOrbital(BaseDataClassWithUnit):
    alpha_energy: Optional[PlainQuantity] = Field(
        default=None, description="alpha orbital energy, unit is `hartree`"
    )
    beta_energy: Optional[PlainQuantity] = Field(
        default=None, description="beta orbital energy, unit is `hartree`"
    )
    alpha_occupancy: Optional[bool] = Field(
        default=None, description="alpha orbital occupancy"
    )
    alpha_symmetry: Optional[str] = Field(
        default=None, description="alpha orbital symmetry"
    )
    beta_occupancy: Optional[bool] = Field(
        default=None, description="beta orbital occupancy"
    )
    beta_symmetry: Optional[str] = Field(
        default=None, description="beta orbital symmetry"
    )
    coefficient: Optional[np.ndarray] = Field(
        default=None, description="coefficient of the orbital"
    )

    def _add_default_units(self) -> None:
        self._default_units.update(
            {
                "alpha_energy": atom_ureg.hartree,
                "beta_energy": atom_ureg.hartree,
            }
        )

    def to_summary_dict(self, **kwargs) -> Dict[str, Any]:
        return {
            f"{key} ({getattr(self, key).units})"
            if isinstance(getattr(self, key), PlainQuantity)
            else key: getattr(self, key).m
            if isinstance(getattr(self, key), PlainQuantity)
            else getattr(self, key)
            for key in self.model_dump(**kwargs).keys()
            if getattr(self, key) is not None
        }


class MolecularOrbitals(BaseDataClassWithUnit):
    __index: int = PrivateAttr(default=0)
    # orbital energies
    electronic_state: Optional[str] = Field(
        default=None, description="electronic state of the molecule"
    )
    alpha_energies: NumpyQuantity = Field(
        default=np.array([]) * atom_ureg.hartree,
        description="alpha orbital energies, unit is `hartree`",
    )
    beta_energies: NumpyQuantity = Field(
        default=np.array([]) * atom_ureg.hartree,
        description="beta orbital energies, unit is `hartree`",
    )
    alpha_occupancies: List[bool | None] = Field(
        default=[], description="alpha orbital occupancies"
    )
    beta_occupancies: List[bool | None] = Field(
        default=[], description="beta orbital occupancies"
    )
    alpha_symmetries: List[str | None] = Field(
        default=[], description="alpha orbital symmetries"
    )
    beta_symmetries: List[str | None] = Field(
        default=[], description="beta orbital symmetries"
    )
    coefficients: List[np.ndarray | None] = Field(
        default=[], description="coefficients of the orbitals"
    )

    def _add_default_units(self) -> None:
        self._default_units.update(
            {"alpha_energies": atom_ureg.hartree, "beta_energies": atom_ureg.hartree}
        )

    @computed_field(description="HOMO orbital idx")
    @property
    def HOMO_id(self) -> Optional[int]:
        for i, alpha_occ in enumerate(self.alpha_occupancies):
            if alpha_occ is None:
                continue
            if i > 0 and not alpha_occ and self.alpha_occupancies[i - 1]:
                return i - 1
        return None

    @computed_field(description="LUMO orbital idx")
    @property
    def LUMO_id(self) -> Optional[int]:
        if self.HOMO_id is None:
            return None
        return self.HOMO_id + 1

    @computed_field(description="beta HOMO orbital idx")
    @property
    def beta_HOMO_id(self) -> Optional[int]:
        for i, beta_occ in enumerate(self.beta_occupancies):
            if beta_occ is None:
                continue
            if i > 0 and not beta_occ and self.beta_occupancies[i - 1]:
                return i - 1
        return None

    @computed_field(description="beta LUMO orbital idx")
    @property
    def beta_LUMO_id(self) -> Optional[int]:
        if self.beta_HOMO_id is None:
            return None
        return self.beta_HOMO_id + 1

    @computed_field(description="SOMO orbital idx")
    @property
    def SOMO_ids(self) -> List[int]:
        if len(self.beta_occupancies) == 0:
            return []
        return [
            i
            for i, (alpha_occ, beta_occ) in enumerate(
                zip(self.alpha_occupancies, self.beta_occupancies, strict=True)
            )
            if alpha_occ and not beta_occ
        ]

    @computed_field(description="NHOMO orbital idx")
    @property
    def NHOMO_id(self) -> Optional[int]:
        if self.HOMO_id is None:
            return None
        if self.HOMO_id == 0:
            return None
        return self.HOMO_id - 1

    @computed_field(description="SLUMO orbital idx")
    @property
    def SLUMO_id(self) -> Optional[int]:
        if self.LUMO_id is None:
            return None
        if self.LUMO_id == len(self.alpha_occupancies):
            return None
        return self.LUMO_id + 1

    @computed_field(
        description="HOMO energy",
    )
    @property
    def HOMO_energy(self) -> Optional[PlainQuantity]:
        if self.HOMO_id is None:
            return None
        if len(self.alpha_energies) <= self.HOMO_id:
            return None
        return self.alpha_energies[self.HOMO_id]

    @computed_field(description="LUMO energy")
    @property
    def LUMO_energy(self) -> Optional[PlainQuantity]:
        if self.LUMO_id is None:
            return None
        if len(self.alpha_energies) <= self.LUMO_id:
            return None
        return self.alpha_energies[self.LUMO_id]

    @computed_field(description="NHOMO energy")
    @property
    def NHOMO_energy(self) -> Optional[PlainQuantity]:
        if self.NHOMO_id is None:
            return None
        if len(self.alpha_energies) <= self.NHOMO_id:
            return None
        return self.alpha_energies[self.NHOMO_id]

    @computed_field(description="SLUMO energy")
    @property
    def SLUMO_energy(self) -> Optional[PlainQuantity]:
        if self.SLUMO_id is None:
            return None
        if len(self.alpha_energies) <= self.SLUMO_id:
            return None
        return self.alpha_energies[self.SLUMO_id]

    @computed_field(description="HOMO-LUMO gap")
    @property
    def HOMO_LUMO_gap(self) -> Optional[PlainQuantity]:
        if self.HOMO_energy is None or self.LUMO_energy is None:
            return None
        return self.LUMO_energy - self.HOMO_energy

    def __iter__(self):
        self.__index = 0
        return self

    @computed_field(description="beta HOMO energy")
    @property
    def beta_HOMO_energy(self) -> Optional[PlainQuantity]:
        if self.beta_HOMO_id is None:
            return None
        if len(self.beta_energies) <= self.beta_HOMO_id:
            return None
        return self.beta_energies[self.beta_HOMO_id]

    @computed_field(description="beta LUMO energy")
    @property
    def beta_LUMO_energy(self) -> Optional[PlainQuantity]:
        if self.beta_LUMO_id is None:
            return None
        if len(self.beta_energies) <= self.beta_LUMO_id:
            return None
        return self.beta_energies[self.beta_LUMO_id]

    @computed_field(description="beta HOMO-LUMO gap")
    @property
    def beta_HOMO_LUMO_gap(self) -> Optional[PlainQuantity]:
        if self.beta_HOMO_energy is None or self.beta_LUMO_energy is None:
            return None
        return self.beta_LUMO_energy - self.beta_HOMO_energy

    @computed_field(description="beta NHOMO orbital idx")
    @property
    def beta_NHOMO_id(self) -> Optional[int]:
        if self.beta_HOMO_id is None:
            return None
        if self.beta_HOMO_id == 0:
            return None
        return self.beta_HOMO_id - 1

    @computed_field(description="beta NHOMO energy")
    @property
    def beta_NHOMO_energy(self) -> Optional[PlainQuantity]:
        if self.beta_NHOMO_id is None:
            return None
        if len(self.beta_energies) <= self.beta_NHOMO_id:
            return None
        return self.beta_energies[self.beta_NHOMO_id]

    @computed_field(description="beta SLUMO orbital idx")
    @property
    def beta_SLUMO_id(self) -> Optional[int]:
        if self.beta_LUMO_id is None:
            return None
        if self.beta_LUMO_id == len(self.beta_occupancies):
            return None
        return self.beta_LUMO_id + 1

    @computed_field(description="beta SLUMO energy")
    @property
    def beta_SLUMO_energy(self) -> Union[PlainQuantity, None]:
        if self.beta_SLUMO_id is None:
            return None
        if len(self.beta_energies) <= self.beta_SLUMO_id:
            return None
        return self.beta_energies[self.beta_SLUMO_id]

    def __next__(
        self,
    ) -> MoleculeOrbital:
        if self.__index >= len(self):
            raise StopIteration
        else:
            self.__index += 1
            return self[self.__index - 1]

    @overload
    def __getitem__(self, orbitalIDX: int) -> MoleculeOrbital: ...
    @overload
    def __getitem__(self, orbitalIDX: slice) -> List[MoleculeOrbital]: ...
    @overload
    def __getitem__(self, orbitalIDX: Sequence) -> List[MoleculeOrbital]: ...
    def __getitem__(
        self, orbitalIDX: Union[int, slice, Sequence]
    ) -> Union[MoleculeOrbital, List[MoleculeOrbital]]:
        def get_item(seq: Sequence[Any], idx: int):
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
                    "alpha_energy": get_item(self.alpha_energies, orbitalIDX),  # type: ignore
                    "beta_energy": get_item(self.beta_energies, orbitalIDX),  # type: ignore
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

    @property
    def HOMO(self) -> Optional[MoleculeOrbital]:
        if self.HOMO_id is None:
            return None
        return self[self.HOMO_id]

    @property
    def LUMO(self) -> Optional[MoleculeOrbital]:
        if self.LUMO_id is None:
            return None
        return self[self.LUMO_id]

    @property
    def beta_HOMO(self) -> Optional[MoleculeOrbital]:
        if self.beta_HOMO_id is None:
            return None
        return self[self.beta_HOMO_id]

    @property
    def beta_LUMO(self) -> Optional[MoleculeOrbital]:
        if self.beta_LUMO_id is None:
            return None
        return self[self.beta_LUMO_id]

    @property
    def SOMOs(self) -> List[MoleculeOrbital]:
        return self[self.SOMO_ids]

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

    def to_summary_dict(self, **kwargs) -> Dict[str, Any]:
        return {
            "electronic_state": self.electronic_state,
            "HOMO_energy": None if self.HOMO_energy is None else self.HOMO_energy.m,
            "LUMO_energy": None if self.LUMO_energy is None else self.LUMO_energy.m,
            "HOMO-LUMO_gap": None
            if self.HOMO_LUMO_gap is None
            else self.HOMO_LUMO_gap.m,
        }


class Vibration(BaseDataClassWithUnit):
    frequency: Optional[PlainQuantity] = Field(
        default=None, description="Frequency of each mode, unit is `cm^-1`"
    )
    reduced_mass: Optional[PlainQuantity] = Field(
        default=None, description="Reduced mass of each mode, unit is `amu`"
    )
    force_constant: Optional[PlainQuantity] = Field(
        default=None,
        description="Force constant of each mode, unit is `mdyne/angstrom`",
    )
    IR_intensity: Optional[PlainQuantity] = Field(
        default=None, description="IR intensity of each mode, unit is `km/mol`"
    )
    vibration_mode: NumpyQuantity = Field(
        default=np.array([[]]) * atom_ureg.angstrom,
        description="Vibration mode of each mode, unit is `angstrom`",
    )

    def _add_default_units(self) -> None:
        self._default_units.update(
            {
                "frequency": atom_ureg.cm_1,
                "reduced_mass": atom_ureg.amu,
                "force_constant": atom_ureg.Unit("mdyne/angstrom"),
                "IR_intensity": atom_ureg.Unit("km/mol"),
                "vibration_mode": atom_ureg.angstrom,
            }
        )

    @computed_field
    @property
    def is_imaginary(self) -> bool:
        return bool(self.frequency is not None and self.frequency < 0)

    def to_summary_dict(self, **kwargs) -> Dict[str, Any]:
        return {
            f"{key} ({getattr(self, key).units})"
            if isinstance(getattr(self, key), PlainQuantity)
            else key: getattr(self, key).m
            if isinstance(getattr(self, key), PlainQuantity)
            else getattr(self, key)
            for key in self.model_dump(**kwargs).keys()
            if getattr(self, key) is not None
        }

    def transform_orientation(
        self, transformation_matrix: np.ndarray, inverse: bool = False
    ) -> None:
        if inverse:
            self.vibration_mode = (
                invert_transform_coords(self.vibration_mode.m, transformation_matrix)
                * self.vibration_mode.u
            )
        else:
            self.vibration_mode = (
                transform_coords(self.vibration_mode.m, transformation_matrix)
                * self.vibration_mode.u
            )


class Vibrations(BaseDataClassWithUnit):
    __index: int = PrivateAttr(default=0)
    frequencies: NumpyQuantity = Field(
        default=np.array([]) * atom_ureg.cm_1,
        description="Frequency of each mode, unit is `cm^-1`",
    )
    reduced_masses: NumpyQuantity = Field(
        default=np.array([]) * atom_ureg.amu,
        description="Reduced mass of each mode, unit is `amu`",
    )
    force_constants: NumpyQuantity = Field(
        default=np.array([]) * atom_ureg.mdyne / atom_ureg.angstrom,
        description="Force constant of each mode, unit is `mdyne/angstrom`",
    )
    IR_intensities: NumpyQuantity = Field(
        default=np.array([]) * atom_ureg.km / atom_ureg.mol,
        description="IR intensity of each mode, unit is `km/mol`",
    )
    vibration_modes: List[NumpyQuantity] = Field(
        default=[],
        description="Vibration mode of each mode, unit is `angstrom`",
    )

    def _add_default_units(self) -> None:
        self._default_units.update(
            {
                "frequency": atom_ureg.cm_1,
                "reduced_mass": atom_ureg.amu,
                "force_constant": atom_ureg.Unit("mdyne/angstrom"),
                "IR_intensity": atom_ureg.Unit("km/mol"),
                "vibration_mode": atom_ureg.angstrom,
            }
        )

    def __iter__(self) -> Self:
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

    @overload
    def __getitem__(self, frameID: int) -> Vibration: ...
    @overload
    def __getitem__(self, frameID: slice) -> List[Vibration]: ...
    @overload
    def __getitem__(self, frameID: Sequence) -> List[Vibration]: ...
    def __getitem__(
        self, frameID: int | slice | Sequence
    ) -> Vibration | List[Vibration]:
        item_dict = {}
        if isinstance(frameID, int):
            if len(self.frequencies) > frameID:
                item_dict["frequency"] = self.frequencies[frameID]
            if len(self.reduced_masses) > frameID:
                item_dict["reduced_mass"] = self.reduced_masses[frameID]
            if len(self.force_constants) > frameID:
                item_dict["force_constant"] = self.force_constants[frameID]
            if len(self.IR_intensities) > frameID:
                item_dict["IR_intensity"] = self.IR_intensities[frameID]
            if len(self.vibration_modes) > frameID:
                item_dict["vibration_mode"] = self.vibration_modes[frameID]
            return Vibration.model_validate(item_dict)
        elif isinstance(frameID, slice):
            return [self[idx] for idx in range(*frameID.indices(len(self.frequencies)))]
        else:
            return [self[idx] for idx in frameID]

    def __len__(self) -> int:
        return len(self.frequencies)

    @property
    def num_imaginary(self) -> int:
        return len(self.imaginary_idxs)

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

    def transform_orientation(
        self, transformation_matrix: np.ndarray, inverse: bool = False
    ) -> None:
        if inverse:
            self.vibration_modes = [
                invert_transform_coords(mode.m, transformation_matrix) * mode.u
                for mode in self.vibration_modes
            ]
        else:
            self.vibration_modes = [
                transform_coords(mode.m, transformation_matrix) * mode.u
                for mode in self.vibration_modes
            ]

    def to_summary_dict(self, **kwargs) -> Dict[str, Any]:
        return {"num_imaginary": self.num_imaginary, "num_vibrations": len(self)}


class ChargeSpinPopulations(BaseDataClassWithUnit):
    # charge and spin populations
    mulliken_charges: List[float] = Field(default=[], description="Mulliken charges")
    mulliken_spins: List[float] = Field(
        default=[], description="Mulliken spin densities"
    )
    apt_charges: List[float] = Field(
        default=[], description="Atomic polarizability tensor charges"
    )
    lowdin_charges: List[float] = Field(default=[], description="Lowdin charges")
    hirshfeld_charges: List[float] = Field(default=[], description="Hirshfeld charges")
    hirshfeld_spins: List[float] = Field(default=[], description="Hirshfeld spins")
    hirshfeld_q_cm5: List[float] = Field(
        default=[], description="Hirshfeld charges in cm5"
    )
    npa_charges: List[float] = Field(default=[], description="NPA charges")

    def _add_default_units(self): ...

    @model_validator(mode="after")
    def validate_charge_spin_populations(self) -> Self:
        available_populations = [
            getattr(self, pop)
            for pop in self.model_fields_set
            if len(getattr(self, pop)) > 0
        ]
        assert all(
            len(pop) == len(available_populations[0]) for pop in available_populations
        ), "All populations must have the same length"
        return self

    def to_summary_dict(self, **kwargs) -> Dict[str, Any]:
        return {}


class TotalSpin(BaseDataClassWithUnit):
    spin_square: Union[float, None] = Field(
        default=None, description="Spin square of the molecule"
    )
    spin_quantum_number: Union[float, None] = Field(
        default=None, description="Spin quantum number of the molecule"
    )

    def _add_default_units(self): ...

    def to_summary_dict(self, **kwargs) -> Dict[str, Any]:
        return self.to_unitless_dump(**kwargs)


class Polarizability(BaseDataClassWithUnit):
    # polarizability
    electronic_spatial_extent: Optional[PlainQuantity] = Field(
        default=None, description="Electronic spatial extent, unit is bohr^2"
    )
    isotropic_polarizability: Optional[PlainQuantity] = Field(
        default=None, description="Isotropic polarizability, unit is bohr^3"
    )
    anisotropic_polarizability: Union[PlainQuantity, None] = Field(
        default=None, description="Anisotropic polarizability, unit is bohr^3"
    )
    polarizability_tensor: Union[PlainQuantity, None] = Field(
        default=np.array([]) * atom_ureg.bohr**3, description="Polarizability tensor"
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
    traceless_quadrupole: Union[PlainQuantity, None] = Field(
        default=np.array([]) * atom_ureg.debye * atom_ureg.angstrom,
        description="Traceless quadrupole moment, unit is `debye*angstrom`",
    )
    octapole: Union[PlainQuantity, None] = Field(
        default=np.array([]) * atom_ureg.debye * atom_ureg.angstrom**2,
        description="Octapole moment, unit is `debye*angstrom**2`",
    )
    hexadecapole: Union[PlainQuantity, None] = Field(
        default=np.array([]) * atom_ureg.debye * atom_ureg.angstrom**3,
        description="Hexadecapole moment, unit is `debye*angstrom**3`",
    )

    def _add_default_units(self) -> None:
        self._default_units.update(
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

    def to_summary_dict(self, **kwargs) -> Dict[str, Any]:
        return {
            f"{key} ({getattr(self, key).units})"
            if isinstance(getattr(self, key), PlainQuantity)
            else key: getattr(self, key).m
            if isinstance(getattr(self, key), PlainQuantity)
            else getattr(self, key)
            for key in self.model_dump(**kwargs).keys()
            if getattr(self, key) is not None
        }


class BondOrders(BaseDataClassWithUnit):
    wiberg_bond_order: np.ndarray = Field(
        default=np.array([[]]), description="Wiberg bond order"
    )
    mo_bond_order: np.ndarray = Field(
        default=np.array([[]]), description="MO bond order, ∑[i∈A]∑[j∈B]P(i,j)"
    )
    mayer_bond_order: np.ndarray = Field(
        default=np.array([[]]),
        description="MAYER POPULATION ANALYSIS bond order, ∑[i∈A]∑[j∈B]P(i,j)",
    )
    atom_atom_overlap_bond_order: np.ndarray = Field(
        default=np.array([[]]), description="Atom-atom overlap bond order"
    )
    nbo_bond_order: np.ndarray = Field(
        default=np.array([[]]), description="NBO bond order"
    )

    def _add_default_units(self): ...

    def to_summary_dict(self, **kwargs) -> Dict[str, Any]:
        return {}


class Dispersions(BaseDataClassWithUnit):
    C6AA: Union[PlainQuantity, None] = Field(
        default=None,
        description="Mol. C6AA dispersion, unit is `bohr^6`",
    )
    C8AA: Union[PlainQuantity, None] = Field(
        default=None,
        description="Mol. C8AA dispersion, unit is `bohr^8`",
    )

    def _add_default_units(self) -> None:
        self._default_units.update(
            {
                "C6AA": atom_ureg.bohr**6,
                "C8AA": atom_ureg.bohr**8,
            }
        )

    def to_summary_dict(self, **kwargs) -> Dict[str, Any]:
        return {
            f"{key} ({getattr(self, key).units})"
            if isinstance(getattr(self, key), PlainQuantity)
            else key: getattr(self, key).m
            if isinstance(getattr(self, key), PlainQuantity)
            else getattr(self, key)
            for key in self.model_dump(**kwargs).keys()
            if getattr(self, key) is not None
        }


class SinglePointProperties(BaseDataClassWithUnit):
    """
    Single point properties.
    """

    vip: Union[PlainQuantity, None] = Field(
        default=None, description="Vertical ionization potential, unit is `eV/particle`"
    )
    vea: Union[PlainQuantity, None] = Field(
        default=None, description="Vertical electron affinity, unit is `eV/particle`"
    )
    gei: Union[PlainQuantity, None] = Field(
        default=None, description="Global Electrophilicity Index, unit is `eV/particle`"
    )
    fukui_positive: List[float] = Field(default=[], description="Fukui Index f(+)")
    fukui_negative: List[float] = Field(default=[], description="Fukui Index f(-)")
    fukui_zero: List[float] = Field(default=[], description="Fukui Index f(0)")
    fod: List[float] = Field(
        default=[], description="fractional occupation density population"
    )

    def _add_default_units(self) -> None:
        self._default_units.update(
            {
                "vip": atom_ureg.Unit("eV / particle"),
                "vea": atom_ureg.Unit("eV / particle"),
                "gei": atom_ureg.Unit("eV / particle"),
            }
        )

    def to_summary_dict(self, **kwargs) -> Dict[str, Any]:
        return {}


class GeometryOptimizationStatus(BaseDataClassWithUnit):
    """
    Geometry optimization status.
    """

    geometry_optimized: Union[bool, None] = Field(
        default=None, description="Whether the geometry has been optimized"
    )
    energy_change_threshold: Optional[float] = Field(
        default=None, description="Energy change threshold"
    )
    rms_force_threshold: Optional[float] = Field(
        default=None,
        description="RMS force threshold in internal some programs use gradient, which has the same absolute value",
    )
    max_force_threshold: Optional[float] = Field(
        default=None,
        description="Maximum force threshold in internal some programs use gradient, which has the same absolute value",
    )
    rms_displacement_threshold: Optional[float] = Field(
        default=None, description="RMS displacement threshold in internal"
    )
    max_displacement_threshold: Optional[float] = Field(
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

    def _add_default_units(self): ...
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

    def __vector__(
        self,
    ) -> (
        np.ndarray[Any, np.dtype[np.bool_]]
        | np.ndarray[Any, np.dtype[np.floating[Any]]]
    ):
        return np.abs(
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
        # TODO: more accurate comparison
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

    def to_summary_dict(self, **kwargs) -> Dict[str, Any]:
        return {
            "geometry_optimized": self.geometry_optimized,
            "energy_change_converged": self.energy_change_converged,
            "rms_force_converged": self.rms_force_converged,
            "max_force_converged": self.max_force_converged,
            "rms_displacement_converged": self.rms_displacement_converged,
            "max_displacement_converged": self.max_displacement_converged,
        }


class Status(BaseDataClassWithUnit):
    scf_converged: bool = Field(
        default=False, description="Whether the SCF has converged"
    )
    normal_terminated: bool = Field(
        default=False, description="Whether the calculation has terminated normally"
    )

    def _add_default_units(self): ...

    def to_summary_dict(self, **kwargs) -> Dict[str, Any]:
        return {
            f"{key} ({getattr(self, key).units})"
            if isinstance(getattr(self, key), PlainQuantity)
            else key: getattr(self, key).m
            if isinstance(getattr(self, key), PlainQuantity)
            else getattr(self, key)
            for key in self.model_dump(**kwargs).keys()
            if getattr(self, key) is not None
        }


class ShieldingTensor(BaseDataClassWithUnit):
    """
    Shielding tensor data.
    """

    atom: str = Field(default="", description="Atom element")
    shielding_tensor: NumpyQuantity = Field(
        default=np.zeros((3, 3)) * atom_ureg.ppm,
        description="Shielding tensor, unit is `ppm`",
    )

    @computed_field()
    @property
    def isotropic(self) -> PlainQuantity:
        return np.mean(np.trace(self.shielding_tensor))

    @computed_field()
    @property
    def anisotropy(self) -> PlainQuantity:
        eigenvalues = np.linalg.eigvals(self.shielding_tensor)
        # 按照 Haeberlen 约定排序: |s_zz - s_iso| >= |s_xx - s_iso| >= |s_yy - s_iso|
        s_iso = np.mean(eigenvalues)
        sorted_eigs = sorted(eigenvalues, key=lambda x: abs(x - s_iso), reverse=True)
        return sorted_eigs[0] - (sorted_eigs[1] + sorted_eigs[2]) / 2

    def _add_default_units(self) -> None:
        self._default_units.update({"shielding_tensor": atom_ureg.ppm})

    def to_summary_dict(self, **kwargs) -> Dict[str, Any]:
        return {
            f"{key} ({getattr(self, key).units})"
            if isinstance(getattr(self, key), PlainQuantity)
            else key: getattr(self, key).m
            if isinstance(getattr(self, key), PlainQuantity)
            else getattr(self, key)
            for key in self.model_dump(**kwargs).keys()
            if getattr(self, key) is not None
        }


class NMR(BaseDataClassWithUnit):
    shielding_tensors: List[ShieldingTensor] = Field(
        default=[], description="NMR shielding tensors"
    )
    spin_spin_coupling_k: Optional[NumpyQuantity] = Field(
        default=None,
        description="Spin-spin coupling constant, unit is `Hz`",
    )
    spin_spin_coupling_j: Optional[NumpyQuantity] = Field(
        default=None,
        description="Spin-spin coupling constant, unit is `Hz`",
    )

    def _add_default_units(self) -> None:
        self._default_units.update(
            {
                "spin_spin_coupling_k": atom_ureg.Hz,
                "spin_spin_coupling_j": atom_ureg.Hz,
            }
        )

    def to_summary_dict(self, **kwargs) -> Dict[str, Any]:
        return self.to_unitless_dump(**kwargs)


class ImplicitSolvation(BaseDataClassWithUnit):
    solvent: Optional[str] = Field(
        default=None,
        description="Solvent used in the QM calculation",
    )
    solvent_model: Optional[str] = Field(
        default=None,
        description="Solvent model used in the QM calculation",
    )
    atomic_radii: Optional[str] = Field(
        default=None,
        description="Atomic radii used in the QM calculation",
    )
    solvent_epsilon: Optional[float] = Field(
        default=None,
        description="Solvent dielectric constant used in the QM calculation",
    )
    solvent_epsilon_infinite: Optional[float] = Field(
        default=None,
        description="Solvent epsilon infinite used in the QM calculation",
    )

    def _add_default_units(self) -> None: ...

    def to_summary_dict(self, **kwargs) -> Dict[str, Any]:
        return {
            f"{key} ({getattr(self, key).units})"
            if isinstance(getattr(self, key), PlainQuantity)
            else key: getattr(self, key).m
            if isinstance(getattr(self, key), PlainQuantity)
            else getattr(self, key)
            for key in self.model_dump(**kwargs).keys()
            if getattr(self, key) is not None
        }
