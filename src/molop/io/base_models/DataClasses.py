from collections.abc import Iterator, Sequence
from typing import Any, ClassVar, cast, overload

import numpy as np
import pandas as pd
from pint._typing import UnitLike
from pint.facets.numpy.quantity import NumpyQuantity
from pint.facets.plain import PlainQuantity
from pydantic import Field, PrivateAttr, computed_field, model_validator
from typing_extensions import Self

from molop.unit import atom_ureg
from molop.utils.functions import invert_transform_coords, transform_coords

from .Bases import BaseDataClassWithUnit


def _is_quantity(value: Any) -> bool:
    return isinstance(value, PlainQuantity)


class AtomInInternalCoords(BaseDataClassWithUnit):
    default_units: ClassVar[dict[str, UnitLike]] = {
        "distance": atom_ureg.angstrom,
        "angle": atom_ureg.degree,
        "dihedral": atom_ureg.degree,
    }
    set_default_units: ClassVar[bool] = True

    symbol: str = Field(description="Atom symbol, also support Dummy Atom")
    distance_to_index: int = Field(default=0, description="Index of atom to define the distance")
    distance: PlainQuantity = Field(
        default=0 * atom_ureg.angstrom, description="Distance between two atoms, unit is `angstrom`"
    )
    angle_to_index: int = Field(default=0, description="Index of atom to define the angle")
    angle: PlainQuantity = Field(
        default=0 * atom_ureg.degree, description="Angle between three atoms, unit is `degree`"
    )
    dihedral_to_index: int = Field(
        default=0, description="Index of atom to define the dihedral angle"
    )
    dihedral: PlainQuantity = Field(
        default=0 * atom_ureg.degree,
        description="Dihedral angle between four atoms, unit is `degree`",
    )
    zmat_format: int = Field(
        default=0,
        description="Z-matrix row format selector: 0=standard dihedral, 1=alternate two-angle row",
    )

    is_dummy: bool = Field(default=False, description="Whether the atom is a dummy atom")
    is_ghost: bool = Field(default=False, description="Whether the atom is a ghost atom")


class InternalCoords(BaseDataClassWithUnit):
    atoms: list[AtomInInternalCoords] = Field(
        default_factory=list, description="Atoms in the internal coordinates"
    )

    @classmethod
    def from_cartesian_coords(cls, symbols: Sequence[str], coords: NumpyQuantity) -> Self:
        """
        Create internal coordinates from atom symbols and cartesian coordinates.

        Parameters:
            symbols (Sequence[str]): Atom symbols.
            coords (NumpyQuantity): Cartesian coordinates of the atoms.

        Returns:
            Self: Internal coordinates.
        """
        if len(symbols) != coords.shape[0]:
            raise ValueError(
                f"Number of symbols ({len(symbols)}) does not match number of coordinates ({coords.shape[0]})"
            )

        xyz = coords.to(atom_ureg.angstrom).magnitude
        atoms: list[AtomInInternalCoords] = []

        for i, symbol in enumerate(symbols):
            atom = AtomInInternalCoords(
                symbol=symbol,
                distance_to_index=max(i - 1, 0),
                angle_to_index=max(i - 2, 0),
                dihedral_to_index=max(i - 3, 0),
            )

            if i > 0:
                atom.distance = np.linalg.norm(xyz[i] - xyz[i - 1]) * atom_ureg.angstrom

            if i > 1:
                v1 = xyz[i - 2] - xyz[i - 1]
                v2 = xyz[i] - xyz[i - 1]
                cosang = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
                cosang = np.clip(cosang, -1.0, 1.0)
                atom.angle = np.arccos(cosang) * atom_ureg.radian

            if i > 2:
                p0, p1, p2, p3 = xyz[i - 3], xyz[i - 2], xyz[i - 1], xyz[i]

                b0 = p1 - p0
                b1 = p2 - p1
                b2 = p3 - p2

                b1u = b1 / np.linalg.norm(b1)

                v = b0 - np.dot(b0, b1u) * b1u
                w = b2 - np.dot(b2, b1u) * b1u

                x = np.dot(v, w)
                y = np.dot(np.cross(b1u, v), w)
                atom.dihedral = np.arctan2(y, x) * atom_ureg.radian

            atoms.append(atom)

        return cls(atoms=atoms)

    def get_symbols(self) -> list[str]:
        return [atom.symbol for atom in self.atoms if not atom.is_dummy and not atom.is_ghost]

    def to_cartesian_coords(self) -> NumpyQuantity:
        """
        Convert internal coordinates to cartesian coordinates.

        Returns:
            NumpyQuantity: Cartesian coordinates of the atoms.
        """
        n = len(self.atoms)
        if n == 0:
            return np.zeros((0, 3)) * atom_ureg.angstrom

        coords = np.zeros((n, 3), dtype=float)

        # Atom 0
        coords[0] = np.array([0.0, 0.0, 0.0])

        if n == 1:
            return coords * atom_ureg.angstrom

        # Atom 1
        atom = self.atoms[1]
        r = atom.distance.to(atom_ureg.angstrom).magnitude
        coords[1] = np.array([r, 0.0, 0.0])

        if n == 2:
            return coords * atom_ureg.angstrom

        # Atom 2: xy plane
        atom = self.atoms[2]
        r = atom.distance.to(atom_ureg.angstrom).magnitude
        theta = atom.angle.to(atom_ureg.radian).magnitude

        j = atom.distance_to_index
        k = atom.angle_to_index

        rj = coords[j]
        rk = coords[k]

        e1 = rj - rk
        e1_norm = np.linalg.norm(e1)
        e1 = np.array([1.0, 0.0, 0.0]) if np.isclose(e1_norm, 0.0) else e1 / e1_norm

        # choose one perpendicular direction
        trial = np.array([0.0, 0.0, 1.0])
        if np.linalg.norm(np.cross(e1, trial)) < 1e-8:
            trial = np.array([0.0, 1.0, 0.0])

        e2 = np.cross(trial, e1)
        e2 = e2 / np.linalg.norm(e2)

        coords[2] = rj + r * (-np.cos(theta) * e1 + np.sin(theta) * e2)

        for i in range(3, n):
            atom = self.atoms[i]
            r = atom.distance.to(atom_ureg.angstrom).magnitude
            theta = atom.angle.to(atom_ureg.radian).magnitude

            c = atom.distance_to_index
            b = atom.angle_to_index
            a = atom.dihedral_to_index

            rc = coords[c]
            rb = coords[b]
            ra = coords[a]

            e1 = rc - rb
            e1_norm = np.linalg.norm(e1)
            e1 = np.array([1.0, 0.0, 0.0]) if np.isclose(e1_norm, 0.0) else e1 / e1_norm

            ab = rb - ra
            normal = np.cross(ab, e1)
            normal_norm = np.linalg.norm(normal)
            if np.isclose(normal_norm, 0.0):
                trial = np.array([0.0, 0.0, 1.0])
                if np.linalg.norm(np.cross(e1, trial)) < 1e-8:
                    trial = np.array([0.0, 1.0, 0.0])
                normal = np.cross(e1, trial)
                normal_norm = np.linalg.norm(normal)
            e3 = normal / normal_norm

            e2 = np.cross(e3, e1)
            e2 = e2 / np.linalg.norm(e2)

            if atom.zmat_format == 1:
                uc = rb - rc
                ua = ra - rc
                uc_norm = np.linalg.norm(uc)
                ua_norm = np.linalg.norm(ua)
                uc = np.array([1.0, 0.0, 0.0]) if np.isclose(uc_norm, 0.0) else uc / uc_norm
                ua = np.array([0.0, 1.0, 0.0]) if np.isclose(ua_norm, 0.0) else ua / ua_norm

                q = ua - np.dot(ua, uc) * uc
                q_norm = np.linalg.norm(q)
                q = e2 if np.isclose(q_norm, 0.0) else q / q_norm
                nvec = np.cross(uc, q)
                nvec_norm = np.linalg.norm(nvec)
                nvec = e3 if np.isclose(nvec_norm, 0.0) else nvec / nvec_norm

                cos_theta = np.cos(theta)
                second_angle = atom.dihedral.to(atom_ureg.radian).magnitude
                cos_alpha = np.cos(second_angle)
                cos_ref = float(np.clip(np.dot(uc, ua), -1.0, 1.0))
                sin_ref = float(np.sqrt(max(0.0, 1.0 - cos_ref**2)))

                if np.isclose(sin_ref, 0.0):
                    y = 0.0
                    z_sq = max(0.0, 1.0 - cos_theta**2)
                else:
                    y = (cos_alpha - cos_ref * cos_theta) / sin_ref
                    z_sq = max(0.0, 1.0 - cos_theta**2 - y**2)
                z = float(np.sqrt(z_sq))

                direction = cos_theta * uc + y * q + z * nvec
                dir_norm = np.linalg.norm(direction)
                if np.isclose(dir_norm, 0.0):
                    direction = cos_theta * uc + y * q - z * nvec
                    dir_norm = np.linalg.norm(direction)
                direction = uc if np.isclose(dir_norm, 0.0) else direction / dir_norm
                coords[i] = rc + r * direction
            else:
                phi = atom.dihedral.to(atom_ureg.radian).magnitude
                coords[i] = rc + r * (
                    -np.cos(theta) * e1
                    + np.sin(theta) * np.cos(phi) * e2
                    + np.sin(theta) * np.sin(phi) * e3
                )

        coords = coords[~np.array([atom.is_dummy or atom.is_ghost for atom in self.atoms])]
        return coords * atom_ureg.angstrom

    def to_XYZ_block(self) -> str:
        symbols = self.get_symbols()
        coords = self.to_cartesian_coords().magnitude
        return f"{len(symbols)}\n\n" + "\n".join(
            [
                f"{symbol} {x:.6f} {y:.6f} {z:.6f}"
                for symbol, (x, y, z) in zip(symbols, coords, strict=True)
            ]
        )


class Energies(BaseDataClassWithUnit):
    default_units: ClassVar[dict[str, UnitLike]] = {
        "electronic_energy": atom_ureg.hartree,
        "scf_energy": atom_ureg.hartree,
        "mp2_energy": atom_ureg.hartree,
        "mp3_energy": atom_ureg.hartree,
        "mp4_energy": atom_ureg.hartree,
        "ccsd_energy": atom_ureg.hartree,
    }
    set_default_units: ClassVar[bool] = True

    # energies
    electronic_energy: PlainQuantity | None = Field(
        default=None,
        description="Electronic energy of the molecule, unit is `hartree`",
    )
    scf_energy: PlainQuantity | None = Field(
        default=None, description="SCF energy of the molecule, unit is `hartree`"
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
    ccsd_energy: PlainQuantity | None = Field(
        default=None, description="CCSD energy of the molecule, unit is `hartree`"
    )

    @property
    def energy(self) -> dict[str, PlainQuantity]:
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

    @computed_field(description="Total energy, unit is `hartree`")  # type: ignore[prop-decorator]
    @property
    def total_energy(self) -> PlainQuantity | None:
        keys = list(self.energy.keys())
        if len(keys) > 0:
            return self.energy[keys[0]].to(atom_ureg.hartree)
        else:
            return None

    def to_summary_dict(self, **kwargs) -> dict[tuple[str, str], Any]:
        return {
            ("Energy", f"{key} ({getattr(self, key).units})")
            if isinstance(getattr(self, key), PlainQuantity)
            else ("Energy", key): getattr(self, key).m
            if isinstance(getattr(self, key), PlainQuantity)
            else getattr(self, key)
            for key in self.model_dump(**kwargs)
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

    def to_summary_dict(self, **kwargs) -> dict[tuple[str, str], Any]:
        return {
            ("Thermal", f"{key} ({getattr(self, key).units})")
            if isinstance(getattr(self, key), PlainQuantity)
            else ("Thermal", key): getattr(self, key).m
            if isinstance(getattr(self, key), PlainQuantity)
            else getattr(self, key)
            for key in self.model_dump(**kwargs)
            if getattr(self, key) is not None
        }


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
    alpha_occupancy: bool | None = Field(default=None, description="alpha orbital occupancy")
    alpha_symmetry: str | None = Field(default=None, description="alpha orbital symmetry")
    beta_occupancy: bool | None = Field(default=None, description="beta orbital occupancy")
    beta_symmetry: str | None = Field(default=None, description="beta orbital symmetry")
    coefficient: np.ndarray | None = Field(default=None, description="coefficient of the orbital")

    def to_summary_dict(self, **kwargs) -> dict[tuple[str, str], Any]:
        return {
            ("Orbital", f"{key} ({getattr(self, key).units})")
            if isinstance(getattr(self, key), PlainQuantity)
            else ("Orbital", key): getattr(self, key).m
            if isinstance(getattr(self, key), PlainQuantity)
            else getattr(self, key)
            for key in self.model_dump(**kwargs)
            if getattr(self, key) is not None
        }


class MolecularOrbitals(BaseDataClassWithUnit, Sequence[MoleculeOrbital]):
    _index_: int = PrivateAttr(default=0)

    default_units: ClassVar[dict[str, UnitLike]] = {
        "alpha_energies": atom_ureg.hartree,
        "beta_energies": atom_ureg.hartree,
    }
    set_default_units: ClassVar[bool] = True

    # orbital energies
    electronic_state: str | None = Field(
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
    alpha_occupancies: list[bool | None] = Field(
        default_factory=list, description="alpha orbital occupancies"
    )
    beta_occupancies: list[bool | None] = Field(
        default_factory=list, description="beta orbital occupancies"
    )
    alpha_symmetries: list[str | None] = Field(
        default_factory=list, description="alpha orbital symmetries"
    )
    beta_symmetries: list[str | None] = Field(
        default_factory=list, description="beta orbital symmetries"
    )
    coefficients: list[np.ndarray | None] = Field(
        default_factory=list, description="coefficients of the orbitals"
    )

    @computed_field(description="HOMO orbital idx")  # type: ignore[prop-decorator]
    @property
    def HOMO_id(self) -> int | None:
        for i, alpha_occ in enumerate(self.alpha_occupancies):
            if alpha_occ is None:
                continue
            if i > 0 and not alpha_occ and self.alpha_occupancies[i - 1]:
                return i - 1
        return None

    @computed_field(description="LUMO orbital idx")  # type: ignore[prop-decorator]
    @property
    def LUMO_id(self) -> int | None:
        if self.HOMO_id is None:
            return None
        return self.HOMO_id + 1

    @computed_field(description="beta HOMO orbital idx")  # type: ignore[prop-decorator]
    @property
    def beta_HOMO_id(self) -> int | None:
        for i, beta_occ in enumerate(self.beta_occupancies):
            if beta_occ is None:
                continue
            if i > 0 and not beta_occ and self.beta_occupancies[i - 1]:
                return i - 1
        return None

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
            if alpha_occ and not beta_occ
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

    def __next__(self) -> MoleculeOrbital:
        if self._index_ >= len(self):
            raise StopIteration
        self._index_ += 1
        return self[self._index_ - 1]

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

    def to_summary_dict(self, **kwargs) -> dict[tuple[str, str], Any]:
        return {
            ("Orbitals", "electronic_state"): self.electronic_state,
            ("Orbitals", "HOMO_energy"): None if self.HOMO_energy is None else self.HOMO_energy.m,
            ("Orbitals", "LUMO_energy"): None if self.LUMO_energy is None else self.LUMO_energy.m,
            ("Orbitals", "HOMO-LUMO_gap"): None
            if self.HOMO_LUMO_gap is None
            else self.HOMO_LUMO_gap.m,
        }


# TODO: ready for NaturalAtomicOrbitals
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

    def to_summary_dict(self, **kwargs) -> dict[tuple[str, str], Any]:
        return {
            ("NaturalAtomicOrbital", "element"): self.element,
            ("NaturalAtomicOrbital", "atom_index"): self.atom_index,
            ("NaturalAtomicOrbital", "angular_momentum"): self.angular_momentum,
            ("NaturalAtomicOrbital", "ao_type"): self.ao_type,
            ("NaturalAtomicOrbital", "occupancy"): self.occupancy,
            ("NaturalAtomicOrbital", "energy"): None if self.energy is None else self.energy.m,
        }


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

    def to_summary_dict(self, **kwargs) -> dict[tuple[str, str], Any]:
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

    def to_summary_dict(self, **kwargs) -> dict[tuple[str, str], Any]:
        return {
            ("NaturalBondOrbital", "orbital_type"): self.orbital_type,
            ("NaturalBondOrbital", "sub_index"): self.sub_index,
            ("NaturalBondOrbital", "bonding_atoms_dict"): self.bonding_atoms_dict,
            ("NaturalBondOrbital", "occupancy"): self.occupancy,
            ("NaturalBondOrbital", "energy"): None if self.energy is None else self.energy.m,
            (
                "NaturalBondOrbital",
                "principal_elocalizations_dict",
            ): self.principal_elocalizations_dict,
        }


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

    def to_summary_dict(self, **kwargs) -> dict[tuple[str, str], Any]:
        return {}


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
    vibration_mode: NumpyQuantity = Field(
        default=np.array([[]]) * atom_ureg.angstrom,
        description="Vibration mode of each mode, unit is `angstrom`",
        exclude_if=lambda x: x.shape == (0, 0),
    )

    @computed_field()  # type: ignore[prop-decorator]
    @property
    def is_imaginary(self) -> bool:
        return bool(self.frequency is not None and cast(Any, self.frequency) < 0)

    def to_summary_dict(self, **kwargs) -> dict[tuple[str, str], Any]:
        return {
            ("Vibration", f"{key} ({getattr(self, key).units})")
            if isinstance(getattr(self, key), PlainQuantity)
            else ("Vibration", key): getattr(self, key).m
            if isinstance(getattr(self, key), PlainQuantity)
            else getattr(self, key)
            for key in self.model_dump(**kwargs)
            if getattr(self, key) is not None
        }

    def transform_orientation(
        self, transformation_matrix: np.ndarray, inverse: bool = False
    ) -> None:
        if inverse:
            mode = cast(Any, self.vibration_mode)
            self.vibration_mode = invert_transform_coords(mode.m, transformation_matrix) * mode.u
        else:
            mode = cast(Any, self.vibration_mode)
            self.vibration_mode = transform_coords(mode.m, transformation_matrix) * mode.u


class Vibrations(BaseDataClassWithUnit, Sequence[Vibration]):
    _index_: int = PrivateAttr(default=0)
    default_units: ClassVar[dict[str, UnitLike]] = {
        "frequency": atom_ureg.cm_1,
        "reduced_mass": atom_ureg.amu,
        "force_constant": atom_ureg.Unit("mdyne/angstrom"),
        "IR_intensity": atom_ureg.Unit("km/mol"),
        "vibration_mode": atom_ureg.angstrom,
    }
    set_default_units: ClassVar[bool] = True

    frequencies: NumpyQuantity = Field(
        default=np.array([]) * atom_ureg.cm_1,
        description="Frequency of each mode, unit is `cm^-1`",
        exclude_if=lambda x: len(x) == 0,
    )
    reduced_masses: NumpyQuantity = Field(
        default=np.array([]) * atom_ureg.amu,
        description="Reduced mass of each mode, unit is `amu`",
        exclude_if=lambda x: len(x) == 0,
    )
    force_constants: NumpyQuantity = Field(
        default=np.array([]) * atom_ureg.mdyne / atom_ureg.angstrom,
        description="Force constant of each mode, unit is `mdyne/angstrom`",
        exclude_if=lambda x: len(x) == 0,
    )
    IR_intensities: NumpyQuantity = Field(
        default=np.array([]) * atom_ureg.km / atom_ureg.mol,
        description="IR intensity of each mode, unit is `km/mol`",
        exclude_if=lambda x: len(x) == 0,
    )
    vibration_modes: list[NumpyQuantity] = Field(
        default_factory=list,
        description="Vibration mode of each mode, unit is `angstrom`",
        exclude_if=lambda x: len(x) == 0,
    )

    def __iter__(self) -> Iterator[Vibration]:  # type: ignore[override]
        for i in range(len(self)):
            yield self[i]

    def __next__(self) -> Vibration:
        if self._index_ >= len(self):
            raise StopIteration
        self._index_ += 1
        return self[self._index_ - 1]

    @overload
    def __getitem__(self, frameID: int) -> Vibration: ...
    @overload
    def __getitem__(self, frameID: slice) -> list[Vibration]: ...
    @overload
    def __getitem__(self, frameID: Sequence) -> list[Vibration]: ...
    def __getitem__(self, frameID: int | slice | Sequence) -> Vibration | list[Vibration]:
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
    def imaginary_idxs(self) -> list[int]:
        return [i for i, freq in enumerate(self) if freq.is_imaginary]

    @property
    def imaginary_vibrations(self) -> "Vibrations":
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
            }
        )

    def transform_orientation(
        self, transformation_matrix: np.ndarray, inverse: bool = False
    ) -> None:
        if inverse:
            self.vibration_modes = [
                invert_transform_coords(cast(Any, mode).m, transformation_matrix)
                * cast(Any, mode).u
                for mode in self.vibration_modes
            ]
        else:
            self.vibration_modes = [
                transform_coords(cast(Any, mode).m, transformation_matrix) * cast(Any, mode).u
                for mode in self.vibration_modes
            ]

    def to_summary_dict(self, **kwargs) -> dict[tuple[str, str], Any]:
        return {
            ("Vibration", "num_imaginary"): self.num_imaginary,
            ("Vibration", "num_vibrations"): len(self),
        }


class ChargeSpinPopulations(BaseDataClassWithUnit):
    # charge and spin populations
    mulliken_charges: list[float] = Field(
        default_factory=list, description="Mulliken charges", exclude_if=lambda x: len(x) == 0
    )
    mulliken_spins: list[float] = Field(
        default_factory=list,
        description="Mulliken spin densities",
        exclude_if=lambda x: len(x) == 0,
    )
    apt_charges: list[float] = Field(
        default_factory=list,
        description="Atomic polarizability tensor charges",
        exclude_if=lambda x: len(x) == 0,
    )
    lowdin_charges: list[float] = Field(
        default_factory=list, description="Lowdin charges", exclude_if=lambda x: len(x) == 0
    )
    hirshfeld_charges: list[float] = Field(
        default_factory=list, description="Hirshfeld charges", exclude_if=lambda x: len(x) == 0
    )
    hirshfeld_spins: list[float] = Field(
        default_factory=list, description="Hirshfeld spins", exclude_if=lambda x: len(x) == 0
    )
    hirshfeld_q_cm5: list[float] = Field(
        default_factory=list,
        description="Hirshfeld charges in cm5",
        exclude_if=lambda x: len(x) == 0,
    )
    npa_charges: list[float] = Field(
        default_factory=list, description="NPA charges", exclude_if=lambda x: len(x) == 0
    )
    npa_alpha_spin_densities: list[float] = Field(
        default_factory=list,
        description="NPA alpha spin densities",
        exclude_if=lambda x: len(x) == 0,
    )
    npa_beta_spin_densities: list[float] = Field(
        default_factory=list,
        description="NPA beta spin densities",
        exclude_if=lambda x: len(x) == 0,
    )

    @model_validator(mode="after")
    def validate_charge_spin_populations(self) -> Self:
        available_populations = [
            getattr(self, pop) for pop in self.model_fields_set if len(getattr(self, pop)) > 0
        ]
        assert all(len(pop) == len(available_populations[0]) for pop in available_populations), (
            "All populations must have the same length"
        )
        return self

    def to_summary_dict(self, **kwargs) -> dict[tuple[str, str], Any]:
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

    def to_summary_dict(self, **kwargs) -> dict[tuple[str, str], Any]:
        return {("TotalSpin", key): value for key, value in self.to_unitless_dump(**kwargs).items()}


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
    polarizability_tensor: PlainQuantity | None = Field(
        default=np.array([]) * atom_ureg.bohr**3,
        description="Polarizability tensor",
        exclude_if=lambda x: (x is None) or (len(x) == 0),
    )
    electric_dipole_moment: PlainQuantity | None = Field(
        default=np.array([]) * atom_ureg.debye,
        description="Electric dipole moment, unit is `debye`",
        exclude_if=lambda x: (x is None) or (len(x) == 0),
    )
    dipole: PlainQuantity | None = Field(
        default=np.array([]) * atom_ureg.debye,
        description="Dipole moment, unit is `debye`",
        exclude_if=lambda x: (x is None) or (len(x) == 0),
    )
    quadrupole: PlainQuantity | None = Field(
        default=np.array([]) * atom_ureg.debye * atom_ureg.angstrom,
        description="Quadrupole moment, unit is `debye*angstrom`",
        exclude_if=lambda x: (x is None) or (len(x) == 0),
    )
    traceless_quadrupole: PlainQuantity | None = Field(
        default=np.array([]) * atom_ureg.debye * atom_ureg.angstrom,
        description="Traceless quadrupole moment, unit is `debye*angstrom`",
        exclude_if=lambda x: (x is None) or (len(x) == 0),
    )
    octapole: PlainQuantity | None = Field(
        default=np.array([]) * atom_ureg.debye * atom_ureg.angstrom**2,
        description="Octapole moment, unit is `debye*angstrom**2`",
        exclude_if=lambda x: (x is None) or (len(x) == 0),
    )
    hexadecapole: PlainQuantity | None = Field(
        default=np.array([]) * atom_ureg.debye * atom_ureg.angstrom**3,
        description="Hexadecapole moment, unit is `debye*angstrom**3`",
        exclude_if=lambda x: (x is None) or (len(x) == 0),
    )

    def to_summary_dict(self, **kwargs) -> dict[tuple[str, str], Any]:
        return {
            ("Polarizability", f"{key} ({getattr(self, key).units})")
            if isinstance(getattr(self, key), PlainQuantity)
            else ("Polarizability", key): getattr(self, key).m
            if isinstance(getattr(self, key), PlainQuantity)
            else getattr(self, key)
            for key in self.model_dump(**kwargs)
            if getattr(self, key) is not None
        }


class BondOrders(BaseDataClassWithUnit):
    wiberg_bond_order: np.ndarray = Field(
        default=np.array([[]]),
        description="Wiberg bond order",
        exclude_if=lambda x: x.shape == (0, 0),
    )
    mo_bond_order: np.ndarray = Field(
        default=np.array([[]]),
        description="MO bond order, ∑[i∈A]∑[j∈B]P(i,j)",
        exclude_if=lambda x: x.shape == (0, 0),
    )
    mayer_bond_order: np.ndarray = Field(
        default=np.array([[]]),
        description="MAYER POPULATION ANALYSIS bond order, ∑[i∈A]∑[j∈B]P(i,j)",
        exclude_if=lambda x: x.shape == (0, 0),
    )
    atom_atom_overlap_bond_order: np.ndarray = Field(
        default=np.array([[]]),
        description="Atom-atom overlap bond order",
        exclude_if=lambda x: x.shape == (0, 0),
    )
    nbo_bond_order: np.ndarray = Field(
        default=np.array([[]]),
        description="NBO bond order",
        exclude_if=lambda x: x.shape == (0, 0),
    )
    nbo_bond_order_for_alpha_spin: np.ndarray = Field(
        default=np.array([[]]),
        description="NBO bond order for alpha spin",
        exclude_if=lambda x: x.shape == (0, 0),
    )
    nbo_bond_order_for_beta_spin: np.ndarray = Field(
        default=np.array([[]]),
        description="NBO bond order for beta spin",
        exclude_if=lambda x: x.shape == (0, 0),
    )

    def to_summary_dict(self, **kwargs) -> dict[tuple[str, str], Any]:
        return {}


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

    def to_summary_dict(self, **kwargs) -> dict[tuple[str, str], Any]:
        return {
            ("Dispersion", f"{key} ({getattr(self, key).units})")
            if isinstance(getattr(self, key), PlainQuantity)
            else ("Dispersion", key): getattr(self, key).m
            if isinstance(getattr(self, key), PlainQuantity)
            else getattr(self, key)
            for key in self.model_dump(**kwargs)
            if getattr(self, key) is not None
        }


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

    def to_summary_dict(self, **kwargs) -> dict[tuple[str, str], Any]:
        return {}


class GeometryOptimizationStatus(BaseDataClassWithUnit):
    """
    Geometry optimization status.
    """

    geometry_optimized: bool | None = Field(
        default=None,
        description="Whether the geometry has been optimized",
        exclude_if=lambda x: x is None,
    )
    energy_change_threshold: float | None = Field(
        default=None,
        description="Energy change threshold",
        exclude_if=lambda x: x is None,
    )
    rms_force_threshold: float | None = Field(
        default=None,
        description="RMS force threshold in internal some programs use gradient, which has the same absolute value",
        exclude_if=lambda x: x is None,
    )
    max_force_threshold: float | None = Field(
        default=None,
        description="Maximum force threshold in internal some programs use gradient, which has the same absolute value",
        exclude_if=lambda x: x is None,
    )
    rms_displacement_threshold: float | None = Field(
        default=None,
        description="RMS displacement threshold in internal",
        exclude_if=lambda x: x is None,
    )
    max_displacement_threshold: float | None = Field(
        default=None,
        description="Maximum displacement threshold in internal",
        exclude_if=lambda x: x is None,
    )
    energy_change: float = Field(
        default=float("inf"),
        description="Energy change",
    )
    rms_force: float = Field(
        default=float("inf"),
        description="RMS force some programs use gradient, which has the same absolute value",
    )
    max_force: float = Field(
        default=float("inf"),
        description="Maximum force some programs use gradient, which has the same absolute value",
    )
    rms_displacement: float = Field(default=float("inf"), description="RMS displacement")
    max_displacement: float = Field(default=float("inf"), description="Maximum displacement")

    @computed_field()  # type: ignore[prop-decorator]
    @property
    def energy_change_converged(self) -> bool | None:
        """
        Whether the energy change has converged.
        """
        if self.energy_change_threshold is None:
            return None
        return self.energy_change < self.energy_change_threshold

    @computed_field()  # type: ignore[prop-decorator]
    @property
    def rms_force_converged(self) -> bool | None:
        """
        Whether the RMS force has converged.
        """
        if self.rms_force_threshold is None:
            return None
        return self.rms_force < self.rms_force_threshold

    @computed_field()  # type: ignore[prop-decorator]
    @property
    def max_force_converged(self) -> bool | None:
        """
        Whether the maximum force has converged.
        """
        if self.max_force_threshold is None:
            return None
        return self.max_force < self.max_force_threshold

    @computed_field()  # type: ignore[prop-decorator]
    @property
    def rms_displacement_converged(self) -> bool | None:
        """
        Whether the RMS displacement has converged.
        """
        if self.rms_displacement_threshold is None:
            return None
        return self.rms_displacement < self.rms_displacement_threshold

    @computed_field()  # type: ignore[prop-decorator]
    @property
    def max_displacement_converged(self) -> bool | None:
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
    ) -> np.ndarray[Any, np.dtype[np.bool_]] | np.ndarray[Any, np.dtype[np.floating[Any]]]:
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
        metrics_list = list(metrics)
        df = pd.DataFrame(
            index=pd.Index(metrics_list),
            columns=pd.Index(["value", "threshold", "converged"]),
        )
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

    def to_summary_dict(self, **kwargs) -> dict[tuple[str, str], Any]:
        return {
            (
                "GeometryOptimizationStatus",
                "geometry_optimized",
            ): self.geometry_optimized,
            (
                "GeometryOptimizationStatus",
                "energy_change_converged",
            ): self.energy_change_converged,
            (
                "GeometryOptimizationStatus",
                "rms_force_converged",
            ): self.rms_force_converged,
            (
                "GeometryOptimizationStatus",
                "max_force_converged",
            ): self.max_force_converged,
            (
                "GeometryOptimizationStatus",
                "rms_displacement_converged",
            ): self.rms_displacement_converged,
            (
                "GeometryOptimizationStatus",
                "max_displacement_converged",
            ): self.max_displacement_converged,
        }


class Status(BaseDataClassWithUnit):
    scf_converged: bool = Field(default=False, description="Whether the SCF has converged")
    normal_terminated: bool = Field(
        default=False, description="Whether the calculation has terminated normally"
    )

    def to_summary_dict(self, **kwargs) -> dict[tuple[str, str], Any]:
        return {
            ("Status", key): getattr(self, key).m
            if isinstance(getattr(self, key), PlainQuantity)
            else getattr(self, key)
            for key in self.model_dump(**kwargs)
            if getattr(self, key) is not None
        }


class ShieldingTensor(BaseDataClassWithUnit):
    """
    Shielding tensor data.
    """

    default_units: ClassVar[dict[str, UnitLike]] = {"shielding_tensor": atom_ureg.ppm}
    set_default_units: ClassVar[bool] = True

    atom: str = Field(default="", description="Atom element")
    shielding_tensor: NumpyQuantity = Field(
        default=np.zeros((3, 3)) * atom_ureg.ppm,
        description="Shielding tensor, unit is `ppm`",
    )

    @computed_field()  # type: ignore[prop-decorator]
    @property
    def isotropic(self) -> PlainQuantity:
        return np.mean(np.trace(self.shielding_tensor))

    @computed_field()  # type: ignore[prop-decorator]
    @property
    def anisotropy(self) -> PlainQuantity:
        eigenvalues = np.linalg.eigvals(self.shielding_tensor)
        # 按照 Haeberlen 约定排序: |s_zz - s_iso| >= |s_xx - s_iso| >= |s_yy - s_iso|
        s_iso = np.mean(eigenvalues)
        sorted_eigs = sorted(eigenvalues, key=lambda x: abs(x - s_iso), reverse=True)
        return sorted_eigs[0] - (sorted_eigs[1] + sorted_eigs[2]) / 2

    def to_summary_dict(self, **kwargs) -> dict[tuple[str, str], Any]:
        return {
            (
                "ShieldingTensor",
                key,
            ): getattr(self, key).m
            if isinstance(getattr(self, key), PlainQuantity)
            else getattr(self, key)
            for key in self.model_dump(**kwargs)
            if getattr(self, key) is not None
        }


class NMR(BaseDataClassWithUnit):
    default_units: ClassVar[dict[str, UnitLike]] = {
        "spin_spin_coupling_k": atom_ureg.Hz,
        "spin_spin_coupling_j": atom_ureg.Hz,
    }
    set_default_units: ClassVar[bool] = True

    shielding_tensors: list[ShieldingTensor] = Field(
        default_factory=list, description="NMR shielding tensors"
    )
    spin_spin_coupling_k: NumpyQuantity | None = Field(
        default=None,
        description="Spin-spin coupling constant, unit is `Hz`",
    )
    spin_spin_coupling_j: NumpyQuantity | None = Field(
        default=None,
        description="Spin-spin coupling constant, unit is `Hz`",
    )

    def to_summary_dict(self, **kwargs) -> dict[tuple[str, str], Any]:
        return {
            ("NMR", key): getattr(self, key).m
            if isinstance(getattr(self, key), PlainQuantity)
            else getattr(self, key)
            for key in self.model_dump(**kwargs)
            if getattr(self, key) is not None
        }


class ImplicitSolvation(BaseDataClassWithUnit):
    solvent: str | None = Field(
        default=None,
        description="Solvent used in the QM calculation",
        exclude_if=lambda x: x is None,
    )
    solvent_model: str | None = Field(
        default=None,
        description="Solvent model used in the QM calculation",
        exclude_if=lambda x: x is None,
    )
    atomic_radii: str | None = Field(
        default=None,
        description="Atomic radii used in the QM calculation",
        exclude_if=lambda x: x is None,
    )
    solvent_epsilon: float | None = Field(
        default=None,
        description="Solvent dielectric constant used in the QM calculation",
        exclude_if=lambda x: x is None,
    )
    solvent_epsilon_infinite: float | None = Field(
        default=None,
        description="Solvent epsilon infinite used in the QM calculation",
        exclude_if=lambda x: x is None,
    )

    def to_summary_dict(self, **kwargs) -> dict[tuple[str, str], Any]:
        return {
            ("ImplicitSolvation", key): getattr(self, key).m
            if isinstance(getattr(self, key), PlainQuantity)
            else getattr(self, key)
            for key in self.model_dump(**kwargs)
            if getattr(self, key) is not None
        }
