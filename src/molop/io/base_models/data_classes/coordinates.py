"""Coordinate and internal-coordinate data classes."""

from __future__ import annotations

from collections.abc import Iterator, Sequence
from typing import ClassVar, Generic, TypeVar, overload

import numpy as np
from pint._typing import UnitLike
from pint.facets.numpy.quantity import NumpyQuantity
from pint.facets.plain import PlainQuantity
from pydantic import ConfigDict, Field
from typing_extensions import Self

from molop.io.base_models.Bases import BaseDataClassWithUnit
from molop.unit import atom_ureg


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


CoordinateAtomT = TypeVar("CoordinateAtomT")


class CoordinateContainer(BaseDataClassWithUnit, Generic[CoordinateAtomT]):
    """Generic ordered container for coordinate-like rows."""

    model_config = ConfigDict(arbitrary_types_allowed=True, extra="forbid")

    items: list[CoordinateAtomT] = Field(
        default_factory=list, description="Ordered coordinate-like items stored in the container"
    )

    @classmethod
    def from_atoms(cls, atoms: Sequence[CoordinateAtomT]) -> Self:
        return cls(items=list(atoms))

    @property
    def atoms(self) -> list[CoordinateAtomT]:
        return self.items

    @atoms.setter
    def atoms(self, value: Sequence[CoordinateAtomT]) -> None:
        self.items = list(value)

    def __iter__(self) -> Iterator[CoordinateAtomT]:  # type: ignore[override]
        return iter(self.items)

    def __len__(self) -> int:
        return len(self.items)

    @overload
    def __getitem__(self, index: int) -> CoordinateAtomT: ...

    @overload
    def __getitem__(self, index: slice) -> list[CoordinateAtomT]: ...

    def __getitem__(self, index: int | slice) -> CoordinateAtomT | list[CoordinateAtomT]:
        return self.items[index]

    def append(self, atom: CoordinateAtomT) -> None:
        self.items.append(atom)

    def extend(self, atoms: Sequence[CoordinateAtomT]) -> None:
        self.items.extend(atoms)

    def clear(self) -> None:
        self.items.clear()

    def iter_real_atoms(self) -> Iterator[CoordinateAtomT]:
        for atom in self.items:
            if getattr(atom, "is_dummy", False) or getattr(atom, "is_ghost", False):
                continue
            yield atom

    def real_atoms(self) -> list[CoordinateAtomT]:
        return list(self.iter_real_atoms())

    def get_symbols(self) -> list[str]:
        symbols: list[str] = []
        for atom in self.iter_real_atoms():
            symbol = getattr(atom, "symbol", "")
            if symbol:
                symbols.append(symbol)
        return symbols


class CoordinateParameter(BaseDataClassWithUnit):
    name: str = Field(description="Parameter name")
    raw_value: str = Field(default="", description="Original parameter expression")
    start: float | None = Field(default=None, description="Start or representative value")
    stop: float | None = Field(default=None, description="End value for scans")
    steps: int | None = Field(default=None, description="Number of scan points")

    @property
    def is_scan(self) -> bool:
        return self.stop is not None or self.steps is not None


class CoordinateParameters(CoordinateContainer[CoordinateParameter]):
    items: list[CoordinateParameter] = Field(
        default_factory=list,
        description="Structured coordinate parameters or scan variables",
    )

    @property
    def parameters(self) -> list[CoordinateParameter]:
        return self.items

    @parameters.setter
    def parameters(self, value: Sequence[CoordinateParameter]) -> None:
        self.items = list(value)

    def as_value_map(self) -> dict[str, float]:
        values: dict[str, float] = {}
        for parameter in self.items:
            if parameter.start is not None:
                values[parameter.name] = parameter.start
        return values


class InternalCoords(CoordinateContainer[AtomInInternalCoords]):
    items: list[AtomInInternalCoords] = Field(
        default_factory=list, description="Atoms in the internal coordinates"
    )

    @classmethod
    def from_cartesian_coords(cls, symbols: Sequence[str], coords: NumpyQuantity) -> Self:
        """Create internal coordinates from atom symbols and Cartesian coordinates."""

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

        return cls.from_atoms(atoms)

    def to_cartesian_coords(self) -> NumpyQuantity:
        """Convert internal coordinates to Cartesian coordinates."""

        n = len(self.atoms)
        if n == 0:
            return np.zeros((0, 3)) * atom_ureg.angstrom

        coords = np.zeros((n, 3), dtype=float)
        coords[0] = np.array([0.0, 0.0, 0.0])

        if n == 1:
            return coords * atom_ureg.angstrom

        atom = self.atoms[1]
        r = atom.distance.to(atom_ureg.angstrom).magnitude
        coords[1] = np.array([r, 0.0, 0.0])

        if n == 2:
            return coords * atom_ureg.angstrom

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
            f"{symbol} {x:.6f} {y:.6f} {z:.6f}"
            for symbol, (x, y, z) in zip(symbols, coords, strict=True)
        )


__all__ = [
    "AtomInInternalCoords",
    "CoordinateContainer",
    "CoordinateParameter",
    "CoordinateParameters",
    "InternalCoords",
]
