from collections.abc import Iterator, Sequence
from typing import Any, ClassVar, Generic, TypeVar, cast, overload

import numpy as np
import pandas as pd
from pint._typing import UnitLike
from pint.facets.numpy.quantity import NumpyQuantity
from pint.facets.plain import PlainQuantity
from pydantic import ConfigDict, Field, computed_field, model_validator
from typing_extensions import Self

from molop.unit import atom_ureg
from molop.utils.functions import invert_transform_coords, transform_coords

from .Bases import (
    BaseDataClassWithUnit,
    PropertyBundle,
    PropertyColumnValue,
    PropertyScalarValue,
    PropertyTable,
    PropertyTransition,
    SpectralBand,
    Spectrum,
    TensorProperty,
)


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

        return cls.from_atoms(atoms)

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


class QMResourceRequest(BaseDataClassWithUnit):
    default_units: ClassVar[dict[str, UnitLike]] = {"memory": atom_ureg.megabyte}
    set_default_units: ClassVar[bool] = True

    num_cpu: int | None = Field(default=None, description="Requested number of CPU cores")
    memory: PlainQuantity | None = Field(
        default=None, description="Requested memory, normalized to megabyte"
    )
    raw: str = Field(default="", description="Raw resource directives")
    options: dict[str, Any] = Field(
        default_factory=dict, description="Program-specific resource options"
    )


class QMBasisSet(BaseDataClassWithUnit):
    name: str = Field(default="", description="Basis-set name or expression")
    role: str = Field(
        default="orbital", description="Basis role, such as orbital, auxiliary, or ecp"
    )
    scope: str | None = Field(
        default=None, description="Scope for this basis assignment, such as global or atom"
    )
    atom_indices: list[int] = Field(
        default_factory=list, description="0-based atom indices covered by this basis"
    )
    element_symbols: list[str] = Field(
        default_factory=list, description="Element symbols covered by this basis"
    )
    raw: str = Field(default="", description="Raw basis-set directive")
    options: dict[str, Any] = Field(
        default_factory=dict, description="Program-specific basis-set options"
    )


class QMModelChemistry(BaseDataClassWithUnit):
    method_family: str | None = Field(default=None, description="Canonical method family")
    method: str | None = Field(default=None, description="Concrete method name")
    reference_method: str | None = Field(
        default=None, description="Reference wavefunction or method"
    )
    functional: str | None = Field(
        default=None, description="DFT or double-hybrid functional"
    )
    basis_set: str | None = Field(default=None, description="Primary orbital basis set")
    auxiliary_basis_set: str | None = Field(
        default=None, description="Primary auxiliary or fitting basis set"
    )
    basis_sets: list[QMBasisSet] = Field(
        default_factory=list, description="Structured global or local basis-set assignments"
    )
    dispersion_correction: str | None = Field(
        default=None, description="Empirical or nonlocal dispersion correction"
    )
    solvation_model: str | None = Field(default=None, description="Solvation model")
    solvent: str | None = Field(default=None, description="Solvent name")
    relativistic: str | None = Field(default=None, description="Relativistic treatment")
    spin_treatment: str | None = Field(default=None, description="Restricted/open-shell treatment")
    raw_keywords: str = Field(default="", description="Raw model chemistry keyword source")
    options: dict[str, Any] = Field(
        default_factory=dict, description="Program-specific model chemistry options"
    )


class QMTaskRequest(BaseDataClassWithUnit):
    task_type: str = Field(description="Canonical task type, such as sp, opt, freq, or scan")
    enabled: bool = Field(default=True, description="Whether this task is requested")
    derivative_order: int | None = Field(
        default=None, description="Requested derivative order when applicable"
    )
    target_state: int | None = Field(default=None, description="Target electronic state/root")
    transition_state: bool = Field(default=False, description="Whether TS optimization is requested")
    scan: bool = Field(default=False, description="Whether this task includes a coordinate scan")
    properties: list[str] = Field(default_factory=list, description="Requested properties")
    source_keywords: list[str] = Field(default_factory=list, description="Source keyword tokens")
    source_blocks: list[str] = Field(default_factory=list, description="Source input block names")
    options: dict[str, Any] = Field(
        default_factory=dict, description="Program-specific task options"
    )


class ExcitedStateRequest(BaseDataClassWithUnit):
    enabled: bool = Field(default=False, description="Whether excited-state treatment is requested")
    family: str | None = Field(
        default=None, description="Excited-state method family, such as TDDFT, CIS, or EOM-CCSD"
    )
    reference_method: str | None = Field(default=None, description="Reference method")
    nroots: int | None = Field(default=None, description="Requested number of roots")
    root: int | None = Field(default=None, description="Primary target root")
    secondary_root: int | None = Field(default=None, description="Secondary target root")
    roots: list[int] = Field(default_factory=list, description="Explicit requested roots")
    singlets: bool = Field(default=False, description="Whether singlet states are requested")
    triplets: bool = Field(default=False, description="Whether triplet states are requested")
    spin_flip: bool = Field(default=False, description="Whether spin-flip treatment is requested")
    follow_root: bool = Field(default=False, description="Whether root following is requested")
    properties: list[str] = Field(default_factory=list, description="Requested excited-state properties")
    source_blocks: list[str] = Field(default_factory=list, description="Source input block names")
    options: dict[str, Any] = Field(
        default_factory=dict, description="Program-specific excited-state options"
    )


class ActiveSpace(BaseDataClassWithUnit):
    electrons: int | None = Field(default=None, description="Number of active electrons")
    orbitals: int | None = Field(default=None, description="Number of active orbitals")
    roots: int | None = Field(default=None, description="Number of averaged or targeted roots")
    active_orbitals: list[int] = Field(
        default_factory=list, description="0-based active orbital indices"
    )
    inactive_orbitals: list[int] = Field(
        default_factory=list, description="0-based inactive orbital indices"
    )
    frozen_orbitals: list[int] = Field(
        default_factory=list, description="0-based frozen orbital indices"
    )
    raw: str = Field(default="", description="Raw active-space expression")
    options: dict[str, Any] = Field(
        default_factory=dict, description="Program-specific active-space options"
    )


class MultireferenceStateBlock(BaseDataClassWithUnit):
    multiplicity: int | None = Field(default=None, description="Spin multiplicity")
    irrep: str | None = Field(default=None, description="Irreducible representation label")
    nroots: int | None = Field(default=None, description="Number of roots in this block")
    excitations: str | None = Field(default=None, description="Excitation selection")
    refs: str | None = Field(default=None, description="Reference-space definition")
    active_space: ActiveSpace | None = Field(default=None, description="Active-space summary")
    raw_lines: list[str] = Field(default_factory=list, description="Raw state-block lines")
    options: dict[str, Any] = Field(
        default_factory=dict, description="Program-specific state-block options"
    )


class MultireferenceRequest(BaseDataClassWithUnit):
    enabled: bool = Field(default=False, description="Whether multi-reference treatment is requested")
    method: str | None = Field(default=None, description="Multi-reference method")
    reference_method: str | None = Field(default=None, description="Reference method")
    ci_type: str | None = Field(default=None, description="Configuration interaction type")
    active_space: ActiveSpace | None = Field(default=None, description="Global active-space summary")
    state_blocks: list[MultireferenceStateBlock] = Field(
        default_factory=list, description="Requested spin/symmetry/root blocks"
    )
    thresholds: dict[str, float] = Field(default_factory=dict, description="Numeric thresholds")
    corrections: list[str] = Field(default_factory=list, description="Requested corrections")
    source_blocks: list[str] = Field(default_factory=list, description="Source input block names")
    options: dict[str, Any] = Field(
        default_factory=dict, description="Program-specific multi-reference options"
    )


class ExplicitSolventRequest(BaseDataClassWithUnit):
    enabled: bool = Field(default=False, description="Whether explicit-solvent placement is requested")
    solvent_model: str | None = Field(default=None, description="Underlying implicit solvent model")
    solvent: str | None = Field(default=None, description="Solvent name chosen for explicit placement")
    solvent_file: str | None = Field(default=None, description="Custom solvent file, if provided")
    nsolv: int | None = Field(default=None, description="Number of explicit solvent molecules")
    cluster_mode: str | None = Field(default=None, description="SOLVATOR cluster mode")
    droplet: bool = Field(default=False, description="Whether droplet mode is enabled")
    radius: float | None = Field(default=None, description="Target droplet radius")
    fixsolute: bool = Field(default=True, description="Whether the solute is kept frozen")
    vacuumsearch: bool = Field(default=False, description="Whether vacuum search is enabled")
    randomsolv: bool = Field(default=False, description="Whether random placement is enabled")
    printlevel: str | None = Field(default=None, description="SOLVATOR print level")
    source_blocks: list[str] = Field(default_factory=list, description="Source input block names")
    options: dict[str, Any] = Field(
        default_factory=dict, description="Program-specific explicit-solvent options"
    )


class ElectronicConfiguration(BaseDataClassWithUnit):
    label: str | None = Field(default=None, description="Configuration label")
    coefficient: float | None = Field(default=None, description="CI or configuration coefficient")
    weight: float | None = Field(default=None, description="Configuration weight")
    occupation: list[float] = Field(default_factory=list, description="Orbital occupation pattern")
    orbital_indices: list[int] = Field(default_factory=list, description="0-based orbital indices")
    raw: str = Field(default="", description="Raw configuration text")


class ElectronicState(BaseDataClassWithUnit):
    default_units: ClassVar[dict[str, UnitLike]] = {
        "energy": atom_ureg.hartree,
        "excitation_energy": atom_ureg.eV,
        "transition_dipole": atom_ureg.debye,
    }
    set_default_units: ClassVar[bool] = True

    state_index: int | None = Field(default=None, description="0-based electronic state index")
    root: int | None = Field(default=None, description="Program root index")
    label: str | None = Field(default=None, description="Electronic-state label")
    multiplicity: int | None = Field(default=None, description="Spin multiplicity")
    spin: float | None = Field(default=None, description="Spin quantum number")
    irrep: str | None = Field(default=None, description="Irreducible representation")
    method: str | None = Field(default=None, description="Method used for this state")
    energy: PlainQuantity | None = Field(default=None, description="State total energy")
    excitation_energy: PlainQuantity | None = Field(default=None, description="Excitation energy")
    oscillator_strength: float | None = Field(default=None, description="Oscillator strength")
    transition_dipole: NumpyQuantity | None = Field(
        default=None, description="Transition dipole vector"
    )
    configurations: list[ElectronicConfiguration] = Field(
        default_factory=list, description="Dominant electronic configurations"
    )
    properties: dict[str, Any] = Field(
        default_factory=dict, description="Additional state-specific properties"
    )
    source: str | None = Field(default=None, description="Program/source section")


class ElectronicStates(BaseDataClassWithUnit, Sequence[ElectronicState]):
    states: list[ElectronicState] = Field(
        default_factory=list, description="Electronic states ordered by source output"
    )

    def __iter__(self) -> Iterator[ElectronicState]:  # type: ignore[override]
        return iter(self.states)

    def __len__(self) -> int:
        return len(self.states)

    @overload
    def __getitem__(self, index: int) -> ElectronicState: ...

    @overload
    def __getitem__(self, index: slice) -> list[ElectronicState]: ...

    def __getitem__(self, index: int | slice) -> ElectronicState | list[ElectronicState]:
        return self.states[index]

    @property
    def ground_state(self) -> ElectronicState | None:
        return self.states[0] if self.states else None

    @property
    def excited_states(self) -> list[ElectronicState]:
        return self.states[1:]

    def to_state_table(self) -> PropertyTable:
        columns: dict[str, PropertyColumnValue] = {
            "state_index": [
                idx if state.state_index is None else state.state_index
                for idx, state in enumerate(self.states)
            ],
            "root": [state.root for state in self.states],
            "label": [state.label for state in self.states],
            "multiplicity": [state.multiplicity for state in self.states],
            "spin": [state.spin for state in self.states],
            "irrep": [state.irrep for state in self.states],
            "method": [state.method for state in self.states],
            "oscillator_strength": [state.oscillator_strength for state in self.states],
        }

        energies = [state.energy for state in self.states]
        if any(energy is not None for energy in energies):
            columns["energy"] = np.array(
                [
                    np.nan if energy is None else energy.to(atom_ureg.hartree).magnitude
                    for energy in energies
                ]
            ) * atom_ureg.hartree

        excitation_energies = [state.excitation_energy for state in self.states]
        if any(excitation_energy is not None for excitation_energy in excitation_energies):
            columns["excitation_energy"] = np.array(
                [
                    np.nan
                    if excitation_energy is None
                    else excitation_energy.to(atom_ureg.eV).magnitude
                    for excitation_energy in excitation_energies
                ]
            ) * atom_ureg.eV

        return PropertyTable(
            columns=columns,
            row_labels=[
                state.label or str(state.root) if state.root is not None else str(idx)
                for idx, state in enumerate(self.states)
            ],
            metadata={"source": "ElectronicStates", "index_base": 0},
        )

    def to_transitions(self) -> list[PropertyTransition]:
        transitions: list[PropertyTransition] = []
        ground_state = self.ground_state
        initial_state = ground_state.state_index if ground_state else 0
        for idx, state in enumerate(self.states):
            if state is ground_state:
                continue
            final_state = state.state_index if state.state_index is not None else idx
            transitions.append(
                PropertyTransition(
                    label=state.label,
                    initial_state=initial_state,
                    final_state=final_state,
                    energy=state.excitation_energy,
                    oscillator_strength=state.oscillator_strength,
                    transition_dipole=state.transition_dipole,
                    properties={
                        "root": state.root,
                        "multiplicity": state.multiplicity,
                        "spin": state.spin,
                        "irrep": state.irrep,
                        "method": state.method,
                        "source": state.source,
                    },
                )
            )
        return transitions

    def to_spectrum(self, label: str = "Electronic transitions") -> Spectrum:
        transitions = self.to_transitions()
        return Spectrum(
            label=label,
            x_label="excitation energy",
            y_label="oscillator strength",
            transitions=transitions,
            metadata={"source": "ElectronicStates"},
        )

    def to_property_bundle(self) -> PropertyBundle:
        transitions = self.to_transitions()
        return PropertyBundle(
            tables={"electronic_states": self.to_state_table()} if self.states else {},
            spectra={"electronic_transitions": self.to_spectrum()} if transitions else {},
            transitions={"electronic": transitions} if transitions else {},
            metadata={"source": "ElectronicStates"},
        )


class MultireferenceResult(BaseDataClassWithUnit):
    method: str | None = Field(default=None, description="Multi-reference method")
    reference_method: str | None = Field(default=None, description="Reference method")
    ci_type: str | None = Field(default=None, description="Configuration interaction type")
    active_space: ActiveSpace | None = Field(default=None, description="Active-space summary")
    electronic_states: ElectronicStates | None = Field(
        default=None, description="State/root resolved multi-reference results"
    )
    corrections: dict[str, PlainQuantity | float | str | None] = Field(
        default_factory=dict, description="Energy corrections or named post-CI corrections"
    )
    diagnostics: list[str] = Field(default_factory=list, description="Parser or method diagnostics")
    properties: dict[str, Any] = Field(
        default_factory=dict, description="Additional multi-reference properties"
    )

    def to_property_bundle(self) -> PropertyBundle:
        scalar_properties: dict[str, PropertyScalarValue] = {
            f"{key}_correction": value
            for key, value in self.corrections.items()
            if value is not None
        }
        bundle = PropertyBundle(
            scalar_properties=scalar_properties,
            metadata={
                "source": "MultireferenceResult",
                "method": self.method,
                "reference_method": self.reference_method,
                "ci_type": self.ci_type,
                "diagnostics": self.diagnostics,
            },
        )
        if self.electronic_states is not None:
            states_bundle = self.electronic_states.to_property_bundle()
            bundle.tables.update(states_bundle.tables)
            bundle.spectra.update(states_bundle.spectra)
            bundle.transitions.update(states_bundle.transitions)
        return bundle


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
    }
    set_default_units: ClassVar[bool] = True

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

    @property
    def energy(self) -> dict[str, PlainQuantity]:
        energy_fields = (
            "ccsd_energy",
            "mp5_energy",
            "mp4_energy",
            "mp3_energy",
            "mp2_energy",
            "electronic_energy",
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
            }.items()
            if value is not None
        }
        return PropertyBundle(
            scalar_properties=scalar_properties,
            metadata={"source": "ThermalInformations"},
        )


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
    coefficients: list[np.ndarray | None] = Field(
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
        if isinstance(frameID, slice):
            return [self[idx] for idx in range(*frameID.indices(len(self.frequencies)))]
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
                    label=f"mode {mode_idx + 1}",
                    center=frequency,
                    intensity=intensity,
                    metadata={"mode_index": mode_idx},
                )
            )
        return Spectrum(
            label=label,
            x_label="wavenumber",
            y_label="IR intensity",
            bands=bands,
            metadata={"source": "Vibrations"},
        )

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

    @property
    def population_names(self) -> list[str]:
        return [
            name
            for name in type(self).model_fields
            if name.endswith(("charges", "spins", "densities", "q_cm5"))
            and len(getattr(self, name)) > 0
        ]

    def to_population_table(self, atom_symbols: Sequence[str] | None = None) -> PropertyTable:
        population_lengths = [
            len(getattr(self, population_name)) for population_name in self.population_names
        ]
        if population_lengths:
            num_atoms = population_lengths[0]
        elif atom_symbols is not None:
            num_atoms = len(atom_symbols)
        else:
            num_atoms = 0

        if atom_symbols is not None and len(atom_symbols) != num_atoms:
            raise ValueError(
                f"atom_symbols length {len(atom_symbols)} does not match population length {num_atoms}"
            )

        columns: dict[str, PropertyColumnValue] = {}
        if num_atoms:
            columns["atom_index"] = np.arange(num_atoms)
        if atom_symbols is not None:
            columns["atom_symbol"] = list(atom_symbols)
        for population_name in self.population_names:
            columns[population_name] = getattr(self, population_name)

        return PropertyTable(
            columns=columns,
            row_labels=[str(i) for i in range(num_atoms)],
            metadata={
                "source": "ChargeSpinPopulations",
                "index_base": 0,
                "population_names": self.population_names,
            },
        )

    def to_property_bundle(self, atom_symbols: Sequence[str] | None = None) -> PropertyBundle:
        table = self.to_population_table(atom_symbols=atom_symbols)
        return PropertyBundle(
            tables={"atomic_populations": table} if len(table) else {},
            metadata={"source": "ChargeSpinPopulations"},
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

    def to_property_bundle(self) -> PropertyBundle:
        scalar_properties: dict[str, PropertyScalarValue] = {
            key: value
            for key, value in {
                "spin_square": self.spin_square,
                "spin_quantum_number": self.spin_quantum_number,
            }.items()
            if value is not None
        }
        return PropertyBundle(
            scalar_properties=scalar_properties,
            metadata={"source": "TotalSpin"},
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

    def to_atomic_property_table(
        self, atom_symbols: Sequence[str] | None = None
    ) -> PropertyTable:
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
