"""Transition-state vibration and endpoint analysis.

The public frame models keep the historical convenience methods, while this
module owns the numerical sampling and topology-voting algorithms.  Callers
provide the frame's vibration runner and molecule-building dependencies so the
analysis layer stays independent from parser, animation, and process-pool
implementation details.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable, Sequence
from typing import Any, Literal, Protocol

import numpy as np
import numpy.typing as npt
from rdkit import Chem

from molop.io.base_models.data_classes.vibrations import Vibration, Vibrations
from molop.utils.types import PintArrayNx3, RdMol


VibrationSamplingMethod = Literal["amplitude", "harmonic_potential"]
MoleculeFactory = Callable[..., Any]
CrowdingChecker = Callable[[RdMol], bool]
VibrationRunner = Callable[..., Iterable[Any]]
EndpointSideSampler = Callable[..., RdMol | None]


class _VibrationFrame(Protocol):
    """Minimal frame contract required by the analysis algorithms."""

    @property
    def atoms(self) -> Sequence[int]: ...

    coords: PintArrayNx3
    charge: int
    multiplicity: int
    vibrations: Vibrations | None

    @property
    def atom_symbols(self) -> Sequence[str]: ...


def topology_frequency_key(rdmol: RdMol) -> str | bytes:
    """Return a conformer-independent key for TS endpoint voting."""

    try:
        smiles = Chem.MolToSmiles(rdmol, canonical=True, isomericSmiles=False)
        if smiles:
            return smiles
    except Exception:
        pass
    topology = Chem.Mol(rdmol)
    topology.RemoveAllConformers()
    return topology.ToBinary()


def most_frequent_topology(candidates: Sequence[RdMol], *, side: str) -> RdMol:
    """Select a side's topology mode and retain its largest-amplitude conformer."""

    if not candidates:
        raise ValueError(f"Failed to reconstruct any {side}-space TS endpoint candidates")

    # Candidates arrive in ascending amplitude order. Updating the representative
    # on every hit preserves the largest-amplitude conformer for the winning graph.
    grouped: dict[str | bytes, tuple[int, int, RdMol]] = {}
    for index, rdmol in enumerate(candidates):
        key = topology_frequency_key(rdmol)
        count = grouped[key][0] + 1 if key in grouped else 1
        grouped[key] = (count, index, Chem.Mol(rdmol))

    _, _, representative = max(grouped.values(), key=lambda item: (item[0], item[1]))
    return representative


def sample_vibration_amplitudes(
    min_ratio: float,
    max_ratio: float,
    steps: int,
    sampling_method: VibrationSamplingMethod,
) -> npt.NDArray[np.float64]:
    """Generate positive mode amplitudes for TS endpoint sampling.

    For a harmonic mode, the potential-energy magnitude is proportional to the
    square of the displacement amplitude, including for an imaginary TS mode
    when the curvature's magnitude is used. Therefore equal-energy spacing is
    obtained by spacing squared amplitudes linearly.
    """

    if sampling_method == "amplitude":
        return np.linspace(
            min_ratio,
            max_ratio,
            num=steps,
            endpoint=True,
            dtype=np.float64,
        )
    if sampling_method == "harmonic_potential":
        return np.sqrt(np.linspace(min_ratio**2, max_ratio**2, num=steps, endpoint=True))
    raise ValueError(
        "Unsupported sampling_method: "
        f"{sampling_method!r}. Use 'amplitude' or 'harmonic_potential'."
    )


def bond_change_pairs(first: RdMol, second: RdMol) -> set[tuple[int, int]]:
    """Return atom pairs whose bond presence or type differs between two graphs."""

    if first.GetNumAtoms() != second.GetNumAtoms():
        return set()

    def bond_types(molecule: RdMol) -> dict[tuple[int, int], Chem.BondType]:
        result: dict[tuple[int, int], Chem.BondType] = {}
        for bond in molecule.GetBonds():
            start_atom_idx = bond.GetBeginAtomIdx()
            end_atom_idx = bond.GetEndAtomIdx()
            if start_atom_idx > end_atom_idx:
                start_atom_idx, end_atom_idx = end_atom_idx, start_atom_idx
            result[(start_atom_idx, end_atom_idx)] = bond.GetBondType()
        return result

    first_bonds = bond_types(first)
    second_bonds = bond_types(second)
    return {
        atom_pair
        for atom_pair in first_bonds.keys() | second_bonds.keys()
        if first_bonds.get(atom_pair) != second_bonds.get(atom_pair)
    }


def bond_change_atom_indices(first: RdMol, second: RdMol) -> set[int]:
    """Return atoms incident to bonds that differ between two endpoint graphs."""

    return {
        atom_index for atom_pair in bond_change_pairs(first, second) for atom_index in atom_pair
    }


def _coordinate_rdmol(frame: _VibrationFrame, coordinates: npt.NDArray[np.float64]) -> RdMol | None:
    return Chem.MolFromXYZBlock(
        f"{len(frame.atoms)}\n"
        + f"charge {frame.charge} multiplicity {frame.multiplicity}\n"
        + "\n".join(
            [
                f"{Chem.Atom(atom).GetSymbol():10s}{x:10.5f}{y:10.5f}{z:10.5f}"
                for atom, (x, y, z) in zip(frame.atoms, coordinates, strict=True)
            ]
        )
    )


def generate_vibration_geometries(
    frame: _VibrationFrame,
    vibration_id: int | None = None,
    vibration: Vibration | None = None,
    *,
    ratio: float = 1.75,
    steps: int = 7,
    molecule_factory: MoleculeFactory,
    crowding_checker: CrowdingChecker,
) -> list[Any]:
    """Generate coordinate-only molecules displaced along one vibration mode."""

    if vibration is None:
        if vibration_id is None:
            vibration_id = 0
        if frame.vibrations is None:
            raise ValueError("No vibrations found in this frame")
        if vibration_id < 0 or vibration_id >= len(frame.vibrations):
            raise IndexError(f"Invalid vibration id {vibration_id}")
        vibration = frame.vibrations[vibration_id]
    assert vibration.vibration_mode.m.shape == frame.coords.m.shape, "Invalid vibration mode"

    molecules: list[Any] = []
    for displacement in np.linspace(-ratio, ratio, num=steps, endpoint=True):
        coordinates = np.asarray(
            frame.coords.m - vibration.vibration_mode.m * displacement,
            dtype=float,
        )
        coordinate_rdmol = _coordinate_rdmol(frame, coordinates)
        if coordinate_rdmol is None or not crowding_checker(coordinate_rdmol):
            continue
        molecules.append(
            molecule_factory(
                atom_symbols=frame.atom_symbols,
                coords=coordinates,
                charge=frame.charge,
                multiplicity=frame.multiplicity,
            )
        )
    return molecules


def infer_possible_ts_endpoints(
    frame: _VibrationFrame,
    *,
    show_3D: bool = False,
    min_ratio: float = 0.6,
    max_ratio: float = 1.4,
    steps: int = 8,
    sampling_method: VibrationSamplingMethod = "harmonic_potential",
    vibration_runner: VibrationRunner,
) -> tuple[RdMol, RdMol]:
    """Infer stable-side TS endpoint topologies by voting over mode displacements."""

    if not np.isfinite(min_ratio) or min_ratio <= 0:
        raise ValueError("min_ratio must be a finite value greater than 0")
    if not np.isfinite(max_ratio) or max_ratio < min_ratio:
        raise ValueError("max_ratio must be finite and greater than or equal to min_ratio")
    if steps < 1:
        raise ValueError("steps must be >= 1")

    amplitudes = sample_vibration_amplitudes(
        min_ratio,
        max_ratio,
        steps,
        sampling_method,
    )
    selected_sides: list[RdMol] = []
    for side_name, direction in (("negative", -1.0), ("positive", 1.0)):
        side_candidates: list[RdMol] = []
        for amplitude in amplitudes:
            # With one vibration step, ``ratio`` selects the signed extreme.
            displaced = vibration_runner(ratio=direction * float(amplitude), steps=1)
            for molecule in displaced:
                if (rdmol := molecule.rdmol) is not None and getattr(
                    molecule, "topology_reconstruction_status", None
                ) != "suspicious_fallback":
                    side_candidates.append(rdmol)
                    break
        selected_sides.append(most_frequent_topology(side_candidates, side=side_name))

    # The side with more disconnected fragments is treated as the precursor.
    # Equal fragment counts retain the deterministic negative/positive ordering.
    selected_sides.sort(key=lambda mol: len(Chem.GetMolFrags(mol)), reverse=True)
    reactant_rdmol, product_rdmol = selected_sides
    if not show_3D:
        reactant_rdmol.RemoveAllConformers()
        product_rdmol.RemoveAllConformers()
    return reactant_rdmol, product_rdmol


def resample_ts_endpoint_side(
    frame: _VibrationFrame,
    *,
    side: str,
    direction: float,
    endpoint: RdMol,
    fixed_atom_indices: set[int],
    amplitudes: npt.NDArray[np.floating[Any]],
    vibration_id_resolver: Callable[[], int],
    molecule_factory: MoleculeFactory,
    crowding_checker: CrowdingChecker,
) -> RdMol | None:
    """Resample one TS mode direction while holding reaction-center atoms fixed."""

    if endpoint.GetNumConformers() == 0:
        raise ValueError(
            "Additional TS endpoint sampling requires endpoint conformers; "
            "call possible_pre_post_ts(show_3D=True) first."
        )

    endpoint_coords = np.asarray(endpoint.GetConformer().GetPositions(), dtype=float)
    ts_coords = np.asarray(frame.coords.m, dtype=float)
    vibrations = frame.vibrations
    if vibrations is None:
        raise ValueError("Additional TS endpoint sampling requires vibration data")
    vibration = vibrations[vibration_id_resolver()]
    mode_coords = np.asarray(vibration.vibration_mode.m, dtype=float)
    if endpoint_coords.shape != ts_coords.shape or mode_coords.shape != ts_coords.shape:
        raise ValueError("TS endpoint and vibration coordinates must have matching shapes")

    fixed_indices = sorted(fixed_atom_indices)
    side_candidates: list[RdMol] = []
    for amplitude in amplitudes:
        displaced_coords = np.array(
            ts_coords + direction * mode_coords * float(amplitude),
            dtype=float,
            copy=True,
        )
        displaced_coords[fixed_indices, :] = endpoint_coords[fixed_indices, :]
        coordinate_rdmol = _coordinate_rdmol(frame, displaced_coords)
        if coordinate_rdmol is None or not crowding_checker(coordinate_rdmol):
            continue

        molecule = molecule_factory(
            atom_symbols=frame.atom_symbols,
            coords=displaced_coords,
            charge=frame.charge,
            multiplicity=frame.multiplicity,
        )
        if (rdmol := molecule.rdmol) is not None and getattr(
            molecule, "topology_reconstruction_status", None
        ) != "suspicious_fallback":
            side_candidates.append(rdmol)

    if not side_candidates:
        return None
    return most_frequent_topology(side_candidates, side=f"{side}-additional")


def infer_additional_ts_endpoints(
    frame: _VibrationFrame,
    pre_rdmol: RdMol,
    post_rdmol: RdMol,
    *,
    min_ratio: float = 0.6,
    max_ratio: float = 1.4,
    steps: int = 8,
    sampling_method: VibrationSamplingMethod = "harmonic_potential",
    vibration_id_resolver: Callable[[], int],
    side_sampler: EndpointSideSampler,
) -> tuple[RdMol, RdMol]:
    """Generate additional TS endpoints with reaction-center atoms fixed."""

    if not np.isfinite(min_ratio) or min_ratio <= 0:
        raise ValueError("min_ratio must be a finite value greater than 0")
    if not np.isfinite(max_ratio) or max_ratio < min_ratio:
        raise ValueError("max_ratio must be finite and greater than or equal to min_ratio")
    if steps < 1:
        raise ValueError("steps must be >= 1")
    if pre_rdmol.GetNumAtoms() != post_rdmol.GetNumAtoms() or pre_rdmol.GetNumAtoms() != len(
        frame.atoms
    ):
        raise ValueError("TS endpoints must contain the same atoms as the source frame")

    amplitudes = sample_vibration_amplitudes(
        min_ratio,
        max_ratio,
        steps,
        sampling_method,
    )
    changed_atom_indices = bond_change_atom_indices(pre_rdmol, post_rdmol)
    if not changed_atom_indices:
        return pre_rdmol, post_rdmol

    ts_coords = np.asarray(frame.coords.m, dtype=float)
    vibrations = frame.vibrations
    if vibrations is None:
        raise ValueError("Additional TS endpoint sampling requires vibration data")
    vibration = vibrations[vibration_id_resolver()]
    mode_coords = np.asarray(vibration.vibration_mode.m, dtype=float)
    if mode_coords.shape != ts_coords.shape:
        raise ValueError("TS vibration mode and source coordinates must have matching shapes")

    endpoint_directions: list[float] = []
    for endpoint in (pre_rdmol, post_rdmol):
        if endpoint.GetNumConformers() == 0:
            raise ValueError(
                "Additional TS endpoint sampling requires endpoint conformers; "
                "call possible_pre_post_ts(show_3D=True) first."
            )
        endpoint_coords = np.asarray(endpoint.GetConformer().GetPositions(), dtype=float)
        if endpoint_coords.shape != ts_coords.shape:
            raise ValueError("TS endpoint and source coordinates must have matching shapes")
        projection = float(np.sum((endpoint_coords - ts_coords) * mode_coords))
        if np.isclose(projection, 0.0):
            raise ValueError("Could not determine the vibration direction of a TS endpoint")
        endpoint_directions.append(1.0 if projection > 0.0 else -1.0)

    additional_endpoints: list[RdMol] = []
    for side, direction, endpoint in zip(
        ("pre", "post"),
        endpoint_directions,
        (pre_rdmol, post_rdmol),
        strict=True,
    ):
        additional = side_sampler(
            side=side,
            direction=direction,
            endpoint=endpoint,
            fixed_atom_indices=changed_atom_indices,
            amplitudes=amplitudes,
        )
        additional_endpoints.append(additional if additional is not None else endpoint)

    additional_endpoints.sort(
        key=lambda molecule: len(Chem.GetMolFrags(molecule)),
        reverse=True,
    )
    return additional_endpoints[0], additional_endpoints[1]


__all__ = [
    "CrowdingChecker",
    "EndpointSideSampler",
    "MoleculeFactory",
    "VibrationRunner",
    "VibrationSamplingMethod",
    "bond_change_atom_indices",
    "bond_change_pairs",
    "generate_vibration_geometries",
    "infer_additional_ts_endpoints",
    "infer_possible_ts_endpoints",
    "most_frequent_topology",
    "resample_ts_endpoint_side",
    "sample_vibration_amplitudes",
    "topology_frequency_key",
]
