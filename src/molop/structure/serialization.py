"""Scientific structure serialization helpers."""

from __future__ import annotations

from collections.abc import Iterable, Sequence

from molop.utils.types import RdMol


def xyz_block(
    atom_symbols: Sequence[str],
    coordinates: Iterable[Sequence[float]],
) -> str:
    """Serialize atom symbols and Cartesian coordinates to MolOP's XYZ format."""

    return f"{len(atom_symbols)}\n\n" + "\n".join(
        f"{atom:10s}{x:18.10f}{y:18.10f}{z:18.10f}"
        for atom, (x, y, z) in zip(atom_symbols, coordinates, strict=True)
    )


def rdmol_to_xyz(rdmol: RdMol) -> str:
    """Serialize the first RDKit conformer to MolOP's XYZ format."""

    return xyz_block(
        [atom.GetSymbol() for atom in rdmol.GetAtoms()],
        rdmol.GetConformer().GetPositions(),
    )


__all__ = ["rdmol_to_xyz", "xyz_block"]
