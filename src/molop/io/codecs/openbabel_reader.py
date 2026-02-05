from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import Any, cast

import numpy as np
from openbabel import pybel

from molop.io.codec_exceptions import ParseError
from molop.io.codec_types import ParseResult, ReaderCodec, StructureLevel
from molop.io.logic.coords_frame_models.XYZFileFrame import XYZFileFrameDisk
from molop.io.logic.coords_models.XYZFile import XYZFileDisk
from molop.unit import atom_ureg


_OPENBABEL_FORMAT_ID = "openbabel"
_FALLBACK_FORMATS = ("sdf", "mol", "mol2", "pdb", "xyz", "smi", "cml")
_registered_registries: set[int] = set()


def _candidate_formats(path: str | Path) -> list[str]:
    suffix = Path(path).suffix.lower().lstrip(".")
    candidates: list[str] = []
    if suffix and suffix in pybel.informats:
        candidates.append(suffix)
    for infmt in _FALLBACK_FORMATS:
        if infmt in pybel.informats and infmt not in candidates:
            candidates.append(infmt)
    return candidates


def _best_effort_charge(mol: pybel.Molecule) -> int:
    for attr in ("charge", "formalcharge"):
        if hasattr(mol, attr):
            try:
                return int(getattr(mol, attr))
            except (TypeError, ValueError):
                pass
    try:
        return int(mol.OBMol.GetTotalCharge())
    except Exception:
        return 0


def _best_effort_multiplicity(mol: pybel.Molecule) -> int:
    for attr in ("spin", "multiplicity"):
        if hasattr(mol, attr):
            try:
                value = int(getattr(mol, attr))
                if value:
                    return value
            except (TypeError, ValueError):
                pass
    try:
        value = int(mol.OBMol.GetTotalSpinMultiplicity())
        if value:
            return value
    except Exception:
        pass
    return 1


def _read_first_molecule(path: str, infmt: str) -> pybel.Molecule | None:
    try:
        mol_iter = cast(Iterator[pybel.Molecule], pybel.readfile(infmt, path))
    except Exception:
        return None
    try:
        for mol in mol_iter:
            if mol is None:
                continue
            if mol.atoms:
                return mol
    except Exception:
        return None
    return None


@dataclass
class _OpenBabelFallbackReader:
    format_id: str = _OPENBABEL_FORMAT_ID
    extensions: frozenset[str] = frozenset()
    priority: int = -10_000

    def read(self, path: str | Path, **kwargs: Any) -> ParseResult[object]:
        path_str = str(path)
        for infmt in _candidate_formats(path_str):
            mol = _read_first_molecule(path_str, infmt)
            if mol is None:
                continue
            atoms = [atom.atomicnum for atom in mol.atoms]
            coords = np.array([atom.coords for atom in mol.atoms], dtype=float)
            charge = _best_effort_charge(mol)
            multiplicity = _best_effort_multiplicity(mol)
            frame = XYZFileFrameDisk(
                atoms=atoms,
                coords=coords * atom_ureg.angstrom,
                charge=charge,
                multiplicity=multiplicity,
            )
            frame.file_path = path_str
            file = XYZFileDisk(
                charge=charge,
                multiplicity=multiplicity,
            )
            file.file_path = path_str
            file.append(frame)
            return ParseResult(
                value=file,
                level=StructureLevel.COORDS,
                detected_format=infmt,
            )
        raise ParseError(f"OpenBabel could not parse {path_str}.")


def register_openbabel_fallback_reader(registry=None) -> None:
    from molop.io.codec_registry import Registry, default_registry

    reg: Registry = default_registry if registry is None else registry
    reg_id = id(reg)
    if reg_id in _registered_registries:
        return
    _registered_registries.add(reg_id)

    def factory() -> ReaderCodec:
        return _OpenBabelFallbackReader()

    reg.register_openbabel_fallback(factory)


__all__ = ["register_openbabel_fallback_reader"]
