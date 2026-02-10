"""
Author: TMJ
Date: 2026-02-10 14:13:08
LastEditors: TMJ
LastEditTime: 2026-02-10 15:16:32
Description: 请填写简介
"""

import sys
from collections.abc import Generator
from pathlib import Path

import pytest
from rdkit import Chem

import molop.structure as structure_public_api
from molop.structure import xyz_to_rdmol
from molop.structure.GraphReconstruction import xyz_to_omol_no_metal


sys.path.append(str(Path(__file__).resolve().parent))

from _helpers_structure import build_methane_xyz_block, build_water_xyz_block


@pytest.fixture(autouse=True)
def clear_xyz_cache() -> Generator[None, None, None]:
    xyz_to_omol_no_metal.cache_clear()
    yield
    xyz_to_omol_no_metal.cache_clear()


@pytest.mark.parametrize(
    "invalid_xyz",
    [
        "1\ninvalid coordinate\nH nan 0.0 0.0\n",
        "1\ninvalid coordinate\nH inf 0.0 0.0\n",
    ],
)
def test_xyz_to_rdmol_returns_none_on_invalid_xyz(invalid_xyz: str) -> None:
    assert xyz_to_rdmol(invalid_xyz) is None


@pytest.mark.parametrize(
    ("xyz_block", "expected_atoms", "expected_bonds"),
    [
        (build_water_xyz_block(), 3, 2),
        (build_methane_xyz_block(), 5, 4),
    ],
)
@pytest.mark.parametrize("make_dative", [True, False])
def test_xyz_to_rdmol_returns_valid_mol_for_public_wrapper(
    xyz_block: str,
    expected_atoms: int,
    expected_bonds: int,
    make_dative: bool,
) -> None:
    mol = xyz_to_rdmol(xyz_block, make_dative=make_dative)

    assert mol is not None
    assert mol.GetNumAtoms() == expected_atoms
    assert mol.GetNumBonds() == expected_bonds
    Chem.SanitizeMol(mol)


def test_xyz_to_rdmol_returns_none_when_omol_to_rdmol_returns_none(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(structure_public_api, "omol_to_rdmol", lambda *_args, **_kwargs: None)

    assert xyz_to_rdmol(build_water_xyz_block(), make_dative=False) is None
