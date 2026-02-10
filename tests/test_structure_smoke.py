import sys
from collections.abc import Generator
from pathlib import Path

import pytest
from openbabel import pybel

from molop.structure.GraphReconstruction import xyz_to_omol_no_metal


sys.path.append(str(Path(__file__).resolve().parent))

from _helpers_structure import (
    build_methane_rdmol,
    build_methane_xyz_block,
    build_water_rdmol,
    build_water_xyz_block,
)


@pytest.fixture(autouse=True)
def clear_xyz_cache() -> Generator[None, None, None]:
    xyz_to_omol_no_metal.cache_clear()
    yield
    xyz_to_omol_no_metal.cache_clear()


def test_xyz_helpers_generate_valid_blocks() -> None:
    water_xyz = build_water_xyz_block()
    methane_xyz = build_methane_xyz_block()

    water_lines = water_xyz.splitlines()
    methane_lines = methane_xyz.splitlines()
    assert water_lines[0] == "3"
    assert methane_lines[0] == "5"
    assert len(water_lines) == 5
    assert len(methane_lines) == 7

    water_ob = pybel.readstring("xyz", water_xyz)
    methane_ob = pybel.readstring("xyz", methane_xyz)
    assert water_ob.OBMol.NumAtoms() == 3
    assert methane_ob.OBMol.NumAtoms() == 5


def test_rdkit_helpers_generate_deterministic_conformers() -> None:
    water_mol = build_water_rdmol()
    methane_mol = build_methane_rdmol()

    assert water_mol.GetNumAtoms() == 3
    assert methane_mol.GetNumAtoms() == 5
    assert water_mol.GetNumConformers() == 1
    assert methane_mol.GetNumConformers() == 1

    methane_conf = methane_mol.GetConformer()
    origin = methane_conf.GetAtomPosition(0)
    assert (origin.x, origin.y, origin.z) == pytest.approx((0.0, 0.0, 0.0), abs=1e-12)


def test_structure_backends_and_cache_smoke() -> None:
    recovered = xyz_to_omol_no_metal(build_water_xyz_block())

    assert recovered is not None
    assert recovered.OBMol.NumAtoms() == 3
    assert xyz_to_omol_no_metal.cache_info().currsize >= 1
