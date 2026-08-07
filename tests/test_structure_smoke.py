import sys  # noqa: I001
from pathlib import Path

# Keep MolOP's RDKit initialization ahead of Open Babel's native extension.
# isort: off
import pytest
from rdkit import Chem
from molgr.interface import xyz_to_rdmol
from openbabel import pybel
# isort: on


sys.path.append(str(Path(__file__).resolve().parent))

from _helpers_structure import (
    build_methane_rdmol,
    build_methane_xyz_block,
    build_water_rdmol,
    build_water_xyz_block,
)


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


def test_molgr_structure_recovery_smoke() -> None:
    recovered = xyz_to_rdmol(build_water_xyz_block())

    assert recovered.GetNumAtoms() == 3
    assert recovered.GetNumBonds() == 2
    assert Chem.MolToSmiles(recovered) == "[H]O[H]"
