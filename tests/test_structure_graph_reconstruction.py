from __future__ import annotations

import sys
from collections.abc import Generator
from pathlib import Path

import pytest
from openbabel import openbabel as ob
from openbabel import pybel

import molop.structure.GraphReconstruction as graph_reconstruction
from molop.structure.GraphReconstruction import (
    MetalAtomPosition,
    assign_charge_radical_for_atom,
    assign_radical_dots,
    break_deformed_ene,
    break_one_bond,
    calc_symmetry_penalty,
    calculate_charge_penalty,
    calculate_coulombic_penalty,
    calculate_metal_penalty,
    calculate_physchem_penalty,
    calculate_radical_penalty,
    calculate_shape_quality,
    calculate_tetrahedron_volume,
    clean_carbene_neighbor_unsaturated,
    clean_neighbor_radicals,
    clean_resonances,
    clean_resonances_0,
    clean_resonances_1,
    clean_resonances_2,
    clean_resonances_3,
    clean_resonances_4,
    clean_resonances_5,
    clean_resonances_6,
    clean_resonances_7,
    clean_resonances_8,
    clean_resonances_9,
    clean_resonances_10,
    clean_resonances_11,
    clean_resonances_12,
    clean_resonances_13,
    combine_metal_with_omol,
    eliminate_1_3_dipole,
    eliminate_carbene_neighbor_heteroatom,
    eliminate_carboxyl,
    eliminate_charge_spliting,
    eliminate_CN_in_doubt,
    eliminate_high_positive_charge_atoms,
    eliminate_negative_charges,
    eliminate_NNN,
    eliminate_positive_charges,
    fresh_omol_charge_radical,
    get_deviation_score,
    get_metal_coordination_sphere,
    get_one_step_resonance,
    get_radical_resonances,
    log_omol_infos,
    make_connections,
    omol_score,
    pre_clean,
    process_resonance,
    validate_omol,
    xyz2omol,
    xyz_to_omol_no_metal,
)


sys.path.append(str(Path(__file__).resolve().parent))

from _helpers_structure import build_methane_xyz_block, build_water_xyz_block


@pytest.fixture(autouse=True)
def clear_xyz_cache() -> Generator[None, None, None]:
    xyz_to_omol_no_metal.cache_clear()
    yield
    xyz_to_omol_no_metal.cache_clear()


def test_xyz2omol_happy_path_water() -> None:
    omol = xyz2omol(build_water_xyz_block(), total_charge=0, total_radical_electrons=0)

    assert omol is not None
    assert omol.OBMol.NumAtoms() == 3
    assert omol.OBMol.NumBonds() == 2


def test_xyz2omol_returns_none_for_negative_radical_request() -> None:
    omol = xyz2omol(build_water_xyz_block(), total_charge=0, total_radical_electrons=-1)

    assert omol is None


def test_xyz_to_omol_no_metal_negative_radical_returns_none() -> None:
    omol = xyz_to_omol_no_metal(build_water_xyz_block(), total_charge=0, total_radical_electrons=-1)

    assert omol is None


def test_validate_omol_true_for_recovered_matching_charge_and_radical() -> None:
    recovered = xyz_to_omol_no_metal(
        build_water_xyz_block(), total_charge=0, total_radical_electrons=0
    )

    assert recovered is not None
    assert validate_omol(recovered, total_charge=0, total_radical_electrons=0)


def test_make_connections_adds_plausible_connectivity_for_tiny_molecule() -> None:
    no_bond_no_xyz = "2\nno\nN 0.000000 0.000000 0.000000\nO 1.900000 0.000000 0.000000\n"
    omol = pybel.readstring("xyz", no_bond_no_xyz)

    assert omol.OBMol.NumAtoms() == 2
    assert omol.OBMol.NumBonds() == 0

    connected = make_connections(omol)

    assert connected.OBMol.NumBonds() >= 1
    assert connected.OBMol.GetBond(1, 2) is not None


def test_graph_reconstruction_scoring_and_penalty_helpers_stable_inputs() -> None:
    recovered_water = xyz2omol(build_water_xyz_block(), total_charge=0, total_radical_electrons=0)
    recovered_methane = xyz2omol(
        build_methane_xyz_block(), total_charge=0, total_radical_electrons=0
    )

    assert recovered_water is not None
    assert recovered_methane is not None

    volume = calculate_tetrahedron_volume(
        (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)
    )
    quality = calculate_shape_quality(
        (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)
    )
    quality_degenerate = calculate_shape_quality(
        (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (2.0, 0.0, 0.0), (3.0, 0.0, 0.0)
    )

    assert volume == pytest.approx(1.0 / 6.0)
    assert 0.0 <= quality <= 1.0
    assert quality_degenerate == 0.0

    oxygen_atom = recovered_water.OBMol.GetAtom(1)
    hydrogen_atom = recovered_water.OBMol.GetAtom(2)
    water_bond = recovered_water.OBMol.GetBond(1, 2)

    assert water_bond is not None
    assert get_deviation_score(recovered_water, 1) >= 0.0
    assert isinstance(calc_symmetry_penalty(recovered_water), float)

    oxygen_atom.SetFormalCharge(-1)
    assert calculate_charge_penalty(oxygen_atom) >= 0.0
    oxygen_atom.SetFormalCharge(0)

    oxygen_atom.SetSpinMultiplicity(1)
    assert calculate_radical_penalty(oxygen_atom) >= 2.0
    oxygen_atom.SetSpinMultiplicity(0)

    methane_carbon = recovered_methane.OBMol.GetAtom(1)
    methane_carbon.SetSpinMultiplicity(1)
    assert calculate_radical_penalty(methane_carbon) >= 0.0
    methane_carbon.SetSpinMultiplicity(0)

    oxygen_atom.SetFormalCharge(1)
    hydrogen_atom.SetFormalCharge(1)
    assert calculate_coulombic_penalty(water_bond) == 15.0
    hydrogen_atom.SetFormalCharge(-1)
    assert calculate_coulombic_penalty(water_bond) == -0.5
    hydrogen_atom.SetFormalCharge(0)
    oxygen_atom.SetFormalCharge(0)
    assert calculate_coulombic_penalty(water_bond) == 0.0

    assert isinstance(calculate_physchem_penalty(recovered_water), float)
    assert calculate_metal_penalty(recovered_water) == 0.0
    assert isinstance(omol_score(recovered_water), float)

    near_atoms = get_metal_coordination_sphere(recovered_water, oxygen_atom, cutoff=1.1)
    assert len(near_atoms) >= 1


def test_graph_reconstruction_cleanup_helpers_are_stable_on_small_molecules() -> None:
    water = xyz2omol(build_water_xyz_block(), total_charge=0, total_radical_electrons=0)
    methane = xyz2omol(build_methane_xyz_block(), total_charge=0, total_radical_electrons=0)

    assert water is not None
    assert methane is not None

    working = pre_clean(water.clone)
    working = clean_carbene_neighbor_unsaturated(working)
    working, charge = eliminate_high_positive_charge_atoms(working, given_charge=0)
    working, charge = eliminate_CN_in_doubt(working, given_charge=charge)
    working, charge = eliminate_carboxyl(working, given_charge=charge)
    working, charge = eliminate_carbene_neighbor_heteroatom(working, given_charge=charge)
    working, charge = eliminate_NNN(working, given_charge=charge)
    working = break_deformed_ene(working, given_charge=charge, given_radical=0)
    working, charge = break_one_bond(working, given_charge=charge, given_radical=0)
    working = clean_neighbor_radicals(working)
    working, charge = eliminate_charge_spliting(working, given_charge=charge)

    resonances = get_one_step_resonance(working)
    all_resonances = get_radical_resonances(working)
    assert isinstance(resonances, list)
    assert isinstance(all_resonances, list)

    working, charge = eliminate_1_3_dipole(working, given_charge=charge)
    working, charge = eliminate_positive_charges(working, given_charge=charge)
    working, charge = eliminate_negative_charges(working, given_charge=charge)
    processed, final_charge = process_resonance(working, charge)

    assert processed.OBMol.NumAtoms() == 3
    assert isinstance(final_charge, int)

    methane_working = methane.clone
    methane_working = clean_resonances_0(methane_working)
    methane_working = clean_resonances_1(methane_working)
    methane_working = clean_resonances_2(methane_working)
    methane_working = clean_resonances_3(methane_working)
    methane_working = clean_resonances_4(methane_working)
    methane_working = clean_resonances_5(methane_working)
    methane_working = clean_resonances_6(methane_working)
    methane_working = clean_resonances_7(methane_working)
    methane_working = clean_resonances_8(methane_working)
    methane_working = clean_resonances_9(methane_working)
    methane_working = clean_resonances_10(methane_working)
    methane_working = clean_resonances_11(methane_working)
    methane_working = clean_resonances_12(methane_working)
    methane_working = clean_resonances_13(methane_working)
    methane_working = clean_resonances(methane_working)

    assert methane_working.OBMol.NumAtoms() == 5


def test_required_logging_assign_validate_and_combine_helpers() -> None:
    water = xyz2omol(build_water_xyz_block(), total_charge=0, total_radical_electrons=0)
    assert water is not None

    log_omol_infos(water, "unit-test")

    refreshed = fresh_omol_charge_radical(water.clone)
    assert refreshed.OBMol.NumAtoms() == water.OBMol.NumAtoms()

    ch_omol = pybel.readstring("xyz", "2\nch\nC 0.0 0.0 0.0\nH 1.0 0.0 0.0\n")
    carbon = ch_omol.OBMol.GetAtom(1)
    carbon.SetFormalCharge(0)
    assert assign_radical_dots(carbon) >= 1
    assign_charge_radical_for_atom(carbon)
    assert carbon.GetSpinMultiplicity() >= 1

    bh4 = pybel.readstring(
        "xyz",
        "5\nbh4\nB 0.0 0.0 0.0\nH 0.7 0.7 0.7\nH -0.7 -0.7 0.7\nH -0.7 0.7 -0.7\nH 0.7 -0.7 -0.7\n",
    )
    boron = bh4.OBMol.GetAtom(1)
    assert boron.GetTotalValence() == 4
    assign_charge_radical_for_atom(boron)
    assert boron.GetFormalCharge() == -1

    oxygen = water.OBMol.GetAtom(1)
    oxygen.SetFormalCharge(1)
    assert validate_omol(water, total_charge=0, total_radical_electrons=0) is False
    oxygen.SetFormalCharge(0)

    oxygen.SetSpinMultiplicity(1)
    assert validate_omol(water, total_charge=0, total_radical_electrons=0) is False
    oxygen.SetSpinMultiplicity(0)

    combined = combine_metal_with_omol(water, metal_list=[])
    assert combined.OBMol.NumAtoms() == water.OBMol.NumAtoms()

    one_metal = [
        MetalAtomPosition(
            idx=1,
            symbol="Li",
            element_idx=3,
            valence=1,
            radical_num=0,
            position_x=5.0,
            position_y=0.0,
            position_z=0.0,
        )
    ]
    with_metal = combine_metal_with_omol(water, metal_list=one_metal)
    assert with_metal.OBMol.NumAtoms() == water.OBMol.NumAtoms() + 1


def _build_test_omol_with_bonds(num_atoms: int = 8) -> pybel.Molecule:
    obmol = ob.OBMol()
    for idx in range(num_atoms):
        atom = obmol.NewAtom()
        atom.SetAtomicNum(6)
        atom.SetVector(float(idx), 0.0, 0.0)
    for idx in range(1, num_atoms):
        obmol.AddBond(idx, idx + 1, 1)
    obmol.AddBond(2, 4, 1)
    return pybel.Molecule(obmol)


def _patch_single_smarts_match(
    monkeypatch: pytest.MonkeyPatch, pattern: str, match: tuple[int, ...]
) -> None:
    class _FakeSmarts:
        def __init__(self, pat: str):
            self.pat = pat
            self.called = False

        def findall(self, _omol: pybel.Molecule):
            if self.pat == pattern and not self.called:
                self.called = True
                return [match]
            return []

    monkeypatch.setattr(graph_reconstruction.pybel, "Smarts", lambda pat: _FakeSmarts(pat))


def test_clean_resonance_branches_with_forced_smarts_matches(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    omol = _build_test_omol_with_bonds()

    def reset_state() -> None:
        for atom in ob.OBMolAtomIter(omol.OBMol):
            atom.SetFormalCharge(0)
            atom.SetSpinMultiplicity(0)
        for bond in ob.OBMolBondIter(omol.OBMol):
            bond.SetBondOrder(1)

    reset_state()
    _patch_single_smarts_match(monkeypatch, "[*-]-[*]=[*]~[*+]", (1, 2, 3, 4))
    omol.OBMol.GetAtom(1).SetFormalCharge(-1)
    omol.OBMol.GetAtom(4).SetFormalCharge(1)
    omol.OBMol.GetBond(2, 3).SetBondOrder(2)
    clean_resonances_0(omol)

    reset_state()
    _patch_single_smarts_match(monkeypatch, "[*-]=[*+]=[*+0]", (1, 2, 3))
    omol.OBMol.GetAtom(1).SetFormalCharge(-1)
    omol.OBMol.GetAtom(2).SetFormalCharge(1)
    omol.OBMol.GetBond(1, 2).SetBondOrder(2)
    omol.OBMol.GetBond(2, 3).SetBondOrder(2)
    clean_resonances_1(omol)

    reset_state()
    _patch_single_smarts_match(
        monkeypatch, "[#8]=[#6](-[!-])-[*]=[*]-[#7-,#6-]", (1, 2, 3, 4, 5, 6)
    )
    omol.OBMol.GetBond(1, 2).SetBondOrder(2)
    omol.OBMol.GetBond(2, 4).SetBondOrder(1)
    omol.OBMol.GetBond(4, 5).SetBondOrder(2)
    omol.OBMol.GetBond(5, 6).SetBondOrder(1)
    omol.OBMol.GetAtom(6).SetFormalCharge(-1)
    clean_resonances_2(omol)

    reset_state()
    _patch_single_smarts_match(monkeypatch, "[#7+]=[*]-[*]=[*]-[#8-]", (1, 2, 3, 4, 5))
    omol.OBMol.GetBond(1, 2).SetBondOrder(2)
    omol.OBMol.GetBond(2, 3).SetBondOrder(1)
    omol.OBMol.GetBond(3, 4).SetBondOrder(2)
    omol.OBMol.GetBond(4, 5).SetBondOrder(1)
    clean_resonances_3(omol)

    reset_state()
    _patch_single_smarts_match(monkeypatch, "[#7+,#8+]=[*]-[#6-,#7-,#8-]", (1, 2, 3))
    omol.OBMol.GetBond(1, 2).SetBondOrder(2)
    clean_resonances_4(omol)

    reset_state()
    _patch_single_smarts_match(monkeypatch, "[#7+0,#8+0,#16+0]=[*+0]-[#6-,#7-]", (1, 2, 3))
    omol.OBMol.GetBond(1, 2).SetBondOrder(2)
    omol.OBMol.GetBond(2, 3).SetBondOrder(1)
    omol.OBMol.GetAtom(3).SetFormalCharge(-1)
    clean_resonances_5(omol)

    reset_state()
    _patch_single_smarts_match(monkeypatch, "[#6]=[#6]=[#6-,#7-]", (1, 2, 3))
    omol.OBMol.GetBond(1, 2).SetBondOrder(2)
    omol.OBMol.GetBond(2, 3).SetBondOrder(2)
    clean_resonances_6(omol)

    reset_state()
    _patch_single_smarts_match(
        monkeypatch, "[*-]1-[*](=[*])-[*]=[*]-[*]=[*]1", (1, 2, 3, 4, 5, 6, 7)
    )
    clean_resonances_7(omol)

    reset_state()
    _patch_single_smarts_match(
        monkeypatch, "[*-]1-[*]=[*]-[*](=[*])-[*]=[*]1", (1, 2, 3, 4, 5, 6, 7)
    )
    clean_resonances_8(omol)

    reset_state()
    _patch_single_smarts_match(monkeypatch, "[*+,*+2,*+3]-,=[*-,*-2,*-3]", (1, 2))
    omol.OBMol.GetAtom(1).SetFormalCharge(1)
    omol.OBMol.GetAtom(2).SetFormalCharge(-1)
    clean_resonances_9(omol)

    reset_state()
    _patch_single_smarts_match(monkeypatch, "[*]-[*]=,#[*]-[*]", (1, 2, 3, 4))
    omol.OBMol.GetAtom(1).SetSpinMultiplicity(1)
    omol.OBMol.GetAtom(4).SetSpinMultiplicity(1)
    omol.OBMol.GetBond(2, 3).SetBondOrder(2)
    clean_resonances_10(omol)

    reset_state()
    _patch_single_smarts_match(monkeypatch, "[#7v3+0,#8v2+0,#16v2+0]-,=,:[*+1]", (1, 2))
    clean_resonances_11(omol)

    reset_state()
    _patch_single_smarts_match(
        monkeypatch, "[#7v3+0,#8v2+0,#16v2+0]-,:[*]=,:[*]-,:[*+1]", (1, 2, 3, 4)
    )
    omol.OBMol.GetBond(1, 2).SetBondOrder(1)
    omol.OBMol.GetBond(2, 3).SetBondOrder(2)
    omol.OBMol.GetBond(3, 4).SetBondOrder(1)
    clean_resonances_12(omol)

    reset_state()
    _patch_single_smarts_match(monkeypatch, "[*-]:[*]=[#7+0,#8+0]", (1, 2, 3))
    omol.OBMol.GetBond(1, 2).SetBondOrder(1)
    clean_resonances_13(omol)

    assert omol.OBMol.NumAtoms() == 8
