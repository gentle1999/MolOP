import sys
from pathlib import Path

import pytest
from rdkit import Chem

from molop.structure import FormatConverter as fc


sys.path.append(str(Path(__file__).resolve().parent))

from _helpers_structure import build_methane_rdmol, build_water_rdmol


def test_rdmol_to_omol_roundtrip_atom_count() -> None:
    rdmol = build_water_rdmol()

    omol = fc.rdmol_to_omol(rdmol)

    assert omol.OBMol.NumAtoms() == rdmol.GetNumAtoms()


def test_omol_to_rdmol_by_graph_returns_sanitizable_mol() -> None:
    source = build_methane_rdmol()
    omol = fc.rdmol_to_omol(source)

    converted = fc.omol_to_rdmol_by_graph(omol)

    assert converted is not None
    assert converted.GetNumAtoms() == source.GetNumAtoms()
    assert converted.GetNumBonds() == source.GetNumBonds()
    Chem.SanitizeMol(converted)


def test_validate_rdmol_accepts_expected_charge_and_radical() -> None:
    methane = build_methane_rdmol()

    assert fc.validate_rdmol(methane, total_charge=0, total_radical=0)
    assert not fc.validate_rdmol(methane, total_charge=1, total_radical=0)


def test_validate_rdmol_singlet_normalization_rule() -> None:
    singlet_carbene = Chem.MolFromSmiles("[CH2]")

    assert singlet_carbene is not None
    assert sum(atom.GetNumRadicalElectrons() for atom in singlet_carbene.GetAtoms()) == 2
    assert fc.validate_rdmol(singlet_carbene, total_charge=0, total_radical=0)
    assert not fc.validate_rdmol(singlet_carbene, total_charge=0, total_radical=1)


def test_omol_to_rdmol_uses_molblock_path_when_valid(monkeypatch: pytest.MonkeyPatch) -> None:
    source = build_water_rdmol()
    omol = fc.rdmol_to_omol(source)

    def fail_if_called(_omol: object) -> None:
        raise AssertionError("graph fallback should not be called for valid molblock path")

    monkeypatch.setattr(fc, "omol_to_rdmol_by_graph", fail_if_called)

    converted = fc.omol_to_rdmol(omol, total_charge=0, total_radical=0)

    assert converted is not None
    assert converted.GetNumAtoms() == source.GetNumAtoms()
    assert converted.GetNumBonds() == source.GetNumBonds()


def test_omol_to_rdmol_falls_back_to_graph_when_molblock_invalid(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    source = build_methane_rdmol()
    omol = fc.rdmol_to_omol(source)

    real_mol_from_molblock = fc.Chem.MolFromMolBlock
    call_count = {"n": 0}

    def first_call_invalid(
        molBlock: str,
        sanitize: bool = True,
        removeHs: bool = True,
        strictParsing: bool = True,
    ):
        call_count["n"] += 1
        if call_count["n"] == 1:
            return None
        return real_mol_from_molblock(
            molBlock,
            sanitize=sanitize,
            removeHs=removeHs,
            strictParsing=strictParsing,
        )

    monkeypatch.setattr(fc.Chem, "MolFromMolBlock", first_call_invalid)

    converted = fc.omol_to_rdmol(omol, total_charge=0, total_radical=0)

    assert call_count["n"] >= 2
    assert converted is not None
    assert converted.GetNumAtoms() == source.GetNumAtoms()
    assert converted.GetNumBonds() == source.GetNumBonds()
    Chem.SanitizeMol(converted)


def test_rdmol_to_gjf_geom_has_expected_block_shape() -> None:
    methane = build_methane_rdmol()

    block = fc.rdmol_to_gjf_geom(methane, total_radical=0)

    lines = block.splitlines()
    assert len(lines) == methane.GetNumAtoms() + 1

    charge_text, multiplicity_text = lines[0].split()
    assert int(charge_text) == 0
    assert int(multiplicity_text) == 1

    for idx, line in enumerate(lines[1:]):
        tokens = line.split()
        assert len(tokens) == 4
        assert tokens[0] == methane.GetAtomWithIdx(idx).GetSymbol()
        float(tokens[1])
        float(tokens[2])
        float(tokens[3])


def test_rdmol_to_gjf_connectivity_uses_one_based_indices() -> None:
    methane = build_methane_rdmol()

    block = fc.rdmol_to_gjf_connectivity(methane)

    lines = block.splitlines()
    assert len(lines) == methane.GetNumAtoms()
    atom_count = methane.GetNumAtoms()

    for expected_idx, line in enumerate(lines, start=1):
        tokens = line.split()
        assert tokens[0] == str(expected_idx)
        assert (len(tokens) - 1) % 2 == 0
        for offset in range(1, len(tokens), 2):
            neighbor_idx = int(tokens[offset])
            bond_order = float(tokens[offset + 1])
            assert 1 <= neighbor_idx <= atom_count
            assert bond_order > 0
