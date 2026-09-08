from collections import Counter
from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem, BondType, rdMolDescriptors

from molop.io.logic.gaussian.log.parsers.G16LogFileParser import G16LogFileParserMemory
from molop.structure import replace_multisite_substituent, replace_substituent
from molop.structure.substitution_sites import (
    find_multisite_sites,
    find_sites,
    normalize_query,
)


FIXTURE = Path(__file__).resolve().parent / "test_files" / "g16log" / "2-TS1-Sp.log"
GENERATED_ATOM_PROP = "_molop_generated_substituent"


@pytest.fixture(scope="module")
def ts_parent() -> Chem.Mol:
    parsed = G16LogFileParserMemory().parse(FIXTURE.read_text(encoding="utf-8"))
    frame = parsed.frames[-1]
    molecule = frame.rdmol
    assert molecule is not None
    return Chem.Mol(molecule)


def _element_counts(mol: Chem.Mol) -> Counter[str]:
    return Counter(atom.GetSymbol() for atom in mol.GetAtoms())


def _sites(mol: Chem.Mol, query: str):
    return find_sites(mol, normalize_query(query))


def _assert_sanitized_without_internal_marks(mol: Chem.Mol) -> None:
    Chem.SanitizeMol(mol)
    assert not any(atom.HasProp(GENERATED_ATOM_PROP) for atom in mol.GetAtoms())


def test_ts_fixture_is_a_complex_coordinate_bearing_parent(ts_parent: Chem.Mol) -> None:
    assert ts_parent.GetNumAtoms() == 101
    assert ts_parent.GetNumConformers() == 1
    assert rdMolDescriptors.CalcMolFormula(ts_parent) == "C44H45F3NNiO4PS2"
    assert {atom.GetSymbol() for atom in ts_parent.GetAtoms()} >= {
        "C",
        "F",
        "Ni",
        "N",
        "O",
        "P",
        "S",
    }


def test_ts_parent_exposes_real_pendant_fragments(ts_parent: Chem.Mol) -> None:
    cf3_sites = _sites(ts_parent, "C(F)(F)F")
    tert_butyl_sites = _sites(ts_parent, "C(C)(C)C")
    phenyl_sites = _sites(ts_parent, "c1ccccc1")

    assert len(cf3_sites) == 1
    assert ts_parent.GetAtomWithIdx(cf3_sites[0].scaffold_atom).GetSymbol() == "S"
    assert ts_parent.GetAtomWithIdx(cf3_sites[0].substituent_atom).GetSymbol() == "C"
    assert cf3_sites[0].bond_type == BondType.SINGLE

    assert len(tert_butyl_sites) == 1
    assert ts_parent.GetAtomWithIdx(tert_butyl_sites[0].scaffold_atom).GetSymbol() == "S"
    assert ts_parent.GetAtomWithIdx(tert_butyl_sites[0].substituent_atom).GetSymbol() == "C"

    # Two phenyl groups are attached to phosphorus and one to the chelating
    # carbon. The fused aromatic ligand is correctly rejected as multi-attached.
    assert len(phenyl_sites) == 3
    assert sorted(
        ts_parent.GetAtomWithIdx(site.scaffold_atom).GetSymbol() for site in phenyl_sites
    ) == ["C", "P", "P"]


def test_replace_cf3_with_complex_oxygenated_substituent(ts_parent: Chem.Mol) -> None:
    source = Chem.Mol(ts_parent)
    before = _element_counts(source)

    replaced = replace_substituent(
        source,
        query="C(F)(F)F",
        replacement="[*]COC(=O)C",
        attempt_num=2,
        angle_split=8,
    )

    after = _element_counts(replaced)
    assert after["F"] == before["F"] - 3
    assert after["O"] == before["O"] + 2
    assert after["Ni"] == before["Ni"] == 1
    assert replaced.GetNumConformers() == 1
    _assert_sanitized_without_internal_marks(replaced)


def test_replace_substituent_retains_parent_3d_conformer_geometry(
    ts_parent: Chem.Mol,
) -> None:
    source = Chem.Mol(ts_parent)
    for atom in source.GetAtoms():
        atom.SetIntProp("_test_parent_atom_idx", atom.GetIdx())

    replaced = replace_substituent(
        source,
        query="C(F)(F)F",
        replacement="[*]COC(=O)C",
        attempt_num=2,
        angle_split=8,
    )

    assert source.GetNumConformers() == 1
    assert replaced.GetNumConformers() == 1
    assert replaced.GetConformer().Is3D()

    retained = [
        (atom.GetIntProp("_test_parent_atom_idx"), atom.GetIdx())
        for atom in replaced.GetAtoms()
        if atom.HasProp("_test_parent_atom_idx")
    ]
    source_indices = [source_idx for source_idx, _result_idx in retained]
    result_indices = [result_idx for _source_idx, result_idx in retained]
    source_positions = np.asarray(source.GetConformer().GetPositions())[source_indices]
    result_positions = np.asarray(replaced.GetConformer().GetPositions())[result_indices]

    source_distances = np.linalg.norm(
        source_positions[:, None, :] - source_positions[None, :, :], axis=2
    )
    result_distances = np.linalg.norm(
        result_positions[:, None, :] - result_positions[None, :, :], axis=2
    )
    assert np.max(np.abs(source_positions - result_positions)) == pytest.approx(0.0, abs=1e-6)
    assert np.max(np.abs(source_distances - result_distances)) == pytest.approx(0.0, abs=1e-6)


def test_replace_tert_butyl_with_complex_oxygenated_substituent(
    ts_parent: Chem.Mol,
) -> None:
    source = Chem.Mol(ts_parent)
    before = _element_counts(source)

    replaced = replace_substituent(
        source,
        query="C(C)(C)C",
        replacement="[*]COC",
        attempt_num=2,
        angle_split=8,
    )

    after = _element_counts(replaced)
    assert after["C"] < before["C"]
    assert after["O"] == before["O"] + 1
    assert not _sites(replaced, "C(C)(C)C")
    _assert_sanitized_without_internal_marks(replaced)


def test_replace_all_real_phenyl_ligands_with_nitrile_substituents(
    ts_parent: Chem.Mol,
) -> None:
    source = Chem.Mol(ts_parent)
    query = normalize_query("c1ccccc1")
    assert len(find_sites(source, query)) == 3

    replaced = replace_substituent(
        source,
        query="c1ccccc1",
        replacement="[*]CC#N",
        replace_all=True,
        attempt_num=2,
        angle_split=8,
    )

    assert not find_sites(replaced, query)
    nitrile_bonds = [
        bond
        for bond in replaced.GetBonds()
        if bond.GetBondType() == BondType.TRIPLE
        and {bond.GetBeginAtom().GetSymbol(), bond.GetEndAtom().GetSymbol()} == {"C", "N"}
    ]
    assert len(nitrile_bonds) == 3
    _assert_sanitized_without_internal_marks(replaced)


def test_replace_explicit_nickel_aryl_dative_site(ts_parent: Chem.Mol) -> None:
    replaced = replace_substituent(
        Chem.Mol(ts_parent),
        query=None,
        replacement="[*]N",
        start_idx=0,
        end_idx=1,
        attempt_num=2,
        angle_split=8,
    )

    nickel = next(atom for atom in replaced.GetAtoms() if atom.GetSymbol() == "Ni")
    dative_bonds = [
        bond
        for bond in nickel.GetBonds()
        if bond.GetBondType()
        in {BondType.DATIVE, BondType.DATIVEL, BondType.DATIVER, BondType.DATIVEONE}
    ]
    assert len(dative_bonds) == 4
    assert sum(neighbor.GetSymbol() == "N" for neighbor in nickel.GetNeighbors()) == 1
    _assert_sanitized_without_internal_marks(replaced)


def test_experimental_multisite_replacement_freezes_fused_epoxide_parent() -> None:
    source = Chem.AddHs(Chem.MolFromSmiles("C1CC2OC2C1"))
    assert source is not None
    assert AllChem.EmbedMolecule(source, randomSeed=7) == 0
    for atom in source.GetAtoms():
        atom.SetIntProp("_test_parent_atom_idx", atom.GetIdx())
    source_positions = np.asarray(source.GetConformer().GetPositions())

    sites = find_multisite_sites(source, normalize_query("C1CO1"), 2)
    assert len(sites) == 1
    assert {(anchor.scaffold_atom, anchor.substituent_atom) for anchor in sites[0].anchors} == {
        (1, 2),
        (5, 4),
    }

    replaced = replace_multisite_substituent(
        source,
        query="C1CO1",
        replacement="[*]C1C([*])C2CCC2C1",
        attempt_num=2,
        angle_split=36,
    )

    retained = [
        (atom.GetIntProp("_test_parent_atom_idx"), atom.GetIdx())
        for atom in replaced.GetAtoms()
        if atom.HasProp("_test_parent_atom_idx")
    ]
    assert len(retained) == 9
    result_positions = np.asarray(replaced.GetConformer().GetPositions())
    np.testing.assert_allclose(
        result_positions[[result_idx for _source_idx, result_idx in retained]],
        source_positions[[source_idx for source_idx, _result_idx in retained]],
        atol=1e-6,
    )
    _assert_sanitized_without_internal_marks(replaced)


def test_experimental_multisite_replacement_works_on_real_ts_fixture(
    ts_parent: Chem.Mol,
) -> None:
    source = Chem.Mol(ts_parent)
    for atom in source.GetAtoms():
        atom.SetIntProp("_test_parent_atom_idx", atom.GetIdx())
    source_positions = np.asarray(source.GetConformer().GetPositions())

    replaced = replace_multisite_substituent(
        source,
        query="c1ccccc1",
        replacement="[*]C1CC([*])CCC1",
        attempt_num=2,
        crowding_threshold=0.6,
        angle_split=24,
    )

    retained = [
        (atom.GetIntProp("_test_parent_atom_idx"), atom.GetIdx())
        for atom in replaced.GetAtoms()
        if atom.HasProp("_test_parent_atom_idx")
    ]
    assert len(retained) == 91
    assert replaced.GetNumAtoms() == 107
    result_positions = np.asarray(replaced.GetConformer().GetPositions())
    np.testing.assert_allclose(
        result_positions[[result_idx for _source_idx, result_idx in retained]],
        source_positions[[source_idx for source_idx, _result_idx in retained]],
        atol=1e-6,
    )
    _assert_sanitized_without_internal_marks(replaced)


def test_experimental_multisite_replacement_accepts_disconnected_query_fragments() -> None:
    source = Chem.AddHs(Chem.MolFromSmiles("OCCO"))
    assert source is not None
    assert AllChem.EmbedMolecule(source, randomSeed=31) == 0
    for atom in source.GetAtoms():
        atom.SetIntProp("_test_parent_atom_idx", atom.GetIdx())
    source_positions = np.asarray(source.GetConformer().GetPositions())

    sites = find_multisite_sites(
        source,
        normalize_query("O.O", require_connected=False),
        2,
    )
    assert len(sites) == 1
    assert len(sites[0].anchors) == 2
    assert len(sites[0].query_match) == 2

    replaced = replace_multisite_substituent(
        source,
        query="O.O",
        replacement="[*]CCCC[*]",
        attempt_num=2,
        crowding_threshold=0.1,
        anchor_tolerance=0.5,
    )

    retained = [
        (atom.GetIntProp("_test_parent_atom_idx"), atom.GetIdx())
        for atom in replaced.GetAtoms()
        if atom.HasProp("_test_parent_atom_idx")
    ]
    result_positions = np.asarray(replaced.GetConformer().GetPositions())
    np.testing.assert_allclose(
        result_positions[[result_idx for _source_idx, result_idx in retained]],
        source_positions[[source_idx for source_idx, _result_idx in retained]],
        atol=1e-6,
    )
    assert len(Chem.GetMolFrags(replaced)) == 1
    _assert_sanitized_without_internal_marks(replaced)


def test_experimental_three_anchor_replacement_freezes_parent_coordinates() -> None:
    source = Chem.AddHs(Chem.MolFromSmiles("C1(C)C(C)C1(C)"))
    assert source is not None
    assert AllChem.EmbedMolecule(source, randomSeed=17) == 0
    for atom in source.GetAtoms():
        atom.SetIntProp("_test_parent_atom_idx", atom.GetIdx())
    source_positions = np.asarray(source.GetConformer().GetPositions())

    sites = find_multisite_sites(source, normalize_query("C1CC1"), 3)
    assert len(sites) == 1
    assert len(sites[0].anchors) == 3

    replaced = replace_multisite_substituent(
        source,
        query="C1CC1",
        replacement="[*]C1C([*])C1([*])C2CC2",
        attempt_num=2,
        crowding_threshold=0.3,
    )

    retained = [
        (atom.GetIntProp("_test_parent_atom_idx"), atom.GetIdx())
        for atom in replaced.GetAtoms()
        if atom.HasProp("_test_parent_atom_idx")
    ]
    assert len(retained) == 12
    result_positions = np.asarray(replaced.GetConformer().GetPositions())
    np.testing.assert_allclose(
        result_positions[[result_idx for _source_idx, result_idx in retained]],
        source_positions[[source_idx for source_idx, _result_idx in retained]],
        atol=1e-6,
    )
    _assert_sanitized_without_internal_marks(replaced)


def test_experimental_four_anchor_replacement_freezes_parent_coordinates() -> None:
    source = Chem.AddHs(Chem.MolFromSmiles("C1(C)C(C)C(C)C1(C)"))
    assert source is not None
    assert AllChem.EmbedMolecule(source, randomSeed=23) == 0
    for atom in source.GetAtoms():
        atom.SetIntProp("_test_parent_atom_idx", atom.GetIdx())
    source_positions = np.asarray(source.GetConformer().GetPositions())

    sites = find_multisite_sites(source, normalize_query("C1CCC1"), 4)
    assert len(sites) == 1
    assert len(sites[0].anchors) == 4

    replaced = replace_multisite_substituent(
        source,
        query="C1CCC1",
        replacement="[*]C1C([*])C([*])C1([*])CCCC",
        attempt_num=2,
    )

    retained = [
        (atom.GetIntProp("_test_parent_atom_idx"), atom.GetIdx())
        for atom in replaced.GetAtoms()
        if atom.HasProp("_test_parent_atom_idx")
    ]
    assert len(retained) == 16
    result_positions = np.asarray(replaced.GetConformer().GetPositions())
    np.testing.assert_allclose(
        result_positions[[result_idx for _source_idx, result_idx in retained]],
        source_positions[[source_idx for source_idx, _result_idx in retained]],
        atol=1e-6,
    )
    _assert_sanitized_without_internal_marks(replaced)
