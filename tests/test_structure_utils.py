import pytest
from rdkit import Chem

from molop.structure import utils as structure_utils
from molop.structure.utils import canonical_smiles, estimate_bond_length


def test_canonical_smiles_matches_rdkit_reference() -> None:
    smiles = "C(C)"
    expected = Chem.CanonSmiles(smiles, useChiral=True)

    assert canonical_smiles(smiles) == expected


def test_canonical_smiles_is_stable_for_already_canonical_input() -> None:
    smiles = "CC"

    assert canonical_smiles(smiles) == smiles


def test_canonical_smiles_stable_mode_reaches_fixed_point(monkeypatch: pytest.MonkeyPatch) -> None:
    mapping = {"seed": "step-1", "step-1": "fixed", "fixed": "fixed"}
    monkeypatch.setattr(structure_utils, "_canonical_smiles_once", lambda smiles: mapping[smiles])

    assert canonical_smiles("seed", stable=True) == "fixed"
    assert canonical_smiles("seed", stable=False) == "step-1"


def test_canonical_smiles_stable_mode_breaks_cycles(monkeypatch: pytest.MonkeyPatch) -> None:
    mapping = {"seed": "step-1", "step-1": "step-2", "step-2": "step-1"}
    monkeypatch.setattr(structure_utils, "_canonical_smiles_once", lambda smiles: mapping[smiles])

    assert canonical_smiles("seed", stable=True) == "step-2"


def test_canonical_smiles_rejects_invalid_round_count() -> None:
    with pytest.raises(ValueError, match="max_rounds must be >= 1"):
        canonical_smiles("CC", stable=True, max_rounds=0)


def test_estimate_bond_length_orders_common_cc_bonds() -> None:
    single = estimate_bond_length(6, 6, Chem.rdchem.BondType.SINGLE)
    double = estimate_bond_length(6, 6, Chem.rdchem.BondType.DOUBLE)
    triple = estimate_bond_length(6, 6, Chem.rdchem.BondType.TRIPLE)

    assert isinstance(single, float)
    assert isinstance(double, float)
    assert isinstance(triple, float)
    assert single >= double >= triple


def test_estimate_bond_length_raises_for_unsupported_bond_type() -> None:
    with pytest.raises(ValueError, match="Unsupported bond type"):
        estimate_bond_length(6, 6, Chem.rdchem.BondType.AROMATIC)
