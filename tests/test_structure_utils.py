import pytest
from rdkit import Chem

from molop.structure.utils import canonical_smiles, estimate_bond_length


def test_canonical_smiles_matches_rdkit_reference() -> None:
    smiles = "C(C)"
    expected = Chem.CanonSmiles(smiles, useChiral=True)

    assert canonical_smiles(smiles) == expected


def test_canonical_smiles_is_stable_for_already_canonical_input() -> None:
    smiles = "CC"

    assert canonical_smiles(smiles) == smiles


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
