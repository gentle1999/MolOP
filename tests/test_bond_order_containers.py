"""Extensible bond-order records follow the atomic-population container contract."""

import json

import numpy as np
import pytest
from pydantic import ValidationError

from molop.io.base_models.ChemFileFrame import BaseCalcFrame
from molop.io.base_models.data_classes import BondOrderMatrix as ExportedBondOrderMatrix
from molop.io.base_models.DataClasses import BondOrderMatrix, BondOrders
from molop.unit import atom_ureg


def test_arbitrary_bond_order_schemes_keep_matrix_and_provenance_together():
    matrix = np.array([[np.nan, 0.0], [0.0, np.nan]])
    record = BondOrderMatrix(
        scheme="qtaim_delocalization_index",
        values=matrix,
        spin_channel="alpha",
        source_label="Observed pair index",
        print_threshold=0.05,
        missing_reason="pair not printed",
        metadata={"source_option": "custom"},
    )
    orders = BondOrders(bond_orders={"custom_pair_index": record})
    assert ExportedBondOrderMatrix is BondOrderMatrix
    assert orders["custom_pair_index"] is record
    assert orders.get_bond_order("absent") is None
    assert orders.bond_order_names == ["custom_pair_index"]
    assert len(orders) == 1
    assert set(orders.model_dump()) == {"bond_orders"}
    dumped = orders.to_unitless_dump()
    assert dumped["bond_orders"]["custom_pair_index"]["values"] == [[None, 0], [0, None]]
    assert dumped["bond_orders"]["custom_pair_index"]["print_threshold"] == 0.05
    binary = orders.to_unitless_dump_with_unit_keys(array_mode="ndarray")
    assert np.isnan(binary["bond_orders"]["custom_pair_index"]["values"][0, 0])
    tensor = orders.to_property_bundle().tensor_properties["custom_pair_index"]
    assert tensor.metadata["scheme"] == record.scheme
    assert tensor.metadata["spin_channel"] == "alpha"
    assert tensor.metadata["source_label"] == record.source_label
    assert tensor.metadata["missing_reason"] == record.missing_reason
    assert tensor.metadata["source_metadata"] == record.metadata
    json.dumps(dumped, allow_nan=False)
    json.dumps(orders.to_property_bundle().to_unitless_dump_with_unit_keys(), allow_nan=False)
    assert np.isnan(matrix[0, 0])


def test_bond_order_records_validate_atom_alignment_without_a_scheme_allowlist():
    assert len(BondOrders()) == 0
    with pytest.raises(ValidationError, match="non-empty"):
        BondOrderMatrix(scheme="custom", values=np.empty((0, 0)))
    with pytest.raises(ValidationError, match="same shape"):
        BondOrders(
            bond_orders={
                "first": BondOrderMatrix(scheme="first", values=np.eye(2)),
                "second": BondOrderMatrix(scheme="second", values=np.eye(3)),
            }
        )
    with pytest.raises(ValidationError, match="Extra inputs"):
        BondOrders.model_validate({"mayer_bond_order": np.eye(2)})
    with pytest.raises(ValidationError, match="frame atom count"):
        BaseCalcFrame(
            atoms=[1, 1],
            coords=np.zeros((2, 3)) * atom_ureg.angstrom,
            bond_orders=BondOrders(
                bond_orders={
                    "custom": BondOrderMatrix(scheme="custom", values=np.eye(3)),
                }
            ),
        )


@pytest.mark.parametrize("observed", [0.0, 0.75, np.nan])
def test_rdkit_embedding_discovers_custom_schemes_and_preserves_observed_zero(observed):
    frame = BaseCalcFrame(
        qm_software="Custom",
        atoms=[1, 1],
        coords=np.array([[0, 0, 0], [0, 0, 0.74]]) * atom_ureg.angstrom,
        bond_orders=BondOrders(
            bond_orders={
                "custom_pair_index": BondOrderMatrix(
                    scheme="custom",
                    values=np.array([[np.nan, observed], [observed, np.nan]]),
                ),
            }
        ),
    )
    mol = frame.qm_embedded_rdmol()
    assert mol.GetNumBonds() == 1
    bond = mol.GetBondBetweenAtoms(0, 1)
    name = "CUSTOM_PAIR_INDEX_BY_CUSTOM"
    if np.isnan(observed):
        assert not bond.HasProp(name)
    else:
        assert bond.GetDoubleProp(name) == observed


@pytest.mark.parametrize("nonbonded_value", [-0.15, 0.0008])
def test_rdkit_embedding_does_not_turn_pair_indices_into_chemical_bonds(nonbonded_value):
    values = np.array([[0, 0.9, -0.1], [0.9, 0, nonbonded_value], [-0.1, nonbonded_value, 0]])
    frame = BaseCalcFrame(
        qm_software="Custom",
        atoms=[8, 1, 1],
        coords=np.array([[0, 0, 0], [0.96, 0, 0], [-0.24, 0.93, 0]]) * atom_ureg.angstrom,
        bonds=[(0, 1, 1, 0), (0, 2, 1, 0)],
        formal_charges=[0, 0, 0],
        formal_num_radicals=[0, 0, 0],
        bond_orders=BondOrders(
            bond_orders={"custom": BondOrderMatrix(scheme="custom", values=values)}
        ),
    )
    embedded = frame.qm_embedded_rdmol()
    assert embedded.GetNumBonds() == 2
    assert embedded.GetBondBetweenAtoms(1, 2) is None
    assert embedded.GetBondBetweenAtoms(0, 2).GetDoubleProp("CUSTOM_BY_CUSTOM") == -0.1
    assert frame.bond_orders["custom"].values[1, 2] == nonbonded_value
    assert frame.rdmol.GetNumBonds() == 2
    assert not frame.rdmol.GetBondBetweenAtoms(0, 2).HasProp("CUSTOM_BY_CUSTOM")
