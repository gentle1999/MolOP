from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from molop.io.base_models.Bases import BaseDataClassWithUnit, PropertyBundle
from molop.io.base_models.DataClasses import EnergyObservation, ThermalInformations
from molop.io.logic.gaussian.log.parsers.G16LogFileParser import G16LogFileParserMemory
from molop.unit import atom_ureg


def test_unitless_dump_with_unit_keys_uses_normalized_model_units() -> None:
    thermal = ThermalInformations(
        ZPVE=1 * atom_ureg.hartree / atom_ureg.particle,
        moments_of_inertia=np.array([1.0, 2.0, 3.0]) * atom_ureg.amu * atom_ureg.bohr**2,
    )

    assert thermal.to_unitless_dump_with_unit_keys(exclude_none=True) == {
        "ZPVE (kilocalorie/mole)": thermal.ZPVE.magnitude,
        "moments_of_inertia (bohr^2*unified_atomic_mass_unit)": [1.0, 2.0, 3.0],
    }


def test_unitless_dump_unit_labels_do_not_use_pint_display_strings(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    model = _QuantityList(values=[1.0 * atom_ureg.amu * atom_ureg.bohr**2])

    def _unexpected_unit_string(_: object) -> str:
        raise AssertionError("Pint display formatting must not define the public unit label")

    monkeypatch.setattr(type(atom_ureg.hartree), "__str__", _unexpected_unit_string)

    assert model.to_unitless_dump_with_unit_keys() == {
        "values (bohr^2*unified_atomic_mass_unit)": [1.0]
    }


def test_unitless_dump_with_unit_keys_recurses_into_models_and_mappings() -> None:
    bundle = PropertyBundle(
        scalar_properties={"energy": -1.0 * atom_ureg.hartree, "count": 2},
        metadata={
            "nested": ThermalInformations(S=3 * atom_ureg.calorie / atom_ureg.mol / atom_ureg.K)
        },
    )

    assert bundle.to_unitless_dump_with_unit_keys(exclude_none=True) == {
        "scalar_properties": {"energy (hartree)": -1.0, "count": 2},
        "vector_properties": {},
        "tensor_properties": {},
        "tables": {},
        "spectra": {},
        "series": {},
        "transitions": {},
        "metadata": {"nested": {"S (calorie/(kelvin*mole))": 3}},
    }


class _QuantityList(BaseDataClassWithUnit):
    values: list


def test_unitless_dump_with_unit_keys_annotates_homogeneous_quantity_lists() -> None:
    model = _QuantityList(
        values=[
            np.array([1.0, 2.0]) * atom_ureg.angstrom,
            np.array([3.0, 4.0]) * atom_ureg.angstrom,
        ]
    )

    assert model.to_unitless_dump_with_unit_keys() == {
        "values (angstrom)": [[1.0, 2.0], [3.0, 4.0]]
    }


def test_unitless_dump_with_unit_keys_handles_nested_observation_models() -> None:
    observation = EnergyObservation(
        method="CCSD(T)",
        quantity_semantics="total_energy",
        value=-10.5 * atom_ureg.hartree,
        source_label="E(CCSD(T))",
    )

    assert observation.to_unitless_dump_with_unit_keys() == {
        "method": "CCSD(T)",
        "quantity_semantics": "total_energy",
        "value (hartree)": -10.5,
        "source_label": "E(CCSD(T))",
    }


class _ArrayPayload(BaseDataClassWithUnit):
    quantity_array: Any
    plain_array: np.ndarray
    label_array: np.ndarray


def test_unitless_dump_with_unit_keys_defaults_to_json_safe_array_lists() -> None:
    model = _ArrayPayload(
        quantity_array=np.array([1.0, 2.0]) * atom_ureg.hartree,
        plain_array=np.array([[1, 2], [3, 4]], dtype=np.int64),
        label_array=np.array(["minimum", "ts"]),
    )

    payload = model.to_unitless_dump_with_unit_keys()

    assert payload == {
        "quantity_array (hartree)": [1.0, 2.0],
        "plain_array": [[1, 2], [3, 4]],
        "label_array": ["minimum", "ts"],
    }
    json.dumps(payload)


def test_unitless_dump_with_unit_keys_can_preserve_numeric_ndarray_copies() -> None:
    model = _ArrayPayload(
        quantity_array=np.array([1.0, 2.0]) * atom_ureg.hartree,
        plain_array=np.array([[1, 2], [3, 4]], dtype=np.int64),
        label_array=np.array(["minimum", "ts"]),
    )

    payload = model.to_unitless_dump_with_unit_keys(array_mode="ndarray")

    quantity_array = payload["quantity_array (hartree)"]
    plain_array = payload["plain_array"]
    assert isinstance(quantity_array, np.ndarray)
    assert isinstance(plain_array, np.ndarray)
    assert not np.shares_memory(quantity_array, model.quantity_array.magnitude)
    assert not np.shares_memory(plain_array, model.plain_array)
    assert payload["label_array"] == ["minimum", "ts"]


@pytest.mark.parametrize(
    "value",
    [float("nan"), float("inf"), np.array([1.0, np.nan])],
)
def test_json_array_mode_rejects_non_finite_values(value: Any) -> None:
    model = _NestedSelection(energy=value, ignored=0)

    with pytest.raises(ValueError, match="non-finite"):
        model.to_unitless_dump_with_unit_keys()


def test_complex_arrays_require_the_binary_sidecar_mode() -> None:
    model = _ArrayPayload(
        quantity_array=np.array([1.0]),
        plain_array=np.array([1.0 + 2.0j]),
        label_array=np.array(["complex"]),
    )

    with pytest.raises(ValueError, match="complex arrays"):
        model.to_unitless_dump_with_unit_keys()
    assert np.array_equal(
        model.to_unitless_dump_with_unit_keys(array_mode="ndarray")["plain_array"],
        model.plain_array,
    )


def test_real_frame_dump_is_separate_complete_and_repeatable() -> None:
    fixture = Path(__file__).resolve().parent / "test_files" / "g16log" / "H2O.log"
    chem_file = G16LogFileParserMemory().parse(fixture.read_text())
    file_payload = chem_file.to_unitless_dump_with_unit_keys(exclude_none=True)
    frame = chem_file.frames[0]
    source_atoms = list(frame.atoms)
    source_coords = frame.coords.magnitude.copy()

    assert "frames" not in file_payload
    assert "_frames_" not in file_payload
    assert frame.bonds == []
    assert frame.formal_charges == []
    assert frame.formal_num_radicals == []

    first = frame.to_unitless_dump_with_unit_keys(exclude_none=True)
    second = frame.to_unitless_dump_with_unit_keys(exclude_none=True)

    assert first == second
    json.dumps(file_payload, allow_nan=False)
    json.dumps(first, allow_nan=False)
    assert first["atoms"] == source_atoms == [8, 1, 1]
    np.testing.assert_allclose(first["coords (angstrom)"], source_coords)
    assert first["bonds"] == [(0, 1, 1, 0), (0, 2, 1, 0)]
    assert first["formal_charges"] == [0, 0, 0]
    assert first["formal_num_radicals"] == [0, 0, 0]

    ndarray_payload = frame.to_unitless_dump_with_unit_keys(
        exclude_none=True,
        array_mode="ndarray",
    )
    exported_coords = ndarray_payload["coords (angstrom)"]
    assert isinstance(exported_coords, np.ndarray)
    assert not np.shares_memory(exported_coords, frame.coords.magnitude)


class _NestedSelection(BaseDataClassWithUnit):
    energy: object
    ignored: int


class _SelectionPayload(BaseDataClassWithUnit):
    kept: int
    omitted: None = None
    nested: _NestedSelection


def test_unitless_dump_with_unit_keys_preserves_model_dump_selection() -> None:
    model = _SelectionPayload(
        kept=1,
        nested=_NestedSelection(energy=-2.0 * atom_ureg.hartree, ignored=3),
    )

    payload = model.to_unitless_dump_with_unit_keys(
        include={"kept": True, "omitted": True, "nested": {"energy", "ignored"}},
        exclude={"nested": {"ignored"}},
        exclude_none=True,
    )

    assert payload == {"kept": 1, "nested": {"energy (hartree)": -2.0}}


def test_unitless_dump_with_unit_keys_rejects_unknown_array_mode() -> None:
    model = _QuantityList(values=[])

    with pytest.raises(ValueError, match="array_mode"):
        model.to_unitless_dump_with_unit_keys(array_mode="sidecar")  # type: ignore[arg-type]
