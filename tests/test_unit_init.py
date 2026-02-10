import pytest
from pint.errors import UndefinedUnitError

from molop.unit import atom_ureg, si_ureg, unit_transform


def test_unit_transform_identity_returns_same_quantity_instance():
    value = si_ureg.Quantity(2.5, "meter")

    transformed = unit_transform(value, "meter")

    assert transformed is value
    assert transformed is not None
    assert transformed.magnitude == pytest.approx(2.5)


def test_unit_transform_handles_none_input():
    assert unit_transform(None, "meter") is None


def test_unit_transform_converts_with_numeric_tolerance():
    one_angstrom = atom_ureg.Quantity(1.0, "angstrom")

    converted = unit_transform(one_angstrom, "bohr")
    round_trip = unit_transform(converted, "angstrom")

    assert converted is not None
    assert converted.magnitude == pytest.approx(1.8897261, rel=1e-6)
    assert str(converted.units) == "bohr"
    assert round_trip is not None
    assert round_trip.magnitude == pytest.approx(1.0, rel=1e-12)


def test_unit_transform_raises_meaningful_error_for_invalid_unit():
    value = si_ureg.Quantity(1.0, "meter")

    with pytest.raises(UndefinedUnitError, match="not_a_real_unit"):
        unit_transform(value, "not_a_real_unit")
