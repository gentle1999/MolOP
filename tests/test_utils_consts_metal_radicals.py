import pytest

from molop.utils.consts import get_possible_metal_radicals


def test_get_possible_metal_radicals_fe2_is_set_and_contains_zero() -> None:
    radicals = get_possible_metal_radicals("Fe", 2)
    assert isinstance(radicals, set)
    assert 0 in radicals


def test_get_possible_metal_radicals_raises_when_valence_too_high() -> None:
    with pytest.raises(ValueError):
        get_possible_metal_radicals("Fe", 999)
