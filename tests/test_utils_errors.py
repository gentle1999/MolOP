import pytest

from molop.utils.errors import ArgumentError, MolError, ReactError


def test_argument_error_is_raisable() -> None:
    with pytest.raises(ArgumentError):
        raise ArgumentError("bad argument")


def test_mol_error_is_raisable() -> None:
    with pytest.raises(MolError):
        raise MolError("molecule error")


def test_react_error_is_raisable() -> None:
    with pytest.raises(ReactError):
        raise ReactError("reaction error")
