import numpy as np
import pytest
from pydantic import ValidationError

from molop.io.base_models.DataClasses import (
    AtomInInternalCoords,
    CoordinateContainer,
    InternalCoords,
)
from molop.unit import atom_ureg


def test_coordinate_container_supports_basic_sequence_operations() -> None:
    container = CoordinateContainer(items=[1, 2])

    container.append(3)
    container.extend([4])

    assert len(container) == 4
    assert container[0] == 1
    assert container[1:3] == [2, 3]


def test_coordinate_container_uses_items_as_the_only_model_field() -> None:
    container = CoordinateContainer(items=[1, 2])

    assert container.items == [1, 2]
    assert container.atoms == [1, 2]
    assert container.model_dump() == {"items": [1, 2]}

    with pytest.raises(ValidationError):
        CoordinateContainer(atoms=[1, 2])


def test_coordinate_container_filters_real_atoms_by_convention() -> None:
    container = CoordinateContainer(
        items=[
            AtomInInternalCoords(symbol="H"),
            AtomInInternalCoords(symbol="X", is_dummy=True),
            AtomInInternalCoords(symbol="He", is_ghost=True),
        ]
    )

    assert container.get_symbols() == ["H"]
    assert [atom.symbol for atom in container.real_atoms()] == ["H"]


def test_internal_coords_still_round_trips_cartesian_coords() -> None:
    coords = np.array([[0.0, 0.0, 0.0], [0.7, 0.0, 0.0], [0.0, 1.0, 0.0]]) * atom_ureg.angstrom
    internal = InternalCoords.from_cartesian_coords(["H", "H", "O"], coords)

    assert len(internal) == 3
    assert internal.get_symbols() == ["H", "H", "O"]
    assert tuple(internal.to_cartesian_coords().shape) == (3, 3)
