import numpy as np
import pytest
from pydantic import BaseModel, ConfigDict, ValidationError

from molop.unit import atom_ureg
from molop.utils.types import (
    Array,
    Array4x4,
    ArrayNx3,
    PintArray,
    PintArray3x3,
    PintArray6Or3x3,
    PintSquareMatrix,
)


class ShapedArrays(BaseModel):
    model_config = ConfigDict(arbitrary_types_allowed=True)

    coords: ArrayNx3[np.float64]
    transform: Array4x4[np.float64]
    tensor: PintArray3x3
    packed_or_tensor: PintArray6Or3x3
    square: PintSquareMatrix


class OpenShapeArrays(BaseModel):
    model_config = ConfigDict(arbitrary_types_allowed=True)

    array: Array[np.float64]
    quantity: PintArray


def test_shaped_numpy_and_pint_types_accept_matching_shapes() -> None:
    payload = ShapedArrays(
        coords=np.zeros((4, 3), dtype=np.float64),
        transform=np.eye(4, dtype=np.float64),
        tensor=np.eye(3) * atom_ureg.ppm,
        packed_or_tensor=np.zeros(6) * atom_ureg.bohr**3,
        square=np.zeros((4, 4)) * atom_ureg.Hz,
    )

    assert payload.coords.shape == (4, 3)
    assert payload.transform.shape == (4, 4)
    assert payload.tensor.shape == (3, 3)
    assert payload.packed_or_tensor.shape == (6,)
    assert payload.square.shape == (4, 4)


def test_open_shape_array_types_require_a_non_scalar_array() -> None:
    payload = OpenShapeArrays(
        array=np.zeros((2, 3, 4), dtype=np.float64),
        quantity=np.zeros((2, 3, 4)) * atom_ureg.bohr,
    )

    assert payload.array.shape == (2, 3, 4)
    assert payload.quantity.shape == (2, 3, 4)

    with pytest.raises(ValidationError, match="at least 1 dimensions"):
        OpenShapeArrays(
            array=np.asarray(1.0),
            quantity=np.asarray(1.0) * atom_ureg.bohr,
        )


@pytest.mark.parametrize(
    ("field", "value"),
    (
        ("coords", np.zeros((4, 2), dtype=np.float64)),
        ("transform", np.zeros((3, 3), dtype=np.float64)),
        ("tensor", np.zeros((2, 2)) * atom_ureg.ppm),
        ("packed_or_tensor", np.zeros(5) * atom_ureg.bohr**3),
        ("square", np.zeros((3, 4)) * atom_ureg.Hz),
    ),
)
def test_shaped_numpy_and_pint_types_reject_mismatched_shapes(field: str, value: object) -> None:
    values = {
        "coords": np.zeros((4, 3), dtype=np.float64),
        "transform": np.eye(4, dtype=np.float64),
        "tensor": np.eye(3) * atom_ureg.ppm,
        "packed_or_tensor": np.zeros(6) * atom_ureg.bohr**3,
        "square": np.zeros((4, 4)) * atom_ureg.Hz,
    }
    values[field] = value

    with pytest.raises(ValidationError, match="must have shape"):
        ShapedArrays.model_validate(values)
