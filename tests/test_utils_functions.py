import numpy as np
import pytest
from pydantic import BaseModel

from molop.utils.functions import (
    fill_symmetric_matrix,
    find_rigid_transform,
    invert_transform_coords,
    is_metal,
    merge_models,
    transform_coords,
    verify_transform_reversibility,
)


def test_is_metal():
    assert is_metal(3)  # Li is metal (atomic number 3)
    assert not is_metal(1)  # H is not metal
    assert not is_metal(6)  # C is not metal
    assert is_metal(26)  # Fe is metal


def test_fill_symmetric_matrix():
    # 1D array representing lower triangular part of 3x3 matrix
    # indices: (0,0), (1,0), (1,1), (2,0), (2,1), (2,2)
    one_d = np.array([1, 2, 3, 4, 5, 6])
    expected = np.array([[1, 2, 4], [2, 3, 5], [4, 5, 6]])
    result = fill_symmetric_matrix(one_d)
    np.testing.assert_array_equal(result, expected)

    with pytest.raises(AssertionError):
        fill_symmetric_matrix(np.array([1, 2, 3, 4]))  # Invalid length


def test_rigid_transform():
    # Create a random set of points
    P = np.random.rand(10, 3)

    # Define a known rotation and translation
    theta = np.radians(30)
    c, s = np.cos(theta), np.sin(theta)
    R_true = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
    t_true = np.array([1.0, 2.0, 3.0])

    # Transform P to get Q
    Q = (R_true @ P.T).T + t_true

    # Find transform
    T_calc = find_rigid_transform(P, Q)

    # Check rotation part
    np.testing.assert_allclose(T_calc[:3, :3], R_true, atol=1e-5)
    # Check translation part
    np.testing.assert_allclose(T_calc[:3, 3], t_true, atol=1e-5)


def test_verify_transform_reversibility():
    P = np.random.rand(10, 3)
    # Trivial case: P -> P (identity transform)
    assert verify_transform_reversibility(P, P)

    # Case with shift
    Q = P + np.array([1, 0, 0])
    assert verify_transform_reversibility(P, Q)


def test_transform_coords():
    coords = np.array([[1.0, 0.0, 0.0]])
    T = np.eye(4)
    T[:3, 3] = np.array([1.0, 2.0, 3.0])  # Translate

    transformed = transform_coords(coords, T)
    expected = np.array([[2.0, 2.0, 3.0]])
    np.testing.assert_allclose(transformed, expected)


def test_invert_transform_coords():
    coords = np.array([[2.0, 2.0, 3.0]])
    T = np.eye(4)
    T[:3, 3] = np.array([1.0, 2.0, 3.0])

    inverted = invert_transform_coords(coords, T)
    expected = np.array([[1.0, 0.0, 0.0]])
    np.testing.assert_allclose(inverted, expected)


class MockModel(BaseModel):
    a: int = 1
    b: str = "test"


def test_merge_models():
    m1 = MockModel(a=1, b="original")
    m2 = MockModel(a=2, b="new")

    # Default: do not force update (m1 keeps its values if set, but here they are set)
    # Actually merge_models logic:
    # for key, value in model_2...: if key not in field_dict: field_dict[key] = value
    # field_dict comes from m1.
    # So if m1 has 'a', it keeps 'a'.

    merged = merge_models(m1, m2, force_update=False)
    assert merged.a == 1
    assert merged.b == "original"

    merged_force = merge_models(m1, m2, force_update=True)
    assert merged_force.a == 2
    assert merged_force.b == "new"
