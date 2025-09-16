import contextlib
from typing import Generator

import joblib
import numpy as np
from pydantic import BaseModel
from tqdm import tqdm


def is_metal(number: int):
    if number in (
        1,
        2,
        5,
        6,
        7,
        8,
        9,
        10,
        14,
        15,
        16,
        17,
        18,
        33,
        34,
        35,
        36,
        52,
        53,
        54,
        85,
        86,
        118,
    ):
        return False
    else:
        return True


def fill_symmetric_matrix(one_d_array: np.ndarray) -> np.ndarray:
    n = (np.sqrt(1 + 8 * len(one_d_array)) - 1) / 2
    assert int(n) == n, "Array length is not a valid triangular number"
    n = int(n)
    matrix = np.zeros((n, n), dtype=one_d_array.dtype)
    for i, j in zip(*np.tril_indices(n), strict=True):
        idx = i * (i + 1) // 2 + j
        matrix[i, j] = one_d_array[idx]
        matrix[j, i] = one_d_array[idx]
    return matrix


def find_rigid_transform(P: np.ndarray, Q: np.ndarray) -> np.ndarray:
    """
    using Kabsch algorithm to calculate the rigid transform matrix from point set P to point set Q.
    P and Q must be corresponding, i.e., P[i] corresponds to Q[i].

    Parameters:
        P (np.ndarray): the first point set, a numpy array with shape (N, 3).
        Q (np.ndarray): the second point set, a numpy array with shape (N, 3).

    Returns:
        T (np.ndarray): a 4x4 homogeneous transformation matrix that maps points from P to Q.
    """
    # check input dimensions
    assert (
        P.shape == Q.shape
    ), f"Input point sets must have the same dimensions, got {P} and {Q}"
    if P.shape[1] != 3:
        raise ValueError("Input point sets must be Nx3 matrices")
    centroid_P = np.mean(P, axis=0)
    centroid_Q = np.mean(Q, axis=0)
    P_centered = P - centroid_P
    Q_centered = Q - centroid_Q
    H = np.dot(P_centered.T, Q_centered)
    U, S, Vt = np.linalg.svd(H)
    R = np.dot(Vt.T, U.T)
    # If the determinant of R is negative, then we need to invert it,
    # otherwise we are looking in the wrong direction.
    if np.linalg.det(R) < 0:
        Vt[2, :] *= -1  # invert Vt's last row
        R = np.dot(Vt.T, U.T)
    t = centroid_Q.T - np.dot(R, centroid_P.T)
    T = np.identity(4)
    T[:3, :3] = R  # left-top 3x3 is rotation matrix
    T[:3, 3] = t  # right-top 3x1 is translation vector
    return T


def transform_coords(coords: np.ndarray, T: np.ndarray) -> np.ndarray:
    """
    transform the coordinates using the given transformation matrix.

    Parameters:
        coords (np.ndarray): the coordinates to be transformed, a numpy array with shape (N, 3).
        T (np.ndarray): the 4x4 homogeneous transformation matrix.

    Returns:
        transformed_coords (np.ndarray): the transformed coordinates, a numpy array with shape (N, 3).
    """
    return np.dot(coords, T[:3, :3].T) + T[:3, 3]


def invert_transform_coords(
    orientation_coords: np.ndarray, T: np.ndarray
) -> np.ndarray:
    return np.dot(orientation_coords - T[:3, 3], T[:3, :3])


def merge_models(model_1: BaseModel, model_2: BaseModel, force_update: bool = False):
    assert type(model_1) is type(
        model_2
    ), f"Models must be of the same type, got {type(model_1)} and {type(model_2)}"
    field_dict = model_1.model_dump(exclude_unset=True)
    if force_update:
        field_dict.update(model_2.model_dump(exclude_unset=True))
    else:
        for key, value in model_2.model_dump(exclude_unset=True).items():
            if key not in field_dict:
                field_dict[key] = value
    return model_1.__class__.model_validate(field_dict)


@contextlib.contextmanager
def tqdm_joblib(tqdm_object: tqdm) -> Generator[tqdm, None, None]:
    """Context manager to patch joblib to report into tqdm progress bar given as argument"""

    if not isinstance(tqdm_object, tqdm):
        raise ValueError("tqdm_object must be an instance of tqdm")

    old_batch_callback = joblib.parallel.BatchCompletionCallBack

    class TqdmBatchCompletionCallback(joblib.parallel.BatchCompletionCallBack):
        def __call__(self, *args, **kwargs) -> None:
            tqdm_object.update(n=self.batch_size)
            return super().__call__(*args, **kwargs)

    try:
        joblib.parallel.BatchCompletionCallBack = TqdmBatchCompletionCallback
        yield tqdm_object
    except Exception as e:
        print(f"An error occurred: {e}")
    finally:
        joblib.parallel.BatchCompletionCallBack = old_batch_callback
        tqdm_object.close()
