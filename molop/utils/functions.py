import numpy as np


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


def fill_symmetric_matrix(one_d_array):
    # 计算矩阵的大小
    n = int((-1 + np.sqrt(1 + 8 * len(one_d_array))) // 2)

    # 创建一个n x n的零矩阵
    matrix = np.zeros((n, n), dtype=one_d_array.dtype)

    # 填充下三角矩阵
    for i, j in zip(*np.tril_indices(n)):
        idx = i * (i + 1) // 2 + j
        matrix[i, j] = one_d_array[idx]
        matrix[j, i] = one_d_array[idx]  # 由于矩阵是对称的，填充上三角部分

    return matrix
