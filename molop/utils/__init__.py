"""
Author: TMJ
Date: 2023-11-02 15:36:39
LastEditors: TMJ
LastEditTime: 2024-02-21 21:07:05
Description: 请填写简介
"""


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
