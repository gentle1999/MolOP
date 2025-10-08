"""
Author: TMJ
Date: 2025-07-30 16:44:54
LastEditors: TMJ
LastEditTime: 2025-09-25 10:16:26
Description: 请填写简介
"""

from molop.io.base_models.SearchPattern import MolOPPattern


class XYZPatterns:
    META = MolOPPattern(
        start_pattern=r"^\d+\s*\n",
        end_pattern=r"^[A-Za-z]+(\s*-?\d+\.\d+){3}",
    )
    ATOMS = MolOPPattern(
        content_pattern=r"^\s*([A-Z][a-z]*)(\s*-?\d+\.\d+)(\s*-?\d+\.\d+)(\s*-?\d+\.\d+)",
        content_repeat=0,
    )
