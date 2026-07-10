"""
Author: TMJ
Date: 2025-07-30 16:44:54
LastEditors: TMJ
LastEditTime: 2025-10-30 18:43:56
Description: 请填写简介
"""

from molop.io.base_models.SearchPattern import MolOPPattern


class XYZPatterns:
    META = MolOPPattern(
        start_pattern=r"^\d+\s*\n",
        end_pattern=r"^[A-Za-z]+(\s*-?\d+\.\d+){3}",
    )
    ATOMS = MolOPPattern(
        content_pattern=r"^\s*(?P<symbol>[A-Z][a-z]*)(?P<x>\s*-?\d+\.\d+)"
        r"(?P<y>\s*-?\d+\.\d+)(?P<z>\s*-?\d+\.\d+)",
        content_repeat=0,
    )
    CHARGE_MULTIPLICITY = MolOPPattern(
        content_pattern=r"^\s*charge (?P<charge>[+-]?\d+) multiplicity (?P<multiplicity>\d+)",
        description="The charge and multiplicity of the XYZ calculation. link 1",
    )


xyz_patterns = XYZPatterns()
