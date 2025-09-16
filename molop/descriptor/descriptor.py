"""
Author: TMJ
Date: 2024-10-19 09:57:26
LastEditors: TMJ
LastEditTime: 2025-09-15 22:45:11
Description: 请填写简介
"""

from typing import List


def calc_rdkit_descs(rdmol, desc_names: List[str]):
    from rdkit.Chem import Descriptors
    from rdkit.ML.Descriptors.MoleculeDescriptors import MolecularDescriptorCalculator

    if desc_names is None:
        desc_names = [desc_name[0] for desc_name in Descriptors._descList]
    desc_calc = MolecularDescriptorCalculator(desc_names)
    return dict(zip(desc_names, desc_calc.CalcDescriptors(rdmol), strict=False))
