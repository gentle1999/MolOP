"""
Author: TMJ
Date: 2023-10-30 14:05:04
LastEditors: TMJ
LastEditTime: 2024-03-02 13:22:19
Description: 请填写简介
"""

from typing import TypeVar

import numpy as np
import numpy.typing as npt
from rdkit import Chem

RdMol = Chem.rdchem.Mol
RWMol = Chem.rdchem.RWMol
RdConformer = Chem.rdchem.Conformer
DType = TypeVar("DType", bound=np.generic)
arrayNx3 = npt.NDArray[DType]
arrayN = npt.NDArray[DType]
