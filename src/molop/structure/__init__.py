"""
Author: TMJ
Date: 2025-10-09 15:23:26
LastEditors: TMJ
LastEditTime: 2025-12-23 22:16:40
Description: 请填写简介
"""

from rdkit import Chem

from .FormatConverter import omol_to_rdmol
from .GraphReconstruction import xyz2omol
from .StructureTransformation import make_dative_bonds


def xyz_to_rdmol(
    xyz_block: str,
    total_charge: int = 0,
    total_radical_electrons: int = 0,
    make_dative: bool = True,
) -> Chem.Mol | None:
    """
    Convert an XYZ block to a RDKit molecule.

    Parameters:
        xyz_block (str): The XYZ block of a molecule.
        total_charge (int, optional): The total charge of the molecule. Defaults to 0.
        total_radical_electrons (int, optional): The total radical electrons of the molecule. Defaults to 0.
        make_dative (bool, optional): Whether to make dative bonds. Defaults to True.

    Returns:
        Chem.Mol: The RDKit molecule.
    """
    possible_omol = xyz2omol(xyz_block, total_charge, total_radical_electrons)
    if possible_omol is None:
        return None
    rdmol_opt = omol_to_rdmol(possible_omol, total_charge, total_radical_electrons)
    if rdmol_opt is None:
        return None
    if make_dative:
        return make_dative_bonds(Chem.RWMol(rdmol_opt))
    return rdmol_opt
