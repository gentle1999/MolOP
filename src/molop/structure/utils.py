"""
Author: TMJ
Date: 2025-07-28 12:19:07
LastEditors: TMJ
LastEditTime: 2025-12-14 19:05:25
Description: 请填写简介
"""

from rdkit import Chem


pt = Chem.GetPeriodicTable()
HETEROATOM = (9, 8, 17, 7, 35, 54, 16, 34, 15)

bond_stereo_mapping = {
    "Z": Chem.rdchem.BondStereo.STEREOZ,
    "E": Chem.rdchem.BondStereo.STEREOE,
}
bond_type_mapping = {
    Chem.rdchem.BondType.SINGLE: 1,
    Chem.rdchem.BondType.DATIVE: 1,
    Chem.rdchem.BondType.DATIVEL: 1,
    Chem.rdchem.BondType.DATIVER: 1,
    Chem.rdchem.BondType.DATIVEONE: 1,
    Chem.rdchem.BondType.DOUBLE: 2,
    Chem.rdchem.BondType.TRIPLE: 3,
    Chem.rdchem.BondType.AROMATIC: 1.5,
    Chem.rdchem.BondType.ZERO: 0,
    Chem.rdchem.BondType.ONEANDAHALF: 1.5,
    Chem.rdchem.BondType.TWOANDAHALF: 2.5,
}
bond_list = list(Chem.rdchem.BondType.values.values())
bond_stereo_list = list(Chem.rdchem.BondStereo.values.values())


def estimate_bond_length(
    begin_atomic_num: int,
    end_atomic_num: int,
    bond_type: Chem.rdchem.BondType = Chem.rdchem.BondType.SINGLE,
) -> float:
    """
    Estimate the bond length between two atoms based on their atomic numbers and bond type.
    Parameters:
        begin_atomic_num (int):
            The atomic number of the first atom.
        end_atomic_num (int):
            The atomic number of the second atom.
        bond_type (Chem.rdchem.BondType):
            The bond type between the two atoms.
    Returns:
        float:
            The estimated bond length between the two atoms.
    """
    single_length = pt.GetRcovalent(begin_atomic_num) + pt.GetRcovalent(end_atomic_num)
    if bond_type in (
        Chem.rdchem.BondType.SINGLE,
        Chem.rdchem.BondType.DATIVE,
        Chem.rdchem.BondType.DATIVEL,
        Chem.rdchem.BondType.DATIVER,
        Chem.rdchem.BondType.DATIVEONE,
    ):
        return single_length
    elif bond_type == Chem.rdchem.BondType.DOUBLE:
        return single_length / 1.15
    elif bond_type == Chem.rdchem.BondType.TRIPLE:
        return single_length / 1.15 / 1.1
    else:
        raise ValueError(f"Unsupported bond type {bond_type}.")


def canonical_smiles(smiles: str):
    original_smi = smiles
    viewed_smi = {original_smi: 1}
    while original_smi != (canonical_smi := Chem.CanonSmiles(original_smi, useChiral=True)) and (
        canonical_smi not in viewed_smi or viewed_smi[canonical_smi] < 2
    ):
        original_smi = canonical_smi
        if original_smi not in viewed_smi:
            viewed_smi[original_smi] = 1
        else:
            viewed_smi[original_smi] += 1
    else:
        return original_smi
