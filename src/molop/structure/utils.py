"""
Author: TMJ
Date: 2025-07-28 12:19:07
LastEditors: TMJ
LastEditTime: 2026-04-18 23:48:18
Description: 请填写简介
"""

from functools import lru_cache

from rdkit import Chem


pt = Chem.GetPeriodicTable()
HETEROATOM = (9, 8, 17, 7, 35, 54, 16, 34, 15)

bond_stereo_mapping = {
    "Z": Chem.rdchem.BondStereo.STEREOZ,
    "E": Chem.rdchem.BondStereo.STEREOE,
}
bond_type_mapping: dict[Chem.rdchem.BondType, float] = {
    Chem.rdchem.BondType.SINGLE: 1.0,
    Chem.rdchem.BondType.DATIVE: 1.0,
    Chem.rdchem.BondType.DATIVEL: 1.0,
    Chem.rdchem.BondType.DATIVER: 1.0,
    Chem.rdchem.BondType.DATIVEONE: 1.0,
    Chem.rdchem.BondType.DOUBLE: 2.0,
    Chem.rdchem.BondType.TRIPLE: 3.0,
    Chem.rdchem.BondType.AROMATIC: 1.5,
    Chem.rdchem.BondType.ZERO: 0.0,
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


@lru_cache(maxsize=8192)
def _canonical_smiles_once(smiles: str) -> str:
    return Chem.CanonSmiles(smiles, useChiral=True)


def canonical_smiles(smiles: str, *, stable: bool = False, max_rounds: int = 8) -> str:
    """Return canonical SMILES with an optional fixed-point fallback.

    Parameters
    ----------
    smiles
        Input SMILES string.
    stable
        When ``False`` (default), use a single RDKit canonicalization pass for speed.
        When ``True``, continue canonicalizing until a fixed point is reached, or stop
        when a cycle is detected.
    max_rounds
        Maximum number of canonicalization rounds when ``stable`` is enabled.
    """
    if max_rounds < 1:
        raise ValueError('max_rounds must be >= 1')

    current = _canonical_smiles_once(smiles)
    if not stable:
        return current

    seen = {smiles, current}
    for _ in range(max_rounds - 1):
        nxt = _canonical_smiles_once(current)
        if nxt == current:
            return current
        if nxt in seen:
            return current
        seen.add(nxt)
        current = nxt
    return current
