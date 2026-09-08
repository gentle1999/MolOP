"""Structure validity and geometry-quality checks."""

from __future__ import annotations

import itertools

from rdkit import Chem
from rdkit.Chem import rdForceFieldHelpers

from molop.utils.types import RdMol

from .utils import estimate_bond_length


_SUPPORTED_DISTANCE_BONDS = (
    Chem.rdchem.BondType.SINGLE,
    Chem.rdchem.BondType.DATIVE,
    Chem.rdchem.BondType.DATIVEL,
    Chem.rdchem.BondType.DATIVER,
    Chem.rdchem.BondType.DATIVEONE,
    Chem.rdchem.BondType.DOUBLE,
    Chem.rdchem.BondType.TRIPLE,
)


def check_crowding(mol: RdMol, threshold: float = 0.5) -> bool:
    """Return whether all atom pairs satisfy the threshold distance."""
    distances = Chem.Get3DDistanceMatrix(mol)
    for start_atom, end_atom in itertools.combinations(mol.GetAtoms(), 2):
        bond = mol.GetBondBetweenAtoms(start_atom.GetIdx(), end_atom.GetIdx())
        bond_type = (
            bond.GetBondType()
            if bond is not None and bond.GetBondType() in _SUPPORTED_DISTANCE_BONDS
            else Chem.rdchem.BondType.SINGLE
        )
        if distances[start_atom.GetIdx()][end_atom.GetIdx()] < threshold * estimate_bond_length(
            start_atom.GetAtomicNum(),
            end_atom.GetAtomicNum(),
            bond_type,
        ):
            return False
    return True


def get_crowding_score(mol: RdMol) -> float:
    """Return the negative UFF energy used to rank replacement conformers."""
    Chem.SanitizeMol(mol)
    return -rdForceFieldHelpers.UFFGetMoleculeForceField(mol).CalcEnergy()
