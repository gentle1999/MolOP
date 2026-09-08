"""Compatibility facade for structure topology and substituent editing."""

from rdkit import Chem

from .editing import (
    SubstituentReplacementError,
    replace_multisite_substituent,
    replace_substituent,
)
from .topology import (
    build_mol_from_atoms_and_bonds,
    get_bond_pairs,
    get_formal_charges,
    get_formal_num_radicals,
    get_total_charge,
    get_total_multiplicity,
    get_total_num_radical,
    reset_atom_index,
)
from .validation import check_crowding, get_crowding_score


__all__ = [
    "Chem",
    "build_mol_from_atoms_and_bonds",
    "check_crowding",
    "get_bond_pairs",
    "get_crowding_score",
    "get_formal_charges",
    "get_formal_num_radicals",
    "get_total_charge",
    "get_total_multiplicity",
    "get_total_num_radical",
    "replace_substituent",
    "replace_multisite_substituent",
    "reset_atom_index",
    "SubstituentReplacementError",
]
