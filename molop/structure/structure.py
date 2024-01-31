"""
Author: TMJ
Date: 2023-06-26 14:24:42
LastEditors: TMJ
LastEditTime: 2023-06-27 21:03:36
Description: 请填写简介
"""

import itertools
from copy import deepcopy
from typing import List, Tuple

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolTransforms

from . import geometry
from ..utils.types import RdConformer, RdMol

RDLogger.DisableLog("rdApp.*")

pt = Chem.GetPeriodicTable()

bond_list = [
    Chem.rdchem.BondType.UNSPECIFIED,
    Chem.rdchem.BondType.SINGLE,
    Chem.rdchem.BondType.DOUBLE,
    Chem.rdchem.BondType.TRIPLE,
    Chem.rdchem.BondType.QUADRUPLE,
    Chem.rdchem.BondType.QUINTUPLE,
    Chem.rdchem.BondType.HEXTUPLE,
    Chem.rdchem.BondType.ONEANDAHALF,
    Chem.rdchem.BondType.TWOANDAHALF,
    Chem.rdchem.BondType.THREEANDAHALF,
    Chem.rdchem.BondType.FOURANDAHALF,
    Chem.rdchem.BondType.FIVEANDAHALF,
    Chem.rdchem.BondType.AROMATIC,
    Chem.rdchem.BondType.IONIC,
    Chem.rdchem.BondType.HYDROGEN,
    Chem.rdchem.BondType.THREECENTER,
    Chem.rdchem.BondType.DATIVEONE,
    Chem.rdchem.BondType.DATIVE,
    Chem.rdchem.BondType.DATIVEL,
    Chem.rdchem.BondType.DATIVER,
    Chem.rdchem.BondType.OTHER,
    Chem.rdchem.BondType.ZERO,
]


def get_bond_pairs(mol: RdMol) -> List[Tuple[int, int, int]]:
    """
    Get bond pair of mol.
    """
    return [
        (
            bond.GetBeginAtomIdx(),
            bond.GetEndAtomIdx(),
            bond_list.index(bond.GetBondType()),
        )
        for bond in mol.GetBonds()
    ]


def get_formal_charges(mol: RdMol) -> List[int]:
    """
    Get formal charge of mol.
    """
    return [atom.GetFormalCharge() for atom in mol.GetAtoms()]


def get_formal_spins(mol: RdMol) -> List[int]:
    """
    Get formal spin of mol.
    """
    return [atom.GetNumRadicalElectrons() for atom in mol.GetAtoms()]


def get_resonance_structures(rdmol, flags=0):
    suppl = Chem.ResonanceMolSupplier(rdmol, flags)
    return [mol for mol in suppl]


def structure_score(rdmol):
    s = 0
    s += sum(atom.GetNumRadicalElectrons() for atom in rdmol.GetAtoms())
    s += sum(abs(atom.GetFormalCharge()) for atom in rdmol.GetAtoms())
    return s


def check_mol_equal(rdmol_1, rdmol_2):
    return rdmol_1.HasSubstructMatch(rdmol_2) and rdmol_2.HasSubstructMatch(rdmol_1)


def get_sub_mol(origin_mol: RdMol, scale: List[int]):
    sub_mol = Chem.RWMol(Chem.MolFromSmiles(""))
    for idx in scale:
        atom = origin_mol.GetAtomWithIdx(idx)
        new_atom = Chem.Atom(atom.GetSymbol())
        new_atom.SetFormalCharge(atom.GetFormalCharge())
        sub_mol.AddAtom(new_atom)

    for bond in origin_mol.GetBonds():
        begin_atom = bond.GetBeginAtom()
        end_atom = bond.GetEndAtom()
        if begin_atom.GetIdx() not in scale or end_atom.GetIdx() not in scale:
            continue
        bond_type = bond.GetBondType()
        sub_mol.AddBond(
            begin_atom.GetIdx(),
            end_atom.GetIdx(),
            bond_type,
        )

    return sub_mol


def replace_mol(mol, query_smi: str, replacement_smi: str, bind_idx: int = None):
    query = Chem.MolFromSmarts(query_smi)
    queried_idx_list = mol.GetSubstructMatches(query)
    start = None
    if bind_idx is None:
        for atom in mol.GetAtomWithIdx(queried_idx_list[0][0]).GetNeighbors():
            if (
                atom.GetIdx()
                not in mol.GetSubstructMatch(Chem.MolFromSmarts(query_smi))
                and mol.GetBondBetweenAtoms(
                    atom.GetIdx(), queried_idx_list[0][0]
                ).GetBondType()
                == Chem.BondType.SINGLE
            ):
                start = atom.GetIdx()
                end = queried_idx_list[0][0]
                break
    else:
        for queried_idx in queried_idx_list:
            if (
                mol.GetBondBetweenAtoms(queried_idx[0], bind_idx)
                and mol.GetBondBetweenAtoms(queried_idx[0], bind_idx).GetBondType()
                == Chem.BondType.SINGLE
            ):
                start = bind_idx
                end = queried_idx[0]
                break
    if start is None:
        raise ValueError("Can not replace")

    origin_mol = Chem.RWMol(mol)
    geometry.standard_orient(origin_mol, [start, end])
    replacement = Chem.MolFromSmiles(replacement_smi)
    replacement = Chem.AddHs(replacement)
    AllChem.EmbedMolecule(replacement)
    replacement.GetAtomWithIdx(0).SetNumRadicalElectrons(0)
    rr = fix_geometry(
        replacement,
    )
    new_mol = Chem.ReplaceSubstructs(
        rr,
        Chem.MolFromSmarts("[At]"),
        origin_mol,
        replacementConnectionPoint=start,
    )[0]
    lines_idx = [
        (idx, x)
        for idx, (x, y, z) in enumerate(new_mol.GetConformer().GetPositions())
        if abs(y - 0) <= 0.001 and abs(z - 0) <= 0.001
    ]
    lines_idx.sort(key=lambda x: x[1])
    rmol = Chem.RWMol(new_mol)
    rmol.RemoveBond(lines_idx[0][0], lines_idx[1][0])
    rmol.RemoveAtom(lines_idx[1][0])
    Chem.MolToMolFile(rmol, "temp.sdf")
    idx_list = [
        idx_list for idx_list in rmol.GetSubstructMatches(replacement) if 0 in idx_list
    ]
    for start_atom in rmol.GetAtomWithIdx(0).GetNeighbors():
        if start_atom.GetIdx() not in idx_list[0]:
            start = start_atom.GetIdx()
            end = 0
            break
    Chem.SanitizeMol(rmol)
    rdMolTransforms.SetBondLength(
        rmol.GetConformer(),
        end,
        start,
        pt.GetRcovalent(new_mol.GetAtomWithIdx(start).GetAtomicNum())
        + pt.GetRcovalent(new_mol.GetAtomWithIdx(end).GetAtomicNum()),
    )
    geometry.standard_orient(rmol, [0, 1, 2])
    return rmol


def check_crowding(mol):
    distances = Chem.Get3DDistanceMatrix(mol)
    for start_atom, end_atom in itertools.combinations(mol.GetAtoms(), 2):
        if distances[start_atom.GetIdx()][end_atom.GetIdx()] < 0.75 * (
            pt.GetRcovalent(start_atom.GetAtomicNum())
            + pt.GetRcovalent(end_atom.GetAtomicNum())
        ):
            return False
    return True


def attempt_replacement(
    mol,
    query_smi: str,
    replacement_smi: str,
    bind_idx: int = None,
    replace_all=False,
    attempt_num=10,
):
    query = Chem.MolFromSmarts(query_smi)
    for i in range(attempt_num):
        if replace_all:
            match_num = len(mol.GetSubstructMatches(query))
            new_mol = replace_mol(mol, query_smi, replacement_smi, bind_idx=None)
            if len(mol.GetSubstructMatches(query)) > match_num:
                raise RuntimeError(
                    f"Endless loop: {replacement_smi} contains {query_smi}"
                )
            while new_mol.GetSubstructMatches(query):
                new_mol = replace_mol(
                    new_mol, query_smi, replacement_smi, bind_idx=None
                )
        else:
            new_mol = replace_mol(mol, query_smi, replacement_smi, bind_idx)
        if check_crowding(new_mol):
            return new_mol
    raise RuntimeError(f"replacement {replacement_smi} is too big")


def fix_geometry(replacement):
    r = Chem.RWMol(replacement)
    idx = r.AddAtom(Chem.Atom("At"))
    r.AddBond(0, idx, Chem.BondType.SINGLE)
    r.UpdatePropertyCache()
    AllChem.EmbedMolecule(r)
    geometry.standard_orient(r, [idx, 0])
    # geometry.translate_anchor(r, 0)
    # geometry.translate_mol(r, vec)
    return r
