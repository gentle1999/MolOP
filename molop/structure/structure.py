"""
Including functions related to the structure of molecules
"""

import itertools
from copy import deepcopy
from typing import List, Tuple

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdFMCS
from rdkit.Chem import rdMolTransforms

from . import geometry, structure_recovery
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


def get_under_bonded_number(atom: Chem.Atom) -> int:
    """
    Get the number of under-bonded electrons of atom.
    """
    return (
        pt.GetDefaultValence(atom.GetAtomicNum())
        + atom.GetFormalCharge()
        - atom.GetTotalValence()
    )


def replace_mol(
    mol: Chem.rdchem.Mol, query_smi: str, replacement_smi: str, bind_idx: int = None
) -> Chem.rdchem.Mol:
    query = Chem.MolFromSmarts(query_smi)
    replacement = Chem.MolFromSmiles(replacement_smi)
    if get_under_bonded_number(replacement.GetAtomWithIdx(0)) != 1:
        raise ValueError("Replacement atom should be one bond left")
    queried_idx_list = mol.GetSubstructMatches(query)
    if len(queried_idx_list) == 0:
        raise ValueError("No substruct match found.")
    start = None
    if bind_idx:
        for queried_idx in queried_idx_list:
            if queried_idx[0] == bind_idx:
                unique_queried_idx = queried_idx
    else:
        unique_queried_idx = queried_idx_list[0]
    for atom in mol.GetAtomWithIdx(unique_queried_idx[0]).GetNeighbors():
        if (
            atom.GetIdx() not in mol.GetSubstructMatch(Chem.MolFromSmarts(query_smi))
            and mol.GetBondBetweenAtoms(
                atom.GetIdx(), unique_queried_idx[0]
            ).GetBondType()
            == Chem.BondType.SINGLE
        ):
            start = atom.GetIdx()
            end = unique_queried_idx[0]
            break
    if start is None:
        raise ValueError("Can not replace")

    origin_mol = Chem.RWMol(mol)
    geometry.standard_orient(origin_mol, [start, end])
    rdMolTransforms.SetBondLength(
        origin_mol.GetConformer(),
        start,
        end,
        pt.GetRcovalent(origin_mol.GetAtomWithIdx(start).GetAtomicNum())
        + pt.GetRcovalent(replacement.GetAtomWithIdx(0).GetAtomicNum()),
    )
    origin_mol.RemoveBond(start, end)
    frags = [
        frag
        for frag, frag_idx in zip(
            Chem.GetMolFrags(origin_mol, asMols=True), Chem.GetMolFrags(origin_mol)
        )
        if end not in frag_idx
    ]
    skeleton = frags[0]
    for frag in frags[1:]:
        skeleton = Chem.CombineMols(skeleton, frag)

    replacement = Chem.AddHs(replacement)
    AllChem.EmbedMolecule(replacement)
    rr = fix_geometry(replacement, bind_type=mol.GetAtomWithIdx(start).GetAtomicNum())
    new_mol = Chem.CombineMols(rr, skeleton)
    lines_idx = [
        (idx, x)
        for idx, (x, y, z) in enumerate(new_mol.GetConformer().GetPositions())
        if abs(y - 0) <= 0.001 and abs(z - 0) <= 0.001
    ]
    lines_idx.sort(key=lambda x: x[1])
    rmol = Chem.RWMol(new_mol)
    rmol.AddBond(lines_idx[0][0], 0, Chem.rdchem.BondType.SINGLE)
    rmol.GetAtomWithIdx(0).SetNumRadicalElectrons(
        rmol.GetAtomWithIdx(0).GetNumRadicalElectrons() - 1
    )
    rmol.GetAtomWithIdx(lines_idx[0][0]).SetNumRadicalElectrons(
        rmol.GetAtomWithIdx(lines_idx[0][0]).GetNumRadicalElectrons() - 1
    )
    replacement.GetAtomWithIdx(0).SetNumRadicalElectrons(0)
    idx_list = [
        idx_list for idx_list in rmol.GetSubstructMatches(replacement) if 0 in idx_list
    ]
    for start_atom in rmol.GetAtomWithIdx(0).GetNeighbors():
        if start_atom.GetIdx() not in idx_list[0]:
            start = start_atom.GetIdx()
            end = 0
            break
    Chem.SanitizeMol(rmol)
    geometry.standard_orient(rmol, [0, 1, 2])
    if len([atom for atom in rmol.GetAtoms()]) > 50:
        return rmol
    public_part = Chem.MolFromSmarts(rdFMCS.FindMCS((origin_mol, rmol)).smartsString)
    for indice in mol.GetSubstructMatches(public_part):
        if unique_queried_idx[0] not in indice:
            origin_indice = list(indice) + [unique_queried_idx[0]]
            break
    for indice in rmol.GetSubstructMatches(public_part):
        if 0 not in indice:
            transed_indice = list(indice) + [0]
            break
    zip_indice = list(zip(origin_indice, transed_indice))
    zip_indice.sort(key=lambda x: x[0])
    mapping = [item[1] for item in zip_indice]
    return reset_atom_index(rmol, mapping)


def check_crowding(mol, threshold=0.75):
    distances = Chem.Get3DDistanceMatrix(mol)
    for start_atom, end_atom in itertools.combinations(mol.GetAtoms(), 2):
        if distances[start_atom.GetIdx()][end_atom.GetIdx()] < threshold * (
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


def fix_geometry(replacement, bind_type: int):
    r = Chem.RWMol(replacement)
    idx = r.AddAtom(Chem.Atom("At"))
    r.AddBond(0, idx, Chem.BondType.SINGLE)
    r.UpdatePropertyCache()
    AllChem.EmbedMolecule(r)
    geometry.standard_orient(r, [idx, 0])
    rdMolTransforms.SetBondLength(
        r.GetConformer(),
        idx,
        0,
        pt.GetRcovalent(r.GetAtomWithIdx(0).GetAtomicNum())
        + pt.GetRcovalent(bind_type),
    )
    r.RemoveAtom(idx)
    return r


def reset_atom_index(mol: Chem.rdchem.Mol, mapping: List[int]) -> Chem.rdchem.Mol:
    # drop the bonds
    rwmol = Chem.RWMol()
    # new index
    new_idx = list(mapping) + [
        idx for idx, atom in enumerate(mol.GetAtoms()) if idx not in mapping
    ]
    mol_copy = Chem.RWMol(mol)
    for i, idx in enumerate(new_idx):
        mol_copy.GetConformer().SetAtomPosition(
            i, mol.GetConformer().GetAtomPosition(idx)
        )

    for idx in new_idx:
        rwmol.AddAtom(mol.GetAtomWithIdx(idx))

    bond_pairs = get_bond_pairs(mol)
    for bond_pair in bond_pairs:
        rwmol.AddBond(
            new_idx.index(bond_pair[0]),
            new_idx.index(bond_pair[1]),
            bond_list[bond_pair[2]],
        )
    rwmol.AddConformer(mol_copy.GetConformer(0))
    return rwmol
