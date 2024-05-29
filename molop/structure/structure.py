"""
Including functions related to the structure of molecules
"""

import itertools
from typing import List, Tuple

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdMolTransforms

from ..utils.types import RdMol
from . import geometry

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


def get_formal_num_radicals(mol: RdMol) -> List[int]:
    """
    Get formal spin of mol.
    """
    return [atom.GetNumRadicalElectrons() for atom in mol.GetAtoms()]


def get_resonance_structures(rdmol, flags=0):
    suppl = Chem.ResonanceMolSupplier(rdmol, flags)
    return [mol for mol in suppl]


def structure_score(rdmol):
    return sum(atom.GetNumRadicalElectrons() for atom in rdmol.GetAtoms()) + sum(
        abs(atom.GetFormalCharge()) for atom in rdmol.GetAtoms()
    )


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
    mol: Chem.rdchem.Mol,
    query_smi: str,
    replacement_smi: str,
    bind_idx: int = None,
    threshold=0.75,
    randomSeed=114514,
) -> Chem.rdchem.Mol:
    """
    Replace the query_smi with replacement_smi in mol.
    """
    query = Chem.MolFromSmarts(query_smi)
    replacement = Chem.MolFromSmiles(replacement_smi)
    if get_under_bonded_number(replacement.GetAtomWithIdx(0)) != 1:
        raise ValueError("Replacement atom should be one bond left")
    queried_idx_list = mol.GetSubstructMatches(query)
    if len(queried_idx_list) == 0:
        raise ValueError("No substruct match found.")
    if bind_idx:
        for queried_idx in queried_idx_list:
            if queried_idx[0] == bind_idx:
                unique_queried_idx = queried_idx
                break
        else:
            raise ValueError(f"No substruct mapping to bind_idx {bind_idx}.")
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
    else:
        raise ValueError("Bond between skeleton and replacement is not single.")
    # start: the idx of atom that bind to the substruct to be replaced
    # end: the idx of atom in substruct to be replaced, and bind to the main structure

    origin_mol = Chem.RWMol(mol)
    geometry.standard_orient(origin_mol, [start, end])
    rdMolTransforms.SetBondLength(
        origin_mol.GetConformer(),
        start,
        end,
        pt.GetRcovalent(origin_mol.GetAtomWithIdx(start).GetAtomicNum())
        + pt.GetRcovalent(replacement.GetAtomWithIdx(0).GetAtomicNum()),
    )
    # build skeleton with queried substruct removed
    # considering the original mol may have been split into frags
    origin_mol.RemoveBond(start, end)
    frags = [
        frag
        for frag, frag_idx in zip(
            Chem.GetMolFrags(origin_mol, asMols=True), Chem.GetMolFrags(origin_mol)
        )
        if end not in frag_idx
    ]  # get all frags except the one that contains the end atom
    skeleton = frags[0]
    for frag in frags[1:]:
        skeleton = Chem.CombineMols(skeleton, frag)

    # build replacement conformer
    replacement = Chem.AddHs(replacement)
    replacement.GetAtomWithIdx(0).SetNumRadicalElectrons(0)
    rr = fix_geometry(
        replacement,
        bind_type=mol.GetAtomWithIdx(start).GetAtomicNum(),
        randomSeed=randomSeed,
    )

    # combine skeleton and replacement
    new_mol = Chem.CombineMols(rr, skeleton)
    # find the new index mapping of start and end (end in the replacement now)
    lines_idx = [
        (idx, x)
        for idx, (x, y, z) in enumerate(new_mol.GetConformer().GetPositions())
        if abs(y - 0) <= 0.001 and abs(z - 0) <= 0.001
    ]
    lines_idx.sort(key=lambda x: x[1])
    rmol = Chem.RWMol(new_mol)
    rmol.AddBond(lines_idx[0][0], 0, Chem.rdchem.BondType.SINGLE)
    skeleton.GetAtomWithIdx(
        lines_idx[0][0] - len([atom for atom in rr.GetAtoms()])
    ).SetNumRadicalElectrons(mol.GetAtomWithIdx(start).GetNumRadicalElectrons())
    rmol.GetAtomWithIdx(lines_idx[0][0]).SetNumRadicalElectrons(
        mol.GetAtomWithIdx(start).GetNumRadicalElectrons()
    )
    replacement.GetAtomWithIdx(0).SetNumRadicalElectrons(0)
    # find the idx list of replacement
    idx_list = [
        idx_list for idx_list in rmol.GetSubstructMatches(replacement) if 0 in idx_list
    ]
    for forth_atom in rmol.GetAtomWithIdx(0).GetNeighbors():
        if forth_atom.GetIdx() in idx_list[0]:
            forth_atom_idx = forth_atom.GetIdx()
            break
    for first_atom in rmol.GetAtomWithIdx(lines_idx[0][0]).GetNeighbors():
        if first_atom.GetIdx() not in idx_list[0]:
            first_atom_idx = first_atom.GetIdx()
            break
    Chem.SanitizeMol(rmol)
    best_angle = 0
    best_crowding_score = 0
    for angle in np.linspace(0, 2 * np.pi, 10):
        rdMolTransforms.SetDihedralRad(
            rmol.GetConformer(),
            first_atom_idx,
            lines_idx[0][0],
            0,
            forth_atom_idx,
            angle,
        )
        crowding_score = get_crowding_socre(rmol)
        if crowding_score > best_crowding_score:
            best_angle = angle
            best_crowding_score = crowding_score
    rdMolTransforms.SetDihedralRad(
        rmol.GetConformer(),
        first_atom_idx,
        lines_idx[0][0],
        0,
        forth_atom_idx,
        best_angle,
    )
    geometry.standard_orient(rmol, [0, 1, 2])
    # recover the original atom index of the skeleton
    if len([atom for atom in skeleton.GetAtoms()]) > 50:
        return rmol
    for indice in mol.GetSubstructMatches(skeleton):
        if unique_queried_idx[0] not in indice:
            origin_indice = list(indice) + [unique_queried_idx[0]]
            break
    for indice in rmol.GetSubstructMatches(skeleton):
        if 0 not in indice:
            transed_indice = list(indice) + [0]
            break
    zip_indice = list(zip(origin_indice, transed_indice))
    zip_indice.sort(key=lambda x: x[0])
    mapping = [item[1] for item in zip_indice]
    return reset_atom_index(rmol, mapping)


def check_crowding(mol, threshold=0.6):
    distances = Chem.Get3DDistanceMatrix(mol)
    for start_atom, end_atom in itertools.combinations(mol.GetAtoms(), 2):
        if distances[start_atom.GetIdx()][end_atom.GetIdx()] < threshold * (
            pt.GetRcovalent(start_atom.GetAtomicNum())
            + pt.GetRcovalent(end_atom.GetAtomicNum())
        ):
            return False
    return True


def get_crowding_socre(mol):
    score = 0
    distances = Chem.Get3DDistanceMatrix(mol)
    for start_atom, end_atom in itertools.combinations(mol.GetAtoms(), 2):
        ideal_distance = pt.GetRcovalent(start_atom.GetAtomicNum()) + pt.GetRcovalent(
            end_atom.GetAtomicNum()
        )
        distance = distances[start_atom.GetIdx()][end_atom.GetIdx()]
        score += np.exp(-ideal_distance / distance)
    return score


def attempt_replacement(
    mol,
    query_smi: str,
    replacement_smi: str,
    bind_idx: int = None,
    replace_all=False,
    attempt_num=10,
    crowding_threshold=0.75,
    randomSeed=114514,
):
    random_seeds = list(range(1, attempt_num + 1))
    np.random.RandomState(randomSeed).shuffle(random_seeds)

    query = Chem.MolFromSmarts(query_smi)
    for random_seed in random_seeds:
        if replace_all:
            match_num = len(mol.GetSubstructMatches(query))
            if len(mol.GetSubstructMatches(query)) > match_num:
                raise RuntimeError(
                    f"Endless loop: {replacement_smi} contains {query_smi}"
                )
            new_mol = replace_mol(
                mol,
                query_smi,
                replacement_smi,
                bind_idx=None,
                threshold=crowding_threshold,
                randomSeed=random_seed,
            )
            while new_mol.GetSubstructMatches(query):
                new_mol = replace_mol(
                    new_mol,
                    query_smi,
                    replacement_smi,
                    bind_idx=None,
                    threshold=crowding_threshold,
                    randomSeed=random_seed,
                )
        else:
            new_mol = replace_mol(
                mol,
                query_smi,
                replacement_smi,
                bind_idx=bind_idx,
                threshold=crowding_threshold,
                randomSeed=random_seed,
            )
        if check_crowding(new_mol, threshold=crowding_threshold):
            return new_mol
    else:
        raise RuntimeError(f"replacement {replacement_smi} is too big")


def fix_geometry(replacement, bind_type: int, randomSeed=114514):
    r = Chem.RWMol(replacement)
    idx = r.AddAtom(Chem.Atom("At"))
    r.AddBond(0, idx, Chem.BondType.SINGLE)
    r.UpdatePropertyCache()
    Chem.SanitizeMol(r)
    AllChem.EmbedMolecule(r, randomSeed=randomSeed)
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
    Chem.SanitizeMol(rwmol)
    return rwmol
