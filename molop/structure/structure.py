"""
Including functions related to the structure of molecules
"""

import itertools
from typing import List, Tuple, Union

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdMolTransforms

from molop.structure import geometry
from molop.utils.types import RdMol

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
    mol,
    replacement_mol: RdMol,
    query_mol: Union[RdMol, None] = None,
    bind_idx: int = None,
    threshold=0.75,
    angle_split=10,
    randomSeed=114514,
    start_idx: int = None,
    end_idx: int = None,
):
    """
    Replace the query_mol with replacement_mol in mol.
    """
    if get_under_bonded_number(replacement_mol.GetAtomWithIdx(0)) != 1:
        raise ValueError("Replacement atom should be one bond left")
    if start_idx and end_idx:
        if (
            mol.GetBondBetweenAtoms(start_idx, end_idx)
            and mol.GetBondBetweenAtoms(start_idx, end_idx).GetBondType()
            == Chem.BondType.SINGLE
            and not mol.GetBondBetweenAtoms(start_idx, end_idx).IsInRing()
        ):
            end = end_idx
            start = start_idx
        else:
            raise ValueError("start_idx and end_idx should be bonded with single but not in a ring.")
    else:
        if query_mol is None:
            raise ValueError("start_idx and end_idx should be provided.")
        else:
            query_mol_ = Chem.MolFromSmarts(Chem.MolToSmarts(query_mol))
            queried_idx_list = mol.GetSubstructMatches(query_mol_)
            # print(queried_idx_list)
            if len(queried_idx_list) == 0:
                raise ValueError("No substruct match found.")
            if bind_idx:
                for queried_idx in queried_idx_list:
                    if bind_idx == queried_idx[0]:
                        unique_queried_idx = queried_idx
                        break
                else:
                    raise ValueError(f"No substruct mapping to bind_idx {bind_idx}.")
            else:
                unique_queried_idx = queried_idx_list[0]
            # print(unique_queried_idx)
            for atom in mol.GetAtomWithIdx(unique_queried_idx[0]).GetNeighbors():
                if (
                    atom.GetIdx() not in mol.GetSubstructMatch(query_mol_)
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
        + pt.GetRcovalent(replacement_mol.GetAtomWithIdx(0).GetAtomicNum()),
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
    replacement = Chem.AddHs(replacement_mol, addCoords=True)
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
    best_angle, best_crowding_score = 0, 0
    for angle in np.linspace(0, 2 * np.pi, angle_split):
        rdMolTransforms.SetDihedralRad(
            rmol.GetConformer(),
            first_atom_idx,
            lines_idx[0][0],
            0,
            forth_atom_idx,
            angle,
        )
        crowding_score = get_crowding_socre(rmol, threshold=threshold)
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
        if mol.GetBondBetweenAtoms(start_atom.GetIdx(), end_atom.GetIdx()):
            continue
        if distances[start_atom.GetIdx()][end_atom.GetIdx()] < threshold * (
            pt.GetRcovalent(start_atom.GetAtomicNum())
            + pt.GetRcovalent(end_atom.GetAtomicNum())
        ):
            return False
    return True


def get_crowding_socre(mol, threshold: float):
    score = 0
    distances = Chem.Get3DDistanceMatrix(mol)
    for start_atom, end_atom in itertools.combinations(mol.GetAtoms(), 2):
        if mol.GetBondBetweenAtoms(start_atom.GetIdx(), end_atom.GetIdx()):
            continue
        ideal_distance = pt.GetRcovalent(start_atom.GetAtomicNum()) + pt.GetRcovalent(
            end_atom.GetAtomicNum()
        )
        distance = distances[start_atom.GetIdx()][end_atom.GetIdx()]
        if distance < threshold * ideal_distance:
            return -float("inf")
        score += np.log(distance / ideal_distance)
    return score


def attempt_replacement(
    mol,
    query: Union[str, RdMol],
    replacement: Union[str, RdMol] = None,
    bind_idx: int = None,
    replace_all=False,
    attempt_num=10,
    crowding_threshold=0.75,
    angle_split=10,
    randomSeed=114514,
    start_idx: int = None,
    end_idx: int = None,
):
    random_seeds = list(range(1, attempt_num + 1))
    np.random.RandomState(randomSeed).shuffle(random_seeds)

    if isinstance(query, str):
        query_mol = Chem.MolFromSmarts(query)
    elif isinstance(query, (RdMol, Chem.RWMol)):
        query_mol = Chem.MolFromSmarts(Chem.MolToSmarts(query))
    if isinstance(replacement, str):
        replacement_mol = Chem.MolFromSmiles(replacement)
    elif isinstance(replacement, (RdMol, Chem.RWMol)):
        replacement_mol = replacement

    for random_seed in random_seeds:
        if replace_all:
            if replacement_mol.HasSubstructMatch(query_mol):
                raise RuntimeError(
                    f"Endless loop: replacement '{Chem.MolToSmiles(replacement_mol)}' contains query '{Chem.MolToSmiles(query_mol)}'"
                )

            new_mol = replace_mol(
                mol=mol,
                query_mol=query_mol,
                replacement_mol=replacement_mol,
                bind_idx=None,
                threshold=crowding_threshold,
                angle_split=angle_split,
                randomSeed=random_seed,
            )
            while new_mol.HasSubstructMatch(query_mol):
                # print(new_mol.GetSubstructMatches(query_mol))
                new_mol = replace_mol(
                    mol=new_mol,
                    query_mol=query_mol,
                    replacement_mol=replacement_mol,
                    bind_idx=None,
                    threshold=crowding_threshold,
                    angle_split=angle_split,
                    randomSeed=random_seed,
                )
        else:
            new_mol = replace_mol(
                mol=mol,
                query_mol=query_mol,
                replacement_mol=replacement_mol,
                bind_idx=bind_idx,
                threshold=crowding_threshold,
                angle_split=angle_split,
                randomSeed=random_seed,
                start_idx=start_idx,
                end_idx=end_idx,
            )
        if check_crowding(new_mol, threshold=crowding_threshold):
            return new_mol
    else:
        raise RuntimeError(
            f"replacement '{Chem.MolToSmiles(replacement_mol)}' is too big in threshold {crowding_threshold}, try smaller threshold"
        )


def fix_geometry(replacement: RdMol, bind_type: int, randomSeed=114514):
    r = Chem.RWMol(replacement)
    idx = r.AddAtom(Chem.Atom("At"))
    r.AddBond(0, idx, Chem.BondType.SINGLE)
    r.UpdatePropertyCache()
    Chem.SanitizeMol(r)
    if replacement.GetNumConformers() == 0:
        AllChem.EmbedMolecule(r, randomSeed=randomSeed)
    else:
        cmap = {
            atom_idx: replacement.GetConformer().GetAtomPosition(atom_idx)
            for atom_idx in range(replacement.GetNumAtoms())
        }
        AllChem.EmbedMolecule(r, randomSeed=randomSeed, coordMap=cmap)
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
