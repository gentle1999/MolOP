"""
Author: TMJ
Date: 2023-06-26 14:24:42
LastEditors: TMJ
LastEditTime: 2023-06-27 21:03:36
Description: 请填写简介
"""

import os
import shutil
from copy import deepcopy
from typing import List, Tuple


from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from tqdm import tqdm

from . import geometry
from .types import RdConformer, RdMol

RDLogger.DisableLog("rdApp.*")

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


def get_bond_matrix(mol: RdMol) -> List[List[int]]:
    """
    Get bond matrix of mol.
    """
    adjacency_matrix = Chem.rdmolops.GetAdjacencyMatrix(mol)
    for bond in mol.GetBonds():
        adjacency_matrix[
            bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        ] = bond_list.index(bond.GetBondType())
        adjacency_matrix[
            bond.GetEndAtomIdx(), bond.GetBeginAtomIdx()
        ] = bond_list.index(bond.GetBondType())
    return adjacency_matrix


def get_bond_pair(mol: RdMol) -> List[Tuple[Tuple[int, int], int]]:
    """
    Get bond pair of mol.
    """
    bond_pair = []
    for bond in mol.GetBonds():
        bond_pair.extend(
            (
                bond.GetBeginAtomIdx(),
                bond.GetEndAtomIdx(),
                bond_list.index(bond.GetBondType()),
            )
        )
        bond_pair.extend(
            (
                bond.GetEndAtomIdx(),
                bond.GetBeginAtomIdx(),
                bond_list.index(bond.GetBondType()),
            )
        )
    return bond_pair


def bond_matrix_2_mol(atom_numbers: List[int], bond_matrix: List[List[int]]) -> RdMol:
    # TODO 速度有待优化
    molecule = Chem.RWMol()
    atom_index = []
    for atom_number in atom_numbers:
        atom = Chem.Atom(atom_number)
        molecular_index = molecule.AddAtom(atom)
        atom_index.append(molecular_index)

    # 在原子和原子直接加入指定种类的键
    for index_x, row_vector in enumerate(bond_matrix):
        for index_y, bond in enumerate(row_vector):
            if index_y <= index_x:
                continue
            if bond == 0:
                continue
            else:
                molecule.AddBond(
                    atom_index[index_x], atom_index[index_y], bond_list[bond]
                )
    return molecule.GetMol()


def bond_pairs_2_mol(
    atom_numbers: List[int], bond_pairs: List[int]
) -> RdMol:
    # TODO 速度有待优化
    molecule = Chem.RWMol()
    atom_index = []
    for atom_number in atom_numbers:
        atom = Chem.Atom(atom_number)
        molecular_index = molecule.AddAtom(atom)
        atom_index.append(molecular_index)
    bond_pairs_ = [bond_pairs[i : i + 3] for i in range(0, len(bond_pairs), 3)]
    for atom_x, atom_y, bond in bond_pairs_:
        if atom_x <= atom_y:
            continue
        if bond == 0:
            continue
        else:
            molecule.AddBond(atom_index[atom_x], atom_index[atom_y], bond_list[bond])
    return molecule.GetMol()


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


def replace_mol(mol, query_smi, replacement_smi):
    query = Chem.MolFromSmarts(query_smi)
    queried_idx = mol.GetSubstructMatch(query)
    for atom in mol.GetAtomWithIdx(queried_idx[0]).GetNeighbors():
        if atom.GetIdx() not in mol.GetSubstructMatch(Chem.MolFromSmarts(query_smi)):
            start = atom.GetIdx()
            break
    # print(start)
    geometry.standard_orient(mol, [start, queried_idx[0]])
    replacement = Chem.MolFromSmiles(replacement_smi)
    replacement = Chem.AddHs(replacement)
    AllChem.EmbedMolecule(replacement)
    replacement.GetAtomWithIdx(0).SetNumRadicalElectrons(0)
    rr = fix_geometry(replacement)
    new_mol = Chem.ReplaceSubstructs(
        mol,
        query,
        rr,
    )[0]
    geometry.standard_orient(new_mol, [0, 1, 2])
    return new_mol


def fix_geometry(replacement):
    r = Chem.RWMol(replacement)
    idx = r.AddAtom(Chem.Atom("C"))
    r.AddBond(0, idx, Chem.BondType.SINGLE)
    r.UpdatePropertyCache()
    AllChem.EmbedMolecule(r)
    geometry.standard_orient(r, [idx, 0])
    r.RemoveAtom(idx)
    return r


def run_replacement(
    data,
    templates,
    replacement_mapping,
    structure_path,
    attempt_time=5,
    overwrite=False,
) -> dict:
    smiles_dict = {}
    with tqdm(total=len(data) * 4) as pbar:
        for template_key, template_value in templates.items():
            template_mol = Chem.MolFromMolFile(
                os.path.join(structure_path, template_value["path"]),
                sanitize=False,
                removeHs=False,
                strictParsing=False,
            )

            for data_idx, data_row in data.iterrows():
                complex_path = os.path.join(
                    structure_path,
                    f"{template_value['output_path']}{os.path.sep}{data_row['Unnamed: 0']}-{template_value['reaction_type']}.xyz",
                )
                mol = deepcopy(template_mol)
                max_volume = 0
                flag = True
                for _ in range(attempt_time):
                    temp_mol = deepcopy(template_mol)
                    flag, temp_mol = attempt_replacement(
                        replacement_mapping, template_value, data_row, temp_mol, flag
                    )
                    if not flag:
                        break
                    temp_volume = AllChem.ComputeMolVolume(temp_mol)
                    if temp_volume > max_volume:
                        mol = deepcopy(temp_mol)
                        max_volume = temp_volume

                    if os.path.exists(complex_path) and not overwrite:
                        break
                if not flag:
                    continue
                for key, fragment in template_value["fragments"].items():
                    idx = data[key].drop_duplicates().to_list().index(data_row[key])
                    cut_out_fragments(
                        mol,
                        prefix=fragment["path_prefix"],
                        smarts=fragment["smarts"],
                        idx=idx,
                        reaction_type=template_value["reaction_type"],
                        output_path=os.path.join(
                            structure_path, fragment["output_path"]
                        ),
                        smiles_dict=smiles_dict,
                        overwrite=overwrite,
                    )
                smiles_dict[complex_path] = Chem.MolToSmiles(mol)
                save_mol(
                    mol, ".".join(complex_path.split(".")[:-1]), overwrite=overwrite
                )
                pbar.set_postfix_str(
                    s=f"{data_row['Unnamed: 0']}-{template_value['reaction_type']}"
                )
                pbar.update()
        pbar.update()
    return smiles_dict


def attempt_replacement(replacement_mapping, template_value, data_row, temp_mol, flag):
    for key, value in template_value.items():
        if key in ("path", "output_path", "reaction_type", "fragments"):
            continue
        if data_row[key] not in replacement_mapping[key][value]:
            flag = False
            break
        methods = replacement_mapping[key][value][data_row[key]]
        for method in methods:
            temp_mol = replace_mol(
                temp_mol, method["query_smi"], method["replacement_smi"]
            )

    return flag, temp_mol


def cut_out_fragments(
    mol, prefix, smarts, idx, reaction_type, output_path, smiles_dict, overwrite=False
):
    rxn = AllChem.ReactionFromSmarts(smarts)
    try:
        m = rxn.RunReactants(Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False))[
            0
        ][0]
    except Exception:
        m = rxn.RunReactants([mol])[0][0]

    file_path = os.path.join(output_path, f"{reaction_type}_{prefix}_{idx}")
    smiles_dict[
        os.path.join(file_path, f"{reaction_type}_{prefix}_{idx}.xyz")
    ] = Chem.MolToSmiles(m)
    if not os.path.exists(file_path):
        os.mkdir(file_path)
    file_name = f"{reaction_type}_{prefix}_{idx}"
    save_mol(m, os.path.join(file_path, file_name), overwrite=overwrite)


def save_mol(m, file_path, overwrite=False):
    if (
        os.path.exists(f"{file_path}.sdf")
        and os.path.exists(f"{file_path}.xyz")
        and not overwrite
    ):
        return
    try:
        Chem.MolToMolFile(m, f"{file_path}.sdf")
    except Exception:
        os.remove(f"{file_path}.sdf")
    Chem.MolToXYZFile(m, f"{file_path}.xyz")


def fix_radical(mol):
    radical_idxs = [
        idx for idx, atom in enumerate(mol.GetAtoms()) if atom.GetNumRadicalElectrons()
    ]
    if len(radical_idxs):
        m = Chem.RWMol(mol)
        for radical_idx in radical_idxs:
            m.GetAtomWithIdx(radical_idx).SetNumRadicalElectrons(0)
        m.AddBond(radical_idxs[0], radical_idxs[1], Chem.BondType.SINGLE)
        return m
    return mol


def fix_dative(mol):
    pt = Chem.GetPeriodicTable()
    dative_atom_idxs = [
        idx
        for idx, atom in enumerate(mol.GetAtoms())
        if atom.GetExplicitValence()
        > atom.GetFormalCharge() + pt.GetDefaultValence(atom.GetAtomicNum())
    ]


def fix_mol(mol):
    m = fix_radical(mol)
    m = fix_dative(m)
    return m


"""def load_mol(file_path) -> RdMol:
    # attempt to read sdf first, if false then read xyz
    file_name = ".".join(file_path.split(".")[:-1])
    try:
        mol = Chem.MolFromMolFile(f"{file_name}.sdf", removeHs=False)
    except Exception:
        mol = Chem.MolFromMolFile(
            MolFormatConversion(f"{file_name}.xyz"), removeHs=False
        )
    return fix_radical(mol)"""
