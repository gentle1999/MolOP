"""
Including functions related to the structure of molecules
"""

import itertools
from typing import Generator, List, Optional, Sequence, Tuple, Union

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import rdForceFieldHelpers, rdMolTransforms
from rdkit.Chem.rdDistGeom import EmbedMolecule

from molop.config import moloplogger
from molop.structure import GeometryTransformation
from molop.structure.utils import (
    bond_list,
    bond_stereo_list,
    bond_stereo_mapping,
    bond_type_mapping,
    estimate_bond_length,
    pt,
)
from molop.utils.functions import is_metal
from molop.utils.types import RdMol

DEBUG_TAG = "[STRUCTURE]"

RDLogger.DisableLog("rdApp.*")  # type: ignore


def get_bond_pairs(mol: RdMol) -> List[Tuple[int, int, int, int]]:
    """
    Get the list of bond pairs in the molecule.

    Parameters:
        mol (RdMol):
            The input molecule
    Returns:
        (List[Tuple[int, int, int, int]]):
            The list of bond pairs in the molecule. Each bond pair is represented
            by a tuple of two atom indices and a bond type index and a bond stereo type index.
    """
    return [
        (
            bond.GetBeginAtomIdx(),
            bond.GetEndAtomIdx(),
            bond_list.index(bond.GetBondType()),
            bond_stereo_list.index(bond.GetStereo()),
        )
        for bond in mol.GetBonds()
    ]


def get_formal_charges(mol: RdMol) -> List[int]:
    """
    Get formal charges of mol.

    Parameters:
        mol (RdMol):
            The input molecule

    Returns:
        List[int]:
            The list of formal charges of atoms in the molecule.
    """
    return [atom.GetFormalCharge() for atom in mol.GetAtoms()]


def get_formal_num_radicals(mol: RdMol) -> List[int]:
    """
    Get formal radical numbers of mol.

    Parameters:
        mol (RdMol):
            The input molecule

    Returns:
        List[int]:
            The list of formal radical numbers of atoms in the molecule.
    """
    return [atom.GetNumRadicalElectrons() for atom in mol.GetAtoms()]


def get_total_charge(mol: RdMol) -> int:
    """
    Get total charge of mol.

    Parameters:
        mol (RdMol):
            The input molecule

    Returns:
        int:
            The total charge of the molecule.
    """
    return sum(get_formal_charges(mol))


def get_total_num_radical(mol: RdMol) -> int:
    """
    Get total radical number of mol.

    Parameters:
        mol (RdMol):
            The input molecule

    Returns:
        int:
            The total radical number of the molecule.
    """
    return sum(get_formal_num_radicals(mol))


def get_total_multiplicity(mol: RdMol) -> int:
    """
    Get total spin multiplicity of mol.

    Parameters:
        mol (RdMol):
            The input molecule

    Returns:
        int:
            The total spin multiplicity of the molecule.
    """
    return get_total_num_radical(mol) + 1


def get_resonance_structures(rdmol: RdMol, flags=0) -> Generator[RdMol, None, None]:
    """
    Get all resonance structures of rdmol.

    Parameters:
        rdmol (RdMol):
            The input molecule
        flags (int):
            Flags for resonance structure generation.

    Returns:
        (Generator[RdMol, None, None]):
            The list of resonance structures of the input molecule.
    """
    return (mol for mol in Chem.ResonanceMolSupplier(rdmol, flags))


def structure_score(rdmol: RdMol) -> int:
    """
    Calculate the structure score of rdmol.

    Parameters:
        rdmol (RdMol):
            The input molecule

    Returns:
        int:
            The structure score of the input molecule.
    """
    return sum(atom.GetNumRadicalElectrons() for atom in rdmol.GetAtoms()) + sum(
        abs(atom.GetFormalCharge()) for atom in rdmol.GetAtoms()
    )


def check_mol_equal(rdmol_1: RdMol, rdmol_2: RdMol) -> bool:
    """
    Check if two molecules are equal. Two molecules are equal if the molecule A contains molecule B
    while the molecule B contains molecule A.

    Parameters:
        rdmol_1 (RdMol):
            The first input molecule
        rdmol_2 (RdMol):
            The second input molecule

    Returns:
        bool:
            True if the two molecules are equal, False otherwise.
    """
    return rdmol_1.HasSubstructMatch(rdmol_2) and rdmol_2.HasSubstructMatch(rdmol_1)


def build_mol_from_atoms_and_bonds(
    atoms: Sequence[int | str],
    bonds: Sequence[Tuple[int, int, int, int]],
    formal_charges: Optional[Sequence[int]] = None,
    formal_num_radicals: Optional[Sequence[int]] = None,
    *,
    coords: Optional[Sequence[Tuple[float, float, float]]] = None,
) -> RdMol:
    """
    Build a molecule from atoms and bonds.

    Parameters:
        atoms (Sequence[int]):
            The list of atom indices.
        bonds (Sequence[Tuple[int, int, int, int]]):
            The list of bond pairs. Each bond pair is represented by a tuple of two atom indices
            and a bond type index and a bond stereo type index.
        formal_charges (Optional[Sequence[int]]):
            The list of formal charges of atoms. If None, all atoms are assigned 0 formal charge.
        formal_num_radicals (Optional[Sequence[int]]):
            The list of formal radical numbers of atoms. If None, all atoms are assigned 0 radical.
        coords (Optional[Sequence[Tuple[float, float, float]]]):
            The list of coordinates of atoms. If None, return a molecule without coordinates.

    Returns:
        RdMol:
            The molecule built from the input atoms and bonds.
    """
    if formal_charges is None:
        formal_charges = [0] * len(atoms)
    if formal_num_radicals is None:
        formal_num_radicals = [0] * len(atoms)
    if coords is None:
        mol = Chem.RWMol()
        for element, formal_charge, formal_num_radical in zip(
            atoms, formal_charges, formal_num_radicals, strict=True
        ):
            atom = Chem.Atom(
                element if isinstance(element, str) else pt.GetElementSymbol(element)
            )
            atom.SetFormalCharge(formal_charge)
            atom.SetNumRadicalElectrons(formal_num_radical)
            mol.AddAtom(atom)
    else:
        temp_xyz = (
            f"{len(atoms)}\n\n"
            + "\n".join(
                f"{element if isinstance(element, str) else pt.GetElementSymbol(element)} {x:15.6f} {y:15.6f} {z:15.6f}"
                for element, (x, y, z) in zip(atoms, coords, strict=True)
            )
            + "\n"
        )
        mol = Chem.RWMol(Chem.MolFromXYZBlock(temp_xyz))
        if mol is None:
            raise ValueError("Invalid input coordinates.")
        for atom_idx, (formal_charge, formal_num_radical) in enumerate(
            zip(formal_charges, formal_num_radicals, strict=True)
        ):
            atom = mol.GetAtomWithIdx(atom_idx)
            atom.SetFormalCharge(formal_charge)
            atom.SetNumRadicalElectrons(formal_num_radical)
    for bond_indexes in bonds:
        begin_atom_idx, end_atom_idx, bond_type_idx, bond_stereo_idx = bond_indexes
        bond_type = bond_list[bond_type_idx]
        bond_stereo = bond_stereo_list[bond_stereo_idx]
        mol.AddBond(
            beginAtomIdx=begin_atom_idx, endAtomIdx=end_atom_idx, order=bond_type
        )
    Chem.SanitizeMol(mol, catchErrors=True)
    Chem.SetAromaticity(mol)
    Chem.DetectBondStereochemistry(mol)
    Chem.SetBondStereoFromDirections(mol)
    Chem.AssignStereochemistryFrom3D(mol)
    Chem.AssignCIPLabels(mol)
    for bond_indexes in bonds:
        begin_atom_idx, end_atom_idx, bond_type_idx, bond_stereo_idx = bond_indexes
        mol.GetBondBetweenAtoms(begin_atom_idx, end_atom_idx).SetStereo(bond_stereo)
    return mol.GetMol()


def get_sub_mol(origin_mol: RdMol, scale: List[int]) -> RdMol:
    """
    Get a sub-molecule of origin_mol with the given scale.

    Parameters:
        origin_mol (RdMol):
            The input molecule
        scale (List[int]):
            The list of atom indices to be included in the sub-molecule.

    Returns:
        RdMol:
            The sub-molecule of origin_mol with the given scale.
    """
    atoms = [origin_mol.GetAtomWithIdx(idx).GetAtomicNum() for idx in scale]
    formal_charges = [origin_mol.GetAtomWithIdx(idx).GetFormalCharge() for idx in scale]
    formal_num_radicals = [
        origin_mol.GetAtomWithIdx(idx).GetNumRadicalElectrons() for idx in scale
    ]
    bonds = [
        (
            bond.GetBeginAtomIdx(),
            bond.GetEndAtomIdx(),
            bond_list[bond.GetBondType()],
            bond_stereo_list[bond.GetStereo()],
        )
        for bond in origin_mol.GetBonds()
        if bond.GetBeginAtomIdx() in scale and bond.GetEndAtomIdx() in scale
    ]
    if origin_mol.GetNumConformers() > 0:
        conf = origin_mol.GetConformer()
        coords = [
            (
                conf.GetAtomPosition(idx).x,
                conf.GetAtomPosition(idx).y,
                conf.GetAtomPosition(idx).z,
            )
            for idx in scale
        ]
    else:
        coords = None
    return build_mol_from_atoms_and_bonds(
        atoms, bonds, formal_charges, formal_num_radicals, coords=coords
    )


def transform_replacement_index(
    mol: RdMol,
    bond_tag: Chem.rdchem.BondType,
    *,
    relative_idx: int = 0,
    absolute_idx: Optional[int] = None,
) -> RdMol:
    """
    Transform the index of the atom in the molecule to let the first atom to be radical atom.

    Parameters:
        mol (RdMol):
            The input molecule
        bond_tag (Chem.rdchem.BondType):
            The bond type of the bond to be replaced.
        relative_idx (int):
            The relative index of the radical atom in the molecule to be transformed to the first atom.
        absolute_idx (Union[int, None]):
            Priority is higher than relative_idx.
            The absolute index of the radical atom in the molecule to be transformed to the first atom.
            If None, the function will try to find the first atom in the molecule that is
            a radical atom.

    Returns:
        RdMol: The transformed molecule.
    """
    if bond_tag in (Chem.rdchem.BondType.SINGLE,):
        radical_idxs: List[int] = [
            atom.GetIdx()
            for atom in mol.GetAtoms()
            if atom.GetNumRadicalElectrons() >= 1
        ]
        if absolute_idx is None:
            if relative_idx >= len(radical_idxs):
                raise ValueError("Relative index is out of range.")
            absolute_idx = radical_idxs[relative_idx]
        elif absolute_idx not in radical_idxs:
            raise ValueError("Absolute index is not a radical atom.")
        return reset_atom_index(mol, [absolute_idx])
    elif bond_tag in (
        Chem.rdchem.BondType.DATIVE,
        Chem.rdchem.BondType.DATIVEL,
        Chem.rdchem.BondType.DATIVER,
        Chem.rdchem.BondType.DATIVEONE,
    ):
        idxs: List[int] = [atom.GetIdx() for atom in mol.GetAtoms()]
        if absolute_idx is None:
            if relative_idx >= len(idxs):
                raise ValueError("Relative index is out of range.")
            absolute_idx = idxs[relative_idx]
        elif absolute_idx not in idxs:
            raise ValueError("Absolute index is not a legal atom.")
        return reset_atom_index(mol, [absolute_idx])
    elif bond_tag in (Chem.rdchem.BondType.DOUBLE,):
        radical_idxs = [
            atom.GetIdx()
            for atom in mol.GetAtoms()
            if atom.GetNumRadicalElectrons() >= 2
        ]
        if absolute_idx is None:
            if relative_idx >= len(radical_idxs):
                raise ValueError("Relative index is out of range.")
            absolute_idx = radical_idxs[relative_idx]
        elif absolute_idx not in radical_idxs:
            raise ValueError("Absolute index is not a radical atom.")
        return reset_atom_index(mol, [absolute_idx])
    elif bond_tag in (Chem.rdchem.BondType.TRIPLE,):
        radical_idxs = [
            atom.GetIdx()
            for atom in mol.GetAtoms()
            if atom.GetNumRadicalElectrons() >= 3
        ]
        if absolute_idx is None:
            if relative_idx >= len(radical_idxs):
                raise ValueError("Relative index is out of range.")
            absolute_idx = radical_idxs[relative_idx]
        elif absolute_idx not in radical_idxs:
            raise ValueError("Absolute index is not a radical atom.")
        return reset_atom_index(mol, [absolute_idx])
    else:
        raise ValueError(f"Unsupported bond type: {bond_tag}.")


def replace_mol(
    mol: RdMol,
    replacement_mol: RdMol,
    query_mol: Union[RdMol, None] = None,
    bind_idx: Union[int, None] = None,
    threshold=0.75,
    angle_split=10,
    randomSeed=114514,
    start_idx: Union[int, None] = None,
    end_idx: Union[int, None] = None,
    *,
    replacement_relative_idx: int = 0,
    replacement_absolute_idx: Union[int, None] = None,
    prefer_ZE: str = "Z",
) -> RdMol:
    """
    Replace the query_mol with replacement_mol in mol.

    Parameters:
        mol (RdMol):
            The input molecule
        replacement_mol (RdMol):
            The replacement molecule
        query_mol (Union[RdMol, None]):
            The query molecule to be replaced.
        bind_idx (int):
            The index of the atom in query_mol to be bound to the replacement_mol.
            If None, the function will try to find the first atom in query_mol that
            is not in the same fragment as the end atom of the bond to be replaced.
        threshold (float):
            The threshold for the similarity between the query_mol and the substructure
            match of replacement_mol in mol.
        angle_split (int):
            The number of angles to split the replacement_mol.
        randomSeed (int):
            The random seed for the conformer generation.
        start_idx (int):
            The index of the atom in the main structure to be replaced. If None, the
            function will try to find the first atom in the main structure that is
            bonded to the end atom of the bond to be replaced.
        end_idx (int):
            The index of the atom in the main structure that is bonded to the start
            atom of the bond to be replaced. If None, the function will try to find
            the first atom in the main structure that is bonded to the start atom of
            the bond to be replaced.
        replacement_relative_idx (int):
            The relative index of the radical atom in the replacement molecule to be
            transformed to the first atom.
        replacement_absolute_idx (Union[int, None]):
            Priority is higher than replacement_relative_idx.
            The absolute index of the radical atom in the replacement molecule to be
            transformed to the first atom.
            If None, the function will try to find the first atom in the replacement
            molecule that is a radical atom.
        prefer_ZE (str):
            The preferred stereochemistry of the bond to be replaced.
            If "Z", the function will try to replace the bond with Z stereochemistry.
            If "E", the function will try to replace the bond with E stereochemistry.
            only works for bond type DOUBLE.

    Returns:
        RdMol:
            The new molecule with the replacement.
    """
    if start_idx and end_idx:
        if (
            mol.GetBondBetweenAtoms(start_idx, end_idx)
            and not mol.GetBondBetweenAtoms(start_idx, end_idx).IsInRing()
        ):
            end = end_idx
            start = start_idx
            bond_tag = mol.GetBondBetweenAtoms(start, end).GetBondType()
        else:
            raise ValueError(
                "start_idx and end_idx should be bonded but not in a ring."
            )
    else:
        if query_mol is None:
            raise ValueError("start_idx and end_idx or query_mol should be provided.")
        else:
            query_mol_ = Chem.MolFromSmarts(Chem.MolToSmarts(query_mol))
            queried_idx_list = mol.GetSubstructMatches(query_mol_, uniquify=False)
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
            for atom in mol.GetAtomWithIdx(unique_queried_idx[0]).GetNeighbors():
                if atom.GetIdx() not in unique_queried_idx:
                    start = atom.GetIdx()
                    end = unique_queried_idx[0]
                    bond_tag = mol.GetBondBetweenAtoms(start, end).GetBondType()
                    break
    assert any(atom.GetNumRadicalElectrons() for atom in replacement_mol.GetAtoms()), (
        "replacement_mol should have at least one radical atom."
    )

    moloplogger.debug(
        f"{DEBUG_TAG} | Initial check of structure replacement passed,"
        f" start: {start} and end: {end} atom found. Bond tag: {bond_tag}."
    )

    # start: the idx of atom that bind to the substruct to be replaced
    # end: the idx of atom in substruct to be replaced, and bind to the main structure
    origin_mol = Chem.RWMol(mol)

    if bond_tag in (
        Chem.rdchem.BondType.SINGLE,
        Chem.rdchem.BondType.DATIVE,
        Chem.rdchem.BondType.DATIVEL,
        Chem.rdchem.BondType.DATIVER,
        Chem.rdchem.BondType.DATIVEONE,
    ):
        GeometryTransformation.standard_orient(origin_mol, [start, end])
        # build skeleton with queried substruct removed
        # considering the original mol may have been split into frags
        skeleton, new_start = get_skeleton(origin_mol=origin_mol, start=start, end=end)

        # build replacement conformer
        replacement = build_replacement(
            replacement_mol=replacement_mol,
            bond_tag=bond_tag,
            randomSeed=randomSeed,
            replacement_relative_idx=replacement_relative_idx,
            replacement_absolute_idx=replacement_absolute_idx,
            start_atomicnum=mol.GetAtomWithIdx(start).GetAtomicNum(),
        )
        # combine skeleton and replacement
        new_start, rmol = combine_skeleton_replacement(
            skeleton=skeleton,
            replacement=replacement,
            start=new_start,
            bond_tag=bond_tag,
        )
        set_best_dihedral(
            rmol=rmol,
            threshold=threshold,
            angle_split=angle_split,
            start=new_start,
            end=0,
        )
    elif bond_tag in (Chem.rdchem.BondType.DOUBLE,):
        # rotate the double bond to put on a plane (if necessary)
        for atom in origin_mol.GetAtomWithIdx(start).GetNeighbors():
            if atom.GetIdx() != end:
                neighbor_idx = atom.GetIdx()
                GeometryTransformation.standard_orient(
                    origin_mol, [start, end, neighbor_idx]
                )
                break
        else:
            GeometryTransformation.standard_orient(origin_mol, [start, end])
        skeleton, new_start = get_skeleton(origin_mol=origin_mol, start=start, end=end)
        replacement = build_replacement(
            replacement_mol=replacement_mol,
            bond_tag=bond_tag,
            randomSeed=randomSeed,
            replacement_relative_idx=replacement_relative_idx,
            replacement_absolute_idx=replacement_absolute_idx,
            start_atomicnum=mol.GetAtomWithIdx(start).GetAtomicNum(),
        )
        anchors = [
            atom.GetIdx()
            for atom in replacement.GetAtomWithIdx(0).GetNeighbors()
            if atom.GetIdx() != 0
        ]
        replacements = []
        if len(anchors) == 0:
            replacement_ = Chem.RWMol(replacement)
            idx = replacement_.AddAtom(Chem.Atom(1))
            GeometryTransformation.standard_orient(replacement_, [idx, 0])
            replacement_.RemoveAtom(idx)
            replacements.append(replacement_)
        else:
            for anchor in anchors:
                replacement_ = Chem.RWMol(replacement)
                idx = replacement_.AddAtom(Chem.Atom(1))
                GeometryTransformation.standard_orient(replacement_, [idx, 0, anchor])
                replacement_.RemoveAtom(idx)
                replacements.append(replacement_)
        if len(replacements) == 0:
            raise ValueError("No valid replacement found.")
        for replacement in replacements:
            new_start, rmol = combine_skeleton_replacement(
                skeleton=skeleton,
                replacement=replacement,
                start=new_start,
                bond_tag=bond_tag,
            )
            temp_mol = Chem.MolFromMolBlock(
                Chem.MolToV3KMolBlock(rmol, includeStereo=False), removeHs=False
            )
            Chem.SanitizeMol(temp_mol)
            Chem.DetectBondStereochemistry(temp_mol)
            if (
                temp_mol.GetBondBetweenAtoms(new_start, 0).GetStereo()
                == bond_stereo_mapping[prefer_ZE]
            ):
                break
        for anchor in rmol.GetAtomWithIdx(0).GetNeighbors():
            if anchor.GetIdx() != new_start:
                set_best_dihedral(
                    rmol=rmol,
                    threshold=threshold,
                    angle_split=angle_split,
                    start=0,
                    end=anchor.GetIdx(),
                )
    elif bond_tag in (Chem.rdchem.BondType.TRIPLE,):
        GeometryTransformation.standard_orient(origin_mol, [start, end])
        skeleton, new_start = get_skeleton(origin_mol=origin_mol, start=start, end=end)
        # build replacement conformer
        replacement = build_replacement(
            replacement_mol=replacement_mol,
            bond_tag=bond_tag,
            randomSeed=randomSeed,
            replacement_relative_idx=replacement_relative_idx,
            replacement_absolute_idx=replacement_absolute_idx,
            start_atomicnum=mol.GetAtomWithIdx(start).GetAtomicNum(),
        )
        # combine skeleton and replacement
        new_start, rmol = combine_skeleton_replacement(
            skeleton=skeleton,
            replacement=replacement,
            start=new_start,
            bond_tag=bond_tag,
        )
        for anchor in rmol.GetAtomWithIdx(0).GetNeighbors():
            if anchor.GetIdx() != new_start:
                set_best_dihedral(
                    rmol=rmol,
                    threshold=threshold,
                    angle_split=angle_split,
                    start=0,
                    end=anchor.GetIdx(),
                )
    else:
        raise ValueError(f"Unsupported bond type: {bond_tag}.")
    GeometryTransformation.standard_orient(rmol, [0, 1])
    # recover the original atom index of the skeleton
    atom_idxs = list(range(rmol.GetNumAtoms()))
    mapping = (
        atom_idxs[replacement.GetNumAtoms() :] + atom_idxs[: replacement.GetNumAtoms()]
    )
    return make_dative_bonds(reset_atom_index(rmol, mapping))


def set_best_dihedral(
    rmol: Chem.RWMol,
    threshold: float,
    angle_split: int,
    start: int,
    end: int,
):
    assert rmol.GetBondBetweenAtoms(start, end) is not None, "The bond is not found."
    assert rmol.GetBondBetweenAtoms(start, end).GetBondType() in (
        Chem.rdchem.BondType.SINGLE,
        Chem.rdchem.BondType.DATIVE,
        Chem.rdchem.BondType.DATIVEL,
        Chem.rdchem.BondType.DATIVER,
        Chem.rdchem.BondType.DATIVEONE,
    ), "The selected bond type should be rotatable."
    first_atom_idx, forth_atom_idx = None, None
    for first_atom in rmol.GetAtomWithIdx(start).GetNeighbors():
        if first_atom.GetIdx() != end:
            first_atom_idx = first_atom.GetIdx()
            break

    for forth_atom in rmol.GetAtomWithIdx(end).GetNeighbors():
        if forth_atom.GetIdx() != start:
            forth_atom_idx = forth_atom.GetIdx()
            break
    if first_atom_idx and forth_atom_idx:
        moloplogger.debug(
            f"{DEBUG_TAG} | Rotating the bond between {start} and {end} to the best angle."
        )
        Chem.SanitizeMol(rmol)
        best_angle, best_crowding_score = 0, float("-inf")
        for angle in np.linspace(0, 2 * np.pi, angle_split):
            rdMolTransforms.SetDihedralRad(
                rmol.GetConformer(),
                first_atom_idx,
                start,
                end,
                forth_atom_idx,
                angle,
            )
            crowding_score = get_crowding_socre(rmol)
            if crowding_score > best_crowding_score and check_crowding(rmol, threshold):
                best_angle = angle
                best_crowding_score = crowding_score
                moloplogger.debug(
                    f"{DEBUG_TAG} | The best angle is {best_angle} with crowding score {best_crowding_score}."
                )
        rdMolTransforms.SetDihedralRad(
            rmol.GetConformer(),
            first_atom_idx,
            start,
            end,
            forth_atom_idx,
            best_angle,
        )


def combine_skeleton_replacement(
    skeleton: Union[Chem.RWMol, RdMol],
    replacement: Union[Chem.RWMol, RdMol],
    start: int,
    bond_tag: Chem.rdchem.BondType,
):
    new_mol = Chem.CombineMols(replacement, skeleton)
    # find the new index mapping of start and end (end in the replacement now)
    new_start = start + replacement.GetNumAtoms()
    rmol = Chem.RWMol(new_mol)
    rmol.AddBond(0, new_start, bond_tag)
    rmol.GetAtomWithIdx(new_start).SetNumRadicalElectrons(
        max(
            0,
            rmol.GetAtomWithIdx(new_start).GetNumRadicalElectrons()
            - bond_type_mapping[bond_tag],
        )
    )
    rmol.GetAtomWithIdx(0).SetNumRadicalElectrons(
        max(
            0,
            rmol.GetAtomWithIdx(0).GetNumRadicalElectrons()
            - bond_type_mapping[bond_tag],
        )
    )
    moloplogger.debug(
        f"{DEBUG_TAG} | Skeleton and Replacement combination of structure replacement passed"
    )

    return new_start, rmol


def get_skeleton(origin_mol: Chem.RWMol, start: int, end: int):
    # Avoiding Aromaticity Determination Problems Associated with Hetero-Hydrogen Bond
    # Substitution in Heterocyclic Aromatic Hydrocarbons
    Chem.Kekulize(origin_mol)
    bond_type = origin_mol.GetBondBetweenAtoms(start, end).GetBondType()
    original_frags = Chem.GetMolFrags(origin_mol, asMols=True)
    origin_mol.RemoveBond(start, end)
    origin_mol.GetAtomWithIdx(start).SetNumRadicalElectrons(
        origin_mol.GetAtomWithIdx(start).GetNumRadicalElectrons()
        + bond_type_mapping[bond_type]
    )
    origin_mol.GetAtomWithIdx(end).SetNumRadicalElectrons(
        origin_mol.GetAtomWithIdx(end).GetNumRadicalElectrons()
        + bond_type_mapping[bond_type]
    )
    # get all frags
    frags, frag_idxs = (
        Chem.GetMolFrags(origin_mol, asMols=True),
        Chem.GetMolFrags(origin_mol),
    )
    # check if the query_mol is in a ring
    if len(frags) <= len(original_frags):
        raise RuntimeError("query_mol should not be in a ring")
    ske_frags = [
        (frag, frag_idx)
        for frag, frag_idx in zip(frags, frag_idxs, strict=True)
        if end not in frag_idx
    ]  # get all frags except the one that contains the end atom
    skeleton = ske_frags[0][0]
    init_mapping = list(ske_frags[0][1])
    for frag, frag_idx in ske_frags[1:]:
        skeleton = Chem.CombineMols(skeleton, frag)
        init_mapping.extend(frag_idx)
    moloplogger.debug(f"{DEBUG_TAG} | initial mapping: {init_mapping}.")
    new_start = init_mapping.index(start)
    mapping = [x[0] for x in sorted(enumerate(init_mapping), key=lambda x: x[1])]
    skeleton = reset_atom_index(skeleton, mapping)
    moloplogger.debug(
        f"{DEBUG_TAG} | Skeleton initialization of structure replacement passed"
        f" with mapping: {mapping}. SMILES: {Chem.MolToSmiles(skeleton)}."
    )
    return skeleton, new_start


def build_replacement(
    replacement_mol: Chem.Mol,
    bond_tag: Chem.rdchem.BondType,
    randomSeed: int,
    start_atomicnum: int,
    replacement_relative_idx: int = 0,
    replacement_absolute_idx: Optional[int] = None,
):
    replacement = Chem.AddHs(replacement_mol, addCoords=True)
    if replacement.GetNumConformers() == 0:
        EmbedMolecule(replacement, randomSeed=randomSeed)
    replacement = transform_replacement_index(
        replacement,
        bond_tag=bond_tag,
        relative_idx=replacement_relative_idx,
        absolute_idx=replacement_absolute_idx,
    )
    replacement.GetAtomWithIdx(0).SetNumRadicalElectrons(0)
    rr = fix_geometry(
        replacement,
        bond_tag=bond_tag,
        bind_type=start_atomicnum,
        randomSeed=randomSeed,
    )
    moloplogger.debug(
        f"{DEBUG_TAG} | Replacement initialization of structure replacement passed"
    )

    return rr


def check_crowding(mol: RdMol, threshold=0.6) -> bool:
    """
    Check if the molecule is crowded.

    Parameters:
        mol (RdMol):
            The input molecule.
        threshold (float):
            The threshold of crowding. If too crowded
            `d(a-b) < threshold * (R(a)+R(b))`, return False.

    Returns:
        bool: True if the molecule is not crowded, False otherwise.
    """
    distances = Chem.Get3DDistanceMatrix(mol)
    for start_atom, end_atom in itertools.combinations(mol.GetAtoms(), 2):
        if distances[start_atom.GetIdx()][
            end_atom.GetIdx()
        ] < threshold * estimate_bond_length(
            start_atom.GetAtomicNum(), end_atom.GetAtomicNum()
        ):
            return False
    return True


def get_crowding_socre(mol: RdMol) -> float:
    """
    Calculate the crowding score of the input molecule.

    Parameters:
        mol (RdMol):
            The input molecule.
    Returns:
        float: The crowding score of the input molecule.
    """
    Chem.SanitizeMol(mol)
    return -rdForceFieldHelpers.UFFGetMoleculeForceField(mol).CalcEnergy()


def attempt_replacement(
    mol,
    query: Union[str, RdMol],
    replacement: Union[str, RdMol],
    bind_idx: Optional[int] = None,
    replace_all=False,
    attempt_num=10,
    crowding_threshold=0.75,
    angle_split=10,
    randomSeed=114514,
    start_idx: Optional[int] = None,
    end_idx: Optional[int] = None,
    *,
    replacement_relative_idx: int = 0,
    replacement_absolute_idx: Union[int, None] = None,
    prefer_ZE: str = "Z",
):
    """
    Attempt to replace the query with the replacement in the input molecule.

    Parameters:
        mol (RdMol):
            The input molecule.
        query (str | RdMol):
            The SMARTS or Mol object to query the substituent in the original molecule.
        replacement (str | RdMol):
            The SMARTS or Mol object of new substituent.
        bind_idx (int):
            The index of the atom to bind the new substituent. The default is None, which means
            to replace the first legal atom in original molecule.
            If specified, try to replace the legal substruct where the atom in it. User should
            meke sure the atom is legal.
            Detail example in (Repalce Substituent)[Repalce Substituent]
        replace_all (bool):
            If True, replace all the substituent queried in the original molecule.
        crowding_threshold (float):
            The threshold of crowding. If the new substituent is too crowded
            `d(a-b) < threshold * (R(a)+R(b))`, the substitution will be rejected.
        angle_split (int):
            Decide how many equal parts of 360° you want to divide. The larger the number the finer
            the rotation will be attempted but the slower the calculation will be.
        randomSeed (int):
            The random seed.
        start_idx (int):
            If both `start_idx` and `end_idx` are specified, simply ignore the `query`, break the
            key between `start_idx` and `end_idx` and replace the base group where `end_idx` is located
        end_idx (int):
            If both `start_idx` and `end_idx` are specified, simply ignore the `query`, break the
            key between `start_idx` and `end_idx` and replace the base group where `end_idx` is located
        replacement_relative_idx (int):
            The relative index of the radical atom in the replacement molecule to be
            transformed to the first atom.
        replacement_absolute_idx (Union[int, None]):
            Priority is higher than replacement_relative_idx.
            The absolute index of the radical atom in the replacement molecule to be
            transformed to the first atom.
            If None, the function will try to find the first atom in the replacement
            molecule that is a radical atom.
        prefer_ZE (str):
            The preferred stereochemistry of the bond to be replaced.
            If "Z", the function will try to replace the bond with Z stereochemistry.
            If "E", the function will try to replace the bond with E stereochemistry.
            only works for bond type DOUBLE.

    Returns:
        RdMol: The new molecule with the replacement.
    """
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
                replacement_relative_idx=replacement_relative_idx,
                replacement_absolute_idx=replacement_absolute_idx,
                prefer_ZE=prefer_ZE,
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
                    replacement_relative_idx=replacement_relative_idx,
                    replacement_absolute_idx=replacement_absolute_idx,
                    prefer_ZE=prefer_ZE,
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
                replacement_relative_idx=replacement_relative_idx,
                replacement_absolute_idx=replacement_absolute_idx,
                prefer_ZE=prefer_ZE,
            )
        if check_crowding(new_mol, threshold=crowding_threshold):
            Chem.SanitizeMol(new_mol)
            return new_mol
    else:
        raise RuntimeError(
            f"replacement '{Chem.MolToSmiles(replacement_mol)}' is too big in threshold {crowding_threshold}, try smaller threshold"
        )


def fix_geometry(
    replacement: RdMol,
    bond_tag: Chem.rdchem.BondType,
    bind_type: int,
    randomSeed=114514,
):
    """
    Fix the geometry of the replacement molecule with the given bind_type.

    Parameters:
        replacement (RdMol):
            The input replacement molecule.
        bind_type (int):
            The atomic number of the atom that binds to the replacement.
        randomSeed (int):
            The random seed for the embedding.

    Returns:
        RdMol: The fixed replacement molecule.
    """
    r = Chem.RWMol(replacement)
    if bond_tag in (Chem.rdchem.BondType.SINGLE,):
        idx = r.AddAtom(Chem.Atom("At"))
        r.AddBond(0, idx, Chem.BondType.SINGLE)
    elif bond_tag in (Chem.rdchem.BondType.DOUBLE,):
        idx = r.AddAtom(Chem.Atom("Po"))
        r.AddBond(0, idx, Chem.BondType.DOUBLE)
    elif bond_tag in (Chem.rdchem.BondType.TRIPLE,):
        idx = r.AddAtom(Chem.Atom("Bi"))
        r.AddBond(0, idx, Chem.BondType.TRIPLE)
    elif bond_tag in (Chem.rdchem.BondType.DATIVE,):
        idx = r.AddAtom(Chem.Atom("Rf"))
        r.AddBond(0, idx, Chem.BondType.DATIVE)
    else:
        raise ValueError(f"Unsupported bond type: {bond_tag}")
    r.UpdatePropertyCache()
    Chem.SanitizeMol(r)
    EmbedMolecule(r, randomSeed=randomSeed)
    GeometryTransformation.standard_orient(r, [idx, 0])
    rdMolTransforms.SetBondLength(
        r.GetConformer(),
        idx,
        0,
        estimate_bond_length(r.GetAtomWithIdx(0).GetAtomicNum(), bind_type, bond_tag),
    )
    r.RemoveAtom(idx)
    return r


def reset_atom_index(mol: Chem.rdchem.Mol, mapping: Sequence[int]) -> Chem.rdchem.Mol:
    """
    Reset the atom index of a molecule according to the given mapping.
    The mapping should be a list of the new index of each atom in the same order as the original molecule.

    e.g. mapping = [2, 0, 1] means the first atom in the new molecule is mapped to the third atom in the original molecule,
    the second atom is mapped to the first atom, and the third atom is mapped to the second atom.

    Parameters:
        mol (Chem.rdchem.Mol): The input molecule.
        mapping (Sequence[int]): The mapping of the atom index.
    Returns:
        Chem.rdchem.Mol: The molecule with the new atom index.
    """
    # new index
    new_idx = list(mapping) + [
        idx for idx, atom in enumerate(mol.GetAtoms()) if idx not in mapping
    ]
    return Chem.RenumberAtoms(mol, new_idx)


def make_dative_bonds(rwmol: Chem.rdchem.RWMol, ratio=1.3) -> Chem.rdchem.RWMol:
    """
    Make dative bonds between the metal atoms and the non-metal atoms.
    Parameters:
        rwmol (Chem.rdchem.RWMol): The editable rdkit molecule.

    Returns:
        Chem.rdchem.RWMol: The editable rdkit molecule with dative bonds.
    """
    Chem.Kekulize(rwmol, clearAromaticFlags=True)
    for atom in rwmol.GetAtoms():
        atom.SetNoImplicit(True)
    exp_ratio = {
        "N": 1.45,
        "O": ratio,
        "S": ratio,
        "P": ratio,
    }
    metal_atoms = [
        atom.GetIdx() for atom in rwmol.GetAtoms() if is_metal(atom.GetAtomicNum())
    ]
    datived_rwmol = []
    for shuffled_metal_atoms in itertools.permutations(metal_atoms, len(metal_atoms)):
        moloplogger.debug(f"{DEBUG_TAG} | Metal atoms order: {shuffled_metal_atoms}")
        temp_rwmol = Chem.RWMol(rwmol)

        # set covalent bonds between metal atoms and negative atoms
        for metal_atom in shuffled_metal_atoms:
            # find negative atoms
            # not BX4-
            # distance <= ratio * estimate_bond_length
            negative_atoms = [
                atom.GetIdx()
                for atom in temp_rwmol.GetAtoms()
                if not (atom.GetAtomicNum() == 5 and atom.GetTotalValence() == 4)
                and atom.GetFormalCharge() < 0
                and atom.GetIdx() != metal_atom
                and temp_rwmol.GetConformer()
                .GetAtomPosition(atom.GetIdx())
                .Distance(temp_rwmol.GetConformer().GetAtomPosition(metal_atom))
                <= ratio
                * estimate_bond_length(
                    temp_rwmol.GetAtomWithIdx(metal_atom).GetAtomicNum(),
                    temp_rwmol.GetAtomWithIdx(atom.GetIdx()).GetAtomicNum(),
                    bond_list[
                        abs(temp_rwmol.GetAtomWithIdx(atom.GetIdx()).GetFormalCharge())
                    ],
                )
            ]
            # sort negative atoms by distance to metal atom
            negative_atoms.sort(
                key=lambda x: temp_rwmol.GetConformer()
                .GetAtomPosition(x)
                .Distance(temp_rwmol.GetConformer().GetAtomPosition(metal_atom))
            )
            moloplogger.debug(f"{DEBUG_TAG} | Negative atoms: {negative_atoms}")
            # add covalent bonds between metal atom and negative atoms
            # until the metal atom has no more positive charge or no more negative atoms
            while (
                # temp_rwmol.GetAtomWithIdx(metal_atom).GetFormalCharge() > 0
                # and negative_atoms
                negative_atoms
            ):
                negative_atom = negative_atoms.pop(0)
                if temp_rwmol.GetBondBetweenAtoms(negative_atom, metal_atom) is None:
                    moloplogger.debug(
                        f"{DEBUG_TAG} | Trying to add dative bond between "
                        f"{metal_atom} and {negative_atom}"
                    )
                    temp_rwmol.AddBond(
                        negative_atom,
                        metal_atom,
                        Chem.BondType.DATIVE,
                    )
        for metal_atom in metal_atoms:
            dative_atoms = [
                idxs[0]
                for idxs in temp_rwmol.GetSubstructMatches(
                    Chem.MolFromSmarts(
                        "[#8v2+0,#8v3+0,#16v2+0,#16v3+0,#16v3+1,#7v3+0,#15v3+0]"
                    )
                )
            ]
            moloplogger.debug(f"{DEBUG_TAG} | possible dative atoms: {dative_atoms}")
            for dative_atom in dative_atoms:
                if not temp_rwmol.GetBondBetweenAtoms(metal_atom, dative_atom):
                    if temp_rwmol.GetConformer().GetAtomPosition(metal_atom).Distance(
                        temp_rwmol.GetConformer().GetAtomPosition(dative_atom)
                    ) < exp_ratio[
                        temp_rwmol.GetAtomWithIdx(dative_atom).GetSymbol()
                    ] * estimate_bond_length(
                        temp_rwmol.GetAtomWithIdx(metal_atom).GetAtomicNum(),
                        temp_rwmol.GetAtomWithIdx(dative_atom).GetAtomicNum(),
                        Chem.rdchem.BondType.DATIVE,
                    ):
                        moloplogger.debug(
                            f"{DEBUG_TAG} | Trying to add dative bond between "
                            f"{metal_atom} and {dative_atom}"
                        )
                        neighbour_atoms = [
                            neighbor.GetIdx()
                            for neighbor in temp_rwmol.GetAtomWithIdx(
                                dative_atom
                            ).GetNeighbors()
                        ]
                        temp_rwmol.AddBond(
                            dative_atom, metal_atom, Chem.BondType.DATIVE
                        )
                        # avoid unusual angle
                        for neighbor_atom in neighbour_atoms:
                            if (
                                rdMolTransforms.GetAngleDeg(
                                    rwmol.GetConformer(),
                                    metal_atom,
                                    dative_atom,
                                    neighbor_atom,
                                )
                                < 100
                            ):
                                temp_rwmol.RemoveBond(dative_atom, metal_atom)
                    Chem.SanitizeMol(temp_rwmol)

        # set covalent bonds between metal atoms and negative atoms that are not metals
        remained_negative_atoms = [
            atom.GetIdx()
            for atom in temp_rwmol.GetAtoms()
            if not (atom.GetAtomicNum() == 5 and atom.GetTotalValence() == 4)
            and atom.GetFormalCharge() < 0
            and not is_metal(atom.GetAtomicNum())
        ]
        moloplogger.debug(
            f"{DEBUG_TAG} | Remained negative atoms: {remained_negative_atoms}"
        )
        for remained_negative_atom in remained_negative_atoms:
            metal_atoms = [
                atom.GetIdx()
                for atom in temp_rwmol.GetAtoms()
                if is_metal(atom.GetAtomicNum())
            ]
            if not metal_atoms:
                continue
            metal_atoms.sort(
                key=lambda x: temp_rwmol.GetConformer()
                .GetAtomPosition(x)
                .Distance(
                    temp_rwmol.GetConformer().GetAtomPosition(remained_negative_atom)
                )
            )
            distance = (
                temp_rwmol.GetConformer()
                .GetAtomPosition(metal_atoms[0])
                .Distance(
                    temp_rwmol.GetConformer().GetAtomPosition(remained_negative_atom)
                )
            )
            tolerance = ratio * estimate_bond_length(
                temp_rwmol.GetAtomWithIdx(metal_atoms[0]).GetAtomicNum(),
                temp_rwmol.GetAtomWithIdx(remained_negative_atom).GetAtomicNum(),
                Chem.rdchem.BondType.DATIVE,
            )
            if distance <= tolerance:
                if (
                    temp_rwmol.GetBondBetweenAtoms(
                        metal_atoms[0], remained_negative_atom
                    )
                    is None
                ):
                    moloplogger.debug(
                        f"{DEBUG_TAG} | Trying to add dative bond between "
                        f"{metal_atoms[0]} and {remained_negative_atom}"
                    )
                    temp_rwmol.AddBond(
                        remained_negative_atom,
                        metal_atoms[0],
                        bond_list[
                            -temp_rwmol.GetAtomWithIdx(
                                remained_negative_atom
                            ).GetFormalCharge()
                        ],
                    )
                    temp_rwmol.GetAtomWithIdx(metal_atoms[0]).SetFormalCharge(
                        temp_rwmol.GetAtomWithIdx(metal_atoms[0]).GetFormalCharge()
                        + temp_rwmol.GetAtomWithIdx(
                            remained_negative_atom
                        ).GetFormalCharge()
                    )
                    temp_rwmol.GetAtomWithIdx(remained_negative_atom).SetFormalCharge(0)
        datived_rwmol.append(temp_rwmol)
    datived_rwmol_smiles = [Chem.MolToSmiles(rwmol) for rwmol in datived_rwmol]
    moloplogger.debug(
        f"{DEBUG_TAG} | Possible resonance structures with dative bonds: \n{datived_rwmol_smiles}"
    )
    return sorted(
        datived_rwmol,
        key=lambda x: sum(abs(atom.GetFormalCharge()) for atom in x.GetAtoms()),
    )[0]
