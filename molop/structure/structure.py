"""
Including functions related to the structure of molecules
"""

import itertools
from typing import List, Sequence, Tuple, Union

import numpy as np
from openbabel import openbabel as ob
from openbabel import pybel
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, SanitizeMol, rdMolTransforms
from rdkit.Chem.rdDistGeom import EmbedMolecule

from molop.logger.logger import moloplogger
from molop.structure import geometry
from molop.utils.functions import is_metal
from molop.utils.types import RdMol

DEBUG_TAG = "[STRUCTURE]"

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


def canonical_smiles(smiles: str):
    original_smi = smiles
    viewed_smi = {original_smi: 1}
    while original_smi != (
        canonical_smi := Chem.CanonSmiles(original_smi, useChiral=True)
    ) and (canonical_smi not in viewed_smi or viewed_smi[canonical_smi] < 2):
        original_smi = canonical_smi
        if original_smi not in viewed_smi:
            viewed_smi[original_smi] = 1
        else:
            viewed_smi[original_smi] += 1
    else:
        return original_smi


def get_bond_pairs(mol: RdMol) -> List[Tuple[int]]:
    """
    Get the list of bond pairs in the molecule.

    Parameters:
        mol (RdMol):
            The input molecule
    Returns:
        List[Tuple[int]]:
            The list of bond pairs in the molecule. Each bond pair is represented
            by a tuple of two atom indices and a bond type index.
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


def get_resonance_structures(rdmol: RdMol, flags=0) -> List[RdMol]:
    """
    Get all resonance structures of rdmol.

    Parameters:
        rdmol (RdMol):
            The input molecule
        flags (int):
            Flags for resonance structure generation.

    Returns:
        List[RdMol]:
            The list of resonance structures of the input molecule.
    """
    suppl = Chem.ResonanceMolSupplier(rdmol, flags)
    return [mol for mol in suppl]


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
    if SanitizeMol(sub_mol):
        return sub_mol.GetMol()
    raise ValueError("Sub-molecule is not valid.")


def get_under_bonded_number(atom: ob.OBAtom) -> int:
    """
    Get the number of atoms under the given atom.
    Suppose the atom is not a metal and follows the Octet rate (Suitable for most small organic molecules)

    Parameters:
        atom (ob.OBAtom): The atom to be checked.
    Returns:
        int: The number of atoms under the given atom.
    """
    if atom.IsMetal():
        return 0
    if atom.GetAtomicNum() <= 2:
        return (
            2
            - pt.GetNOuterElecs(atom.GetAtomicNum())
            - atom.GetTotalValence()
            + atom.GetFormalCharge()
        )
    if atom.GetAtomicNum() == 5 and atom.GetTotalValence() < 4:
        return (
            6
            - -pt.GetNOuterElecs(atom.GetAtomicNum())
            - atom.GetTotalValence()
            + atom.GetFormalCharge()
        )
    return (
        8
        - pt.GetNOuterElecs(atom.GetAtomicNum())
        - atom.GetTotalValence()
        + atom.GetFormalCharge()
    )


def get_radical_number(atom: ob.OBAtom) -> int:
    if atom.IsMetal():
        return 0
    under_bonded_number = get_under_bonded_number(atom)
    if under_bonded_number <= 0:
        return 0
    if atom.GetFormalCharge() > 0:
        return under_bonded_number - atom.GetFormalCharge() * 2
    else:
        return under_bonded_number


def transform_replacement_index(
    mol: RdMol,
    *,
    relative_idx: int = 0,
    absolute_idx: Union[int, None] = None,
) -> RdMol:
    """
    Transform the index of the atom in the molecule to let the first atom to be radical atom.

    Parameters:
        mol (RdMol):
            The input molecule
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
    if mol.GetAtomWithIdx(0).GetNumRadicalElectrons() == 1:
        return mol
    radical_idxs = [
        atom.GetIdx() for atom in mol.GetAtoms() if atom.GetNumRadicalElectrons() == 1
    ]
    if absolute_idx is None:
        absolute_idx = radical_idxs[relative_idx]
    elif absolute_idx not in radical_idxs:
        raise ValueError("Absolute index is not a radical atom.")
    return reset_atom_index(mol, [absolute_idx])


def replace_mol(
    mol: RdMol,
    replacement_mol: RdMol,
    query_mol: Union[RdMol, None] = None,
    bind_idx: int = None,
    threshold=0.75,
    angle_split=10,
    randomSeed=114514,
    start_idx: int = None,
    end_idx: int = None,
    *,
    replacement_relative_idx: int = 0,
    replacement_absolute_idx: Union[int, None] = None,
) -> RdMol:
    """
    Replace the query_mol with replacement_mol in mol.

    Parameters:
        mol (RdMol):
            The input molecule
        replacement_mol (RdMol):
            The replacement molecule
        query_mol (Union[RdMol, None]):
            The query molecule to be replaced. If None, the function will try to find
            the first substructure match of replacement_mol in mol.
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
            bonded to the end atom of the bond to be replaced and has a single bond.
        end_idx (int):
            The index of the atom in the main structure that is bonded to the start
            atom of the bond to be replaced. If None, the function will try to find
            the first atom in the main structure that is bonded to the start atom of
            the bond to be replaced and has a single bond.
        replacement_relative_idx (int):
            The relative index of the radical atom in the replacement molecule to be
            transformed to the first atom.
        replacement_absolute_idx (Union[int, None]):
            Priority is higher than replacement_relative_idx.
            The absolute index of the radical atom in the replacement molecule to be
            transformed to the first atom.
            If None, the function will try to find the first atom in the replacement
            molecule that is a radical atom.

    Returns:
        RdMol:
            The new molecule with the replacement.
    """
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
            raise ValueError(
                "start_idx and end_idx should be bonded with single but not in a ring."
            )
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
    moloplogger.debug(f"{DEBUG_TAG}: Initial check of structure replacement passed")

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
    frags, frag_idxs = (
        Chem.GetMolFrags(origin_mol, asMols=True),
        Chem.GetMolFrags(origin_mol),
    )
    ske_frags = [
        (frag, frag_idx)
        for frag, frag_idx in zip(frags, frag_idxs)
        if end not in frag_idx
    ]  # get all frags except the one that contains the end atom
    skeleton = ske_frags[0][0]
    init_mapping = list(ske_frags[0][1])
    for frag, frag_idx in ske_frags[1:]:
        skeleton = Chem.CombineMols(skeleton, frag)
        init_mapping.extend(frag_idx)
    mapping = list(
        map(lambda x: x[0], sorted(list(enumerate(init_mapping)), key=lambda x: x[1]))
    )
    skeleton = reset_atom_index(skeleton, mapping)
    moloplogger.debug(
        f"{DEBUG_TAG}: Skeleton initialization of structure replacement passed"
    )

    # build replacement conformer
    replacement = Chem.AddHs(replacement_mol, addCoords=True)
    if replacement.GetNumConformers() == 0:
        EmbedMolecule(replacement, randomSeed=randomSeed)
    replacement = transform_replacement_index(
        replacement,
        relative_idx=replacement_relative_idx,
        absolute_idx=replacement_absolute_idx,
    )
    replacement.GetAtomWithIdx(0).SetNumRadicalElectrons(0)
    rr = fix_geometry(
        replacement,
        bind_type=mol.GetAtomWithIdx(start).GetAtomicNum(),
        randomSeed=randomSeed,
    )
    moloplogger.debug(
        f"{DEBUG_TAG}: Replacement initialization of structure replacement passed"
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
    moloplogger.debug(
        f"{DEBUG_TAG}: Skeleton and Replacement combination of structure replacement passed"
    )

    # find the idx list of replacement
    idx_list = [atom.GetIdx() for atom in rr.GetAtoms()]
    if replacement.GetNumAtoms() >= 2:
        for forth_atom in rmol.GetAtomWithIdx(0).GetNeighbors():
            if forth_atom.GetIdx() in idx_list:
                forth_atom_idx = forth_atom.GetIdx()
                break
        for first_atom in rmol.GetAtomWithIdx(lines_idx[0][0]).GetNeighbors():
            if first_atom.GetIdx() not in idx_list:
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
    atom_idxs = list(range(rmol.GetNumAtoms()))
    mapping = (
        atom_idxs[replacement.GetNumAtoms() :] + atom_idxs[: replacement.GetNumAtoms()]
    ) 
    return make_dative_bonds(reset_atom_index(rmol, mapping))


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
        if mol.GetBondBetweenAtoms(start_atom.GetIdx(), end_atom.GetIdx()):
            continue
        if distances[start_atom.GetIdx()][end_atom.GetIdx()] < threshold * (
            pt.GetRcovalent(start_atom.GetAtomicNum())
            + pt.GetRcovalent(end_atom.GetAtomicNum())
        ):
            return False
    return True


def get_crowding_socre(mol: RdMol, threshold: float) -> float:
    """
    Calculate the crowding score of the input molecule.

    Parameters:
        mol (RdMol):
            The input molecule.
        threshold (float):
            The threshold of crowding. If too crowded
    Returns:
        float: The crowding score of the input molecule.
    """
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
            moloplogger.debug(
                f"{DEBUG_TAG}: Distance between {start_atom.GetIdx()} and {end_atom.GetIdx()}"
                f" is {distance} and ideal distance is {ideal_distance}"
            )
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
    *,
    replacement_relative_idx: int = 0,
    replacement_absolute_idx: Union[int, None] = None,
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
    idx = r.AddAtom(Chem.Atom("At"))
    r.AddBond(0, idx, Chem.BondType.SINGLE)
    r.UpdatePropertyCache()
    Chem.SanitizeMol(r)
    cmap = {
        atom_idx: replacement.GetConformer().GetAtomPosition(atom_idx)
        for atom_idx in range(replacement.GetNumAtoms())
    }
    EmbedMolecule(r, randomSeed=randomSeed, coordMap=cmap)
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


def omol_to_rdmol_by_graph(omol: pybel.Molecule) -> Chem.rdchem.Mol:
    """
    Convert a pybel molecule to a rdkit molecule by graph.
    Parameters:
        omol (pybel.Molecule): The pybel molecule to be converted.
    Returns:
        Chem.rdchem.Mol: The rdkit molecule.
    """
    bonds = [
        (bond.GetBeginAtomIdx() - 1, bond.GetEndAtomIdx() - 1, bond.GetBondOrder())
        for bond in ob.OBMolBondIter(omol.OBMol)
    ]
    formal_charges = [atom.GetFormalCharge() for atom in ob.OBMolAtomIter(omol.OBMol)]
    formal_radicals = [
        get_radical_number(atom) for atom in ob.OBMolAtomIter(omol.OBMol)
    ]
    rwmol = Chem.RWMol(Chem.MolFromXYZBlock(omol.write("xyz")))
    for bond in bonds:
        rwmol.AddBond(bond[0], bond[1], bond_list[bond[2]])
    for atom, charge, radical in zip(rwmol.GetAtoms(), formal_charges, formal_radicals):
        atom.SetNoImplicit(True)
        atom.SetFormalCharge(charge)
        atom.SetNumRadicalElectrons(radical)
    rdmol = Chem.MolFromMolBlock(Chem.MolToMolBlock(rwmol), removeHs=False)
    if rdmol is None:
        return None
    return rdmol


def rdmol_to_omol(rdmol: Chem.rdchem.Mol) -> pybel.Molecule:
    """
    Convert a rdkit molecule to a pybel molecule.
    Parameters:
        rdmol (Chem.rdchem.Mol): The rdkit molecule to be converted.
    Returns:
        pybel.Molecule: The pybel molecule.
    """
    return pybel.readstring("sdf", Chem.MolToMolBlock(rdmol, forceV3000=True))


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
        temp_rwmol = Chem.RWMol(rwmol)
        for metal_atom in shuffled_metal_atoms:
            negative_atoms = [
                atom.GetIdx()
                for atom in temp_rwmol.GetAtoms()
                if atom.GetIdx() != metal_atom
                and atom.GetFormalCharge() < 0
                and temp_rwmol.GetConformer()
                .GetAtomPosition(atom.GetIdx())
                .Distance(temp_rwmol.GetConformer().GetAtomPosition(metal_atom))
                <= ratio
                * (
                    pt.GetRcovalent(
                        temp_rwmol.GetAtomWithIdx(metal_atom).GetAtomicNum()
                    )
                    + pt.GetRcovalent(
                        temp_rwmol.GetAtomWithIdx(atom.GetIdx()).GetAtomicNum()
                    )
                )
            ]
            negative_atoms.sort(
                key=lambda x: temp_rwmol.GetConformer()
                .GetAtomPosition(x)
                .Distance(temp_rwmol.GetConformer().GetAtomPosition(metal_atom))
            )
            while (
                temp_rwmol.GetAtomWithIdx(metal_atom).GetFormalCharge() > 0
                and negative_atoms
            ):
                negative_atom = negative_atoms.pop(0)
                temp_rwmol.AddBond(
                    metal_atom,
                    negative_atom,
                    bond_list[
                        -temp_rwmol.GetAtomWithIdx(negative_atom).GetFormalCharge()
                    ],
                )
                temp_rwmol.GetAtomWithIdx(metal_atom).SetFormalCharge(
                    temp_rwmol.GetAtomWithIdx(metal_atom).GetFormalCharge()
                    + temp_rwmol.GetAtomWithIdx(negative_atom).GetFormalCharge()
                )
                temp_rwmol.GetAtomWithIdx(negative_atom).SetFormalCharge(0)
        for metal_atom in metal_atoms:
            for dative_atom in [
                idxs[0]
                for idxs in temp_rwmol.GetSubstructMatches(
                    Chem.MolFromSmarts("[Ov2+0,Sv2+0,Nv3+0,Pv3+0]")
                )
            ]:
                if not temp_rwmol.GetBondBetweenAtoms(metal_atom, dative_atom):
                    if temp_rwmol.GetConformer().GetAtomPosition(metal_atom).Distance(
                        temp_rwmol.GetConformer().GetAtomPosition(dative_atom)
                    ) < exp_ratio[
                        temp_rwmol.GetAtomWithIdx(dative_atom).GetSymbol()
                    ] * (
                        pt.GetRcovalent(
                            temp_rwmol.GetAtomWithIdx(metal_atom).GetAtomicNum()
                        )
                        + pt.GetRcovalent(
                            temp_rwmol.GetAtomWithIdx(dative_atom).GetAtomicNum()
                        )
                    ):
                        temp_rwmol.AddBond(
                            dative_atom, metal_atom, Chem.BondType.DATIVE
                        )
        remained_negative_atoms = [
            atom.GetIdx()
            for atom in temp_rwmol.GetAtoms()
            if atom.GetFormalCharge() < 0 and not is_metal(atom.GetAtomicNum())
        ]
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
            if temp_rwmol.GetConformer().GetAtomPosition(metal_atoms[0]).Distance(
                temp_rwmol.GetConformer().GetAtomPosition(remained_negative_atom)
            ) <= ratio * (
                pt.GetRcovalent(
                    temp_rwmol.GetAtomWithIdx(metal_atoms[0]).GetAtomicNum()
                )
                + pt.GetRcovalent(
                    temp_rwmol.GetAtomWithIdx(remained_negative_atom).GetAtomicNum()
                )
            ):
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
