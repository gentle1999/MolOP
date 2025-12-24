from openbabel import openbabel as ob
from openbabel import pybel
from rdkit import Chem

from molop.config import moloplogger
from molop.structure.utils import bond_list, bond_type_mapping


DEBUG_TAG = "[FORMAT CONVERTER]"


def rdmol_to_omol(rdmol: Chem.rdchem.Mol) -> pybel.Molecule:
    """
    Convert a rdkit molecule to a pybel molecule.
    Parameters:
        rdmol (Chem.rdchem.Mol): The rdkit molecule to be converted.
    Returns:
        pybel.Molecule: The pybel molecule.
    """
    return pybel.readstring("sdf", Chem.MolToMolBlock(rdmol, forceV3000=True))


def omol_to_rdmol_by_graph(omol: pybel.Molecule) -> Chem.rdchem.Mol | None:
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
    formal_radicals = [atom.GetSpinMultiplicity() for atom in ob.OBMolAtomIter(omol.OBMol)]
    rwmol = Chem.RWMol(Chem.MolFromXYZBlock(omol.write("xyz")))
    for bond in bonds:
        rwmol.AddBond(bond[0], bond[1], bond_list[bond[2]])
    for atom, charge, radical in zip(
        rwmol.GetAtoms(), formal_charges, formal_radicals, strict=True
    ):
        atom.SetNoImplicit(True)
        atom.SetFormalCharge(charge)
        atom.SetNumRadicalElectrons(radical)
    return Chem.MolFromMolBlock(Chem.MolToMolBlock(rwmol), removeHs=False)


def validate_rdmol(rdmol: Chem.rdchem.Mol, total_charge=0, total_radical=0) -> bool:
    """
    Check whether the rdmol is valid or not.
    Parameters:
        rdmol (Chem.rdchem.Mol): The rdmol to be checked.
        total_charge (int): The total charge of the target molecule.
        total_radical (int): The total radical number of the target molecule.
    Returns:
        bool: Whether the rdmol is valid or not.
    """
    if Chem.SanitizeMol(rdmol, catchErrors=True) != 0:
        return False
    charge = sum(atom.GetFormalCharge() for atom in rdmol.GetAtoms())
    radical = sum(atom.GetNumRadicalElectrons() for atom in rdmol.GetAtoms())
    # singlet state
    radical_in_singlet = sum(atom.GetNumRadicalElectrons() % 2 for atom in rdmol.GetAtoms())
    if total_radical == radical_in_singlet:
        radical = radical_in_singlet
    moloplogger.debug(
        f"{DEBUG_TAG} | check rdmol, target charge: {total_charge} charge: {charge}, "
        f"target radical: {total_radical} radical: {radical}"
    )
    return charge == total_charge and radical == total_radical


def omol_to_rdmol(omol: pybel.Molecule, total_charge=0, total_radical=0) -> Chem.rdchem.Mol | None:
    """
    Convert a pybel molecule to a rdkit molecule.

    Parameters:
        omol (pybel.Molecule): The pybel molecule to be converted.
        total_charge (int): The total charge of the target molecule.
        total_radical (int): The total radical number of the target molecule.
    Returns:
        Chem.rdchem.Mol: The rdkit molecule.
    """
    moloplogger.debug(f"{DEBUG_TAG} | try tranform by mol")
    if rdmol := Chem.MolFromMolBlock(omol.write("mol"), removeHs=False):  # noqa: SIM102
        if validate_rdmol(rdmol, total_charge, total_radical):
            moloplogger.debug(f"{DEBUG_TAG} | success by graph, mol: {Chem.MolToSmiles(rdmol)}")
            return rdmol
    moloplogger.debug(f"{DEBUG_TAG} | try tranform by graph")
    if rdmol := omol_to_rdmol_by_graph(omol):  # noqa: SIM102
        if validate_rdmol(rdmol, total_charge, total_radical):
            moloplogger.debug(f"{DEBUG_TAG} | success by graph, smiles: {Chem.MolToSmiles(rdmol)}")
            return rdmol
    return None


def rdmol_to_gjf_geom(rdmol: Chem.rdchem.Mol, total_radical=0) -> str:
    """
    Convert a rdkit molecule to a gjf geometry block.

    Parameters:
        rdmol (Chem.rdchem.Mol): The rdkit molecule to be converted.
        total_radical (int): The total radical number of the target molecule.
    Returns:
        str: The gjf geometry block.
    """
    charge = sum(atom.GetFormalCharge() for atom in rdmol.GetAtoms())
    radical = sum(atom.GetNumRadicalElectrons() for atom in rdmol.GetAtoms())
    # singlet state
    radical_in_singlet = sum(atom.GetNumRadicalElectrons() % 2 for atom in rdmol.GetAtoms())
    if total_radical == radical_in_singlet:
        radical = radical_in_singlet
    spin_multiplicity = radical + 1
    return f"{charge} {spin_multiplicity}\n" + "\n".join(
        [
            f"{atom:10s}{x:18.10f}{y:18.10f}{z:18.10f}"
            for atom, (x, y, z) in zip(
                [atom.GetSymbol() for atom in rdmol.GetAtoms()],
                rdmol.GetConformer().GetPositions(),
                strict=True,
            )
        ]
    )


def rdmol_to_gjf_connectivity(rdmol: Chem.rdchem.Mol) -> str:
    """
    Convert a rdkit molecule to a gjf connectivity block.

    Parameters:
        rdmol (Chem.rdchem.Mol): The rdkit molecule to be converted.
    Returns:
        str: The gjf connectivity block.
    """
    lines = []
    for atom_idx in range(rdmol.GetNumAtoms()):
        atom = rdmol.GetAtomWithIdx(atom_idx)
        gjf_idx = atom_idx + 1
        line_parts = [str(gjf_idx)]
        for bond in atom.GetBonds():
            begin_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
            neighbor_idx = end_idx if begin_idx == atom_idx else begin_idx
            gjf_neighbor_idx = neighbor_idx + 1
            btype = bond.GetBondType()
            order = bond_type_mapping.get(btype, "1.0")
            line_parts.append(f"{gjf_neighbor_idx} {order}")
        lines.append(" ".join(line_parts))
    return "\n".join(lines)
