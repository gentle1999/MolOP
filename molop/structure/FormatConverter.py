from typing import Optional

from openbabel import openbabel as ob
from openbabel import pybel
from rdkit import Chem

from molop.config import moloplogger
from molop.structure.utils import bond_list

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


def omol_to_rdmol_by_graph(omol: pybel.Molecule) -> Optional[Chem.rdchem.Mol]:
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
        atom.GetSpinMultiplicity() for atom in ob.OBMolAtomIter(omol.OBMol)
    ]
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
    radical_in_singlet = sum(
        atom.GetNumRadicalElectrons() % 2 for atom in rdmol.GetAtoms()
    )
    if total_radical == radical_in_singlet:
        radical = radical_in_singlet
    moloplogger.debug(
        f"{DEBUG_TAG} | check rdmol, target charge: {total_charge} charge: {charge}, "
        f"target radical: {total_radical} radical: {radical}"
    )
    return charge == total_charge and radical == total_radical


def omol_to_rdmol(
    omol: pybel.Molecule, total_charge=0, total_radical=0
) -> Optional[Chem.rdchem.Mol]:
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
    if rdmol := Chem.MolFromMolBlock(omol.write("mol"), removeHs=False):
        if validate_rdmol(rdmol, total_charge, total_radical):
            moloplogger.debug(
                f"{DEBUG_TAG} | success by graph, mol: {Chem.MolToSmiles(rdmol)}"
            )
            return rdmol
    moloplogger.debug(f"{DEBUG_TAG} | try tranform by graph")
    if rdmol := omol_to_rdmol_by_graph(omol):
        if validate_rdmol(rdmol, total_charge, total_radical):
            moloplogger.debug(
                f"{DEBUG_TAG} | success by graph, smiles: {Chem.MolToSmiles(rdmol)}"
            )
            return rdmol
    return None


def validate_omol(
    omol: pybel.Molecule, total_charge: int, total_radical_electrons: int
) -> bool:
    """
    Check if the final structure is valid.
    Parameters:
        omol (pybel.Molecule): The pybel molecule to be checked.
        total_charge (int): The total charge of the molecule.
        total_radical_electrons (int): The total radical of the molecule.
    Returns:
        bool: True if the final structure is valid, False otherwise.
    """
    if (
        charge_sum := sum(atom.OBAtom.GetFormalCharge() for atom in omol.atoms)
    ) != total_charge:
        moloplogger.debug(
            f"{DEBUG_TAG} | Charge check failed, total charge: "
            f"{total_charge}, actual charge: {charge_sum}"
        )
        return False

    radical_sum = sum(atom.OBAtom.GetSpinMultiplicity() for atom in omol.atoms)
    radical_sum_singlet = sum(
        atom.OBAtom.GetSpinMultiplicity() % 2 for atom in omol.atoms
    )
    if radical_sum_singlet == total_radical_electrons:
        radical_sum = radical_sum_singlet
    if radical_sum != total_radical_electrons:
        moloplogger.debug(
            f"{DEBUG_TAG} | Radical check failed, total radical: "
            f"{total_radical_electrons}, actual radical: {radical_sum}"
        )
        return False
    return True
