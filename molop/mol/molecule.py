"""
Author: TMJ
Date: 2023-03-22 19:14:54
LastEditors: TMJ
LastEditTime: 2023-10-30 15:17:26
Description: 请填写简介
"""
import hashlib
from copy import deepcopy
from typing import Callable, List, Union

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors, Draw

from molop.utils.errors import ArgumentError, MolError, ReactError
from molop.utils.types import RdMol

RDLogger.DisableLog("rdApp.*")


class Molecule:
    """
    Molecule base class.
    Parameters:
        smiles str:
            Optional, default is None\n
            The SMILES string is a string of the molecule's SMILES.\n
            If you pass a SMILES string, it will be converted to a RDKit molecule.
        rd_mol RdMol:
            RdMol, optional, default is None\n
            The RDKit molecule is a molecule object that can be converted to a SMILES string.\n
            If you pass a RDKit molecule, it will be converted to a SMILES string.\n
            If the RDKit molecule is not provided, the SMILES string is generated from the RDKit molecule.
        sanitize bool:
            Optional, default is True\n
            If True, sanitize the RDKit molecule.

    How to use
    ----------
        >>> import mechaevo as me
        >>> mol = me.Molecule(smiles='CC')

    or

        >>> from rdkit import Chem
        >>> import mechaevo as me
        >>> rdmol = Chem.MolFromSmiles('CC')
        >>> mol = me.Molecule(mol=rdmol)
    """

    def __init__(
        self,
        smiles: str = None,
        rd_mol: RdMol = None,
        sanitize: bool = True,
        use_chirality: bool = False,
    ):
        self._is_none = False
        if rd_mol:
            self._rd_mol = rd_mol
            self._smiles = Chem.MolToSmiles(rd_mol)
        elif smiles:
            self._rd_mol = Chem.MolFromSmiles(smiles)
            self._smiles = Chem.MolToSmiles(self._rd_mol)
        else:
            self._rd_mol = None
            self._smiles = None
            self._is_none = True
        self._is_sanitized = False
        if sanitize:
            self.sanitize()
        self._use_chirality = use_chirality
        # TODO update comments
        self._check_chirality()
        if not use_chirality:
            self._remove_chirality()
        self.addHs()
        try:
            self._coords = self._rd_mol.GetConformer().GetPositions()
        except:
            self._coords = None

    @property
    def local_charge(self) -> dict:
        """
        Get the local charge of the molecule.
        Returns:
            Dict of the local charges of the atoms by index.
        """
        if self._is_none:
            return None
        local_charge_dict = {}
        for index in range(len(self)):
            charge = self._rd_mol.GetAtomWithIdx(index).GetFormalCharge()
            if charge != 0:
                local_charge_dict[index] = charge
        return local_charge_dict

    @property
    def global_charge(self) -> int:
        """
        Get the global charge of the molecule.
        Returns:
            Totally formal charge of molecule.
        """
        return None if self._is_none else Chem.GetFormalCharge(self._rd_mol)

    @property
    def local_radical_electrons(self) -> dict:
        if self._is_none:
            return None
        local_radical_electrons_dict = {}
        for index in range(len(self)):
            radical = self._rd_mol.GetAtomWithIdx(index).GetNumRadicalElectrons()
            if radical != 0:
                local_radical_electrons_dict[index] = radical
        return local_radical_electrons_dict

    @property
    def radical_electrons(self) -> int:
        return None if self._is_none else Descriptors.NumRadicalElectrons(self._rd_mol)

    @property
    def multi(self) -> int:
        return None if self._is_none else 2 * self.radical_electrons + 1

    @property
    def formula(self) -> str:
        elements = self.get_elements()
        return "".join(
            [f"{key}{value}" for key, value in elements.items() if value > 0]
        )

    @property
    def is_use_chirality(self):
        return self._use_chirality

    def _check_chirality(self):
        if self._is_none or not self._use_chirality:
            return
        # print(self._smiles)
        if len(Chem.FindMolChiralCenters(self._rd_mol)) > 0:
            self._has_chirality = True

    @property
    def has_chirality(self):
        try:
            return self._has_chirality
        except Exception:
            self._check_chirality()
            return self._has_chirality

    def _remove_chirality(self):
        if self._is_none:
            return
        for idx, chiral_type in Chem.FindMolChiralCenters(
            self._rd_mol, includeUnassigned=True
        ):
            self._rd_mol.GetAtomWithIdx(idx).SetChiralTag(
                Chem.rdchem.ChiralType.CHI_UNSPECIFIED
            )
        self._use_chirality = False

    def get_elements(self):
        pt = Chem.GetPeriodicTable()
        elements = {pt.GetElementSymbol(idx): 0 for idx in range(1, 55)}
        mol = self.copy()
        mol.addHs()
        for atom in mol._rd_mol.GetAtoms():
            elements[atom.GetSymbol()] += 1
        return elements

    def path_name(self, postfix=None) -> str:
        if postfix is not None:
            return f"{self.formula}@{postfix}"
        md5_machine = hashlib.md5()
        md5_machine.update(self._smiles.encode("utf-8"))
        post_fix = md5_machine.hexdigest()[:4]
        return f"{self.formula}@{post_fix}"

    def sanitize(self):
        """Sanitize the molecule."""
        if self._is_none or self._is_sanitized:
            return
        try:
            Chem.SanitizeMol(self._rd_mol)
        except:
            self._rd_mol = Chem.MolFromSmiles(self._smiles, sanitize=False)
            raise MolError(f'"{self._smiles}" can not be sanitized')

        self._smiles = Chem.MolToSmiles(self._rd_mol)
        self._is_sanitized = True

    def addHs(self):
        """Add hydrogens to the molecule."""
        if self._is_none:
            return
        self._rd_mol = Chem.AddHs(self._rd_mol)
        self._smiles = Chem.MolToSmiles(self._rd_mol)

    def removeHs(self):
        """Remove hydrogens to the molecule."""
        if self._is_none:
            return
        self._rd_mol = Chem.RemoveHs(self._rd_mol)
        self._smiles = Chem.MolToSmiles(self._rd_mol)

    @property
    def coords(self):
        """Get the coordinates of the molecule."""
        if self._coords:
            self._coords = self._rd_mol.GetConformer().GetPositions()
            return self._coords
        else:
            self.embedMolecule()
            return self._coords

    def embedMolecule(self, random_state=3407):
        """
        Embed the molecule for the given random seed.
        Args:
            random_state int:
                random state for embedding.
        """
        if self._is_none:
            return
        self.addHs()
        if self.global_charge > 0:
            tmp_rd_mol = deepcopy(self._rd_mol)
            if (
                AllChem.EmbedMolecule(
                    tmp_rd_mol,
                    randomSeed=random_state,
                    enforceChirality=self._use_chirality,
                )
                == -1
            ):
                self._rd_mol = self.embedMolecule_cation(random_state=random_state)
            else:
                self._rd_mol = tmp_rd_mol
        else:
            AllChem.EmbedMolecule(
                self._rd_mol,
                randomSeed=random_state,
                enforceChirality=self._use_chirality,
            )
        self._coords = self._rd_mol.GetConformer().GetPositions()
        self._smiles = Chem.MolToSmiles(self._rd_mol)

    def embedMolecule_cation(self, random_state=3407) -> RdMol:
        """
        Embed the cation for the given random seed. This is to avoid the fragmentation of some cations with special structures under the ETKDG algorithm.
        Args:
            random_state int:
                random state for embedding.
        Returns:
            A cation with a conformation containing a hydrogen atom
        """
        mol_add_H = deepcopy(self._rd_mol)
        index = []
        for i in mol_add_H.GetSubstructMatch(Chem.MolFromSmarts("[C+]")):
            mol_add_H.GetAtomWithIdx(i).SetFormalCharge(0)
            mol_add_H.GetAtomWithIdx(i).SetNumRadicalElectrons(0)
            index.append(i)
        mol_add_H = AllChem.AddHs(mol_add_H)
        AllChem.EmbedMolecule(
            mol_add_H, randomSeed=random_state, enforceChirality=self._use_chirality
        )
        for i in index:
            mol_add_H.GetAtomWithIdx(i).SetFormalCharge(1)
            mol_add_H.GetAtomWithIdx(i).SetNumRadicalElectrons(0)
        return mol_add_H

    def save_xyz(self, path, random_state=3407):
        """
        Save the RDKIT molecules to a .xyz file. Will perform conformational optimization of molecules using embed and preserve the stereo structure.
        Args:
            path str:
                path to the file where the molecules will be written.
            random_state int:
                random state for embedding.
        """
        if self._is_none:
            return None
        if self._coords:
            pass
        mol_add_H = self.copy()
        mol_add_H.addHs()
        mol_add_H.embedMolecule(random_state=random_state)
        Chem.MolToXYZFile(mol_add_H.to_mol(), path)

    def __eq__(self, other: "Molecule") -> bool:
        """
        Check if two molecules are equal. If the two molecules are equal, the Molecule A is the substructure of the Molecule B as well as the Molecule B is the substructure of the Molecule A
        Args:
            other 'Molecule': The other object to compare with.
        Returns:
            True if both molecules are the same.
        """
        if not isinstance(other, (Molecule)):
            return False
        if self._is_none and other._is_none:
            return True
        if self._is_none or other._is_none:
            return False
        return self._rd_mol.HasSubstructMatch(
            other._rd_mol, useChirality=self._use_chirality
        ) and other._rd_mol.HasSubstructMatch(
            self._rd_mol, useChirality=other._use_chirality
        )

    def is_equal(
        self,
        other: Union["Molecule", RdMol, str],
        use_chirality=False,
    ) -> bool:
        if isinstance(other, (Molecule)):
            return self == other
        if isinstance(other, RdMol):
            return self == Molecule(rd_mol=other, use_chirality=use_chirality)
        return (
            self == Molecule(smiles=other, use_chirality=use_chirality)
            if isinstance(other, str)
            else False
        )

    def __lt__(self, other: "Molecule") -> bool:
        """
        Check if the first given object is less than the second.
        Following the rules: 1. Number of heavy atoms; 2. MolWt
        Args:
            other Union['Molecule', 'Mol']: The other object to compare with.
        Returns:
            True if the first object is less than the second.
        """
        if not isinstance(other, (Molecule)):
            raise MolError(
                'Comparisons can only be made with class "Molecule"'
            )
        if self._is_none:
            return True
        elif other._is_none:
            return False
        if self._rd_mol.GetNumAtoms() < other._rd_mol.GetNumAtoms():
            return True
        elif Descriptors.MolWt(self._rd_mol) < Descriptors.MolWt(other._rd_mol):
            return True

    def __len__(self) -> int:
        """
        Returns:
            The number of heavy atoms in the molecule.
        """
        return 0 if self._is_none else self._rd_mol.GetNumAtoms()

    @property
    def display_rd_mol(self) -> RdMol:
        if self._is_none:
            return None
        mol_no_H = Chem.RemoveHs(self._rd_mol)
        for atom in mol_no_H.GetAtoms():
            atom.SetAtomMapNum(0)
        Chem.SanitizeMol(mol_no_H)
        mol_no_H.RemoveAllConformers()
        return mol_no_H

    @property
    def display_smiles(self) -> str:
        return "None" if self._is_none else Chem.MolToSmiles(self.display_rd_mol)

    def __str__(self) -> str:
        """
        Returns:
            The SMILES string of molecule.
        """
        return self.display_smiles

    def __repr__(self) -> str:
        chirality_flag = " chirality" if self._use_chirality else ""
        return f"Melecule({self}{chirality_flag})"

    def __hash__(self) -> int:
        """
        Returns:
            Hash of SMILES string.
        """
        return 0 if self._is_none else hash(self._smiles)

    def get_total_atoms_number(self) -> int:
        """
        Returns:
            The total number of atoms in the molecule (Including hydrogen atoms).
        """
        m = self.copy()
        m.addHs()
        return len(m)

    def draw(self, **kwargs):
        """
        Draw the molecule.
        """
        return None if self._is_none else Draw.MolToImage(self.display_rd_mol, **kwargs)

    def to_smiles(self) -> str:
        """
        Returns:
            The SMILES string.
        """
        return self._smiles

    def to_mol(self) -> RdMol:
        """
        Returns:
            The RDKIT Molecule.
        """
        return self._rd_mol

    def copy(self) -> "Molecule":
        """
        Returns:
            Copy of self.
        """
        return deepcopy(self)

    def refresh(self):
        """
        Use SMILES to regenerate RDKIT mol for the intermediary.
        """
        if self._is_none:
            return
        self._rd_mol = Chem.MolFromSmiles(Chem.MolToSmiles(self._rd_mol))
        self._rd_mol = Chem.AddHs(self._rd_mol)
