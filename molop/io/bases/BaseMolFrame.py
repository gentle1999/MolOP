"""
Author: TMJ
Date: 2024-06-17 20:42:47
LastEditors: TMJ
LastEditTime: 2024-06-17 20:47:40
Description: 请填写简介
"""

import os
from typing import Dict, List, Literal, Sequence, Tuple, Union

import numpy as np
from openbabel import openbabel as ob
from openbabel import pybel
from pint.facets.plain import PlainQuantity
from pydantic import Field, PrivateAttr, computed_field
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds

from molop.io.bases.DataClasses import BaseDataClassWithUnit
from molop.logger.logger import logger
from molop.structure.geometry import get_geometry_info
from molop.structure.structure import (
    bond_list,
    get_bond_pairs,
    get_formal_charges,
    get_formal_num_radicals,
)
from molop.structure.structure_recovery import omol_to_rdmol, xyz_block_to_omol
from molop.unit import atom_ureg
from molop.utils.types import RdMol

pt = Chem.GetPeriodicTable()


class BaseMolFrame(BaseDataClassWithUnit):
    atoms: List[int] = Field(
        default=[],
        description="atom numbers",
        title="Atom numbers",
    )
    coords: Union[PlainQuantity, None] = Field(
        default=np.array([[]]) * atom_ureg.angstrom,
        description="Atom coordinates, unit is `angstrom`",
        title="Atom coordinates",
    )
    standard_coords: Union[PlainQuantity, None] = Field(
        default=np.array([[]]) * atom_ureg.angstrom,
        description="Atom coordinates in standard orientation, default unit is `angstrom`",
    )
    charge: int = Field(default=0, description="Molecule total charge")
    multiplicity: int = Field(default=1, description="Molecule total multiplicity")

    @computed_field
    @property
    def smiles(self) -> str:
        return self.to_standard_SMILES()

    bonds: List[Tuple[int, int, int]] = Field(
        default=[],
        description="Bond information, each bond is represented by a tuple of three integers, "
        "where the first two integers represent the indices of the two atoms involved in the bond, "
        "and the third integer represents the bond order (follow the RDKit convention)",
    )
    formal_charges: List[int] = Field(
        default=[],
        description="Formal charges of each atom",
    )
    formal_num_radicals: List[int] = Field(
        default=[],
        description="Number of radical electrons of each atom",
    )
    _rdmol: RdMol = PrivateAttr(default=None)

    @property
    def atom_symbols(self) -> List[str]:
        """
        Get the atom symbols.
        """
        return [Chem.Atom(atom).GetSymbol() for atom in self.atoms]

    @property
    def total_electrons(self) -> int:
        """
        Get the total electrons.

        Returns:
            The total electrons.
        """
        return sum(Chem.Atom(atom).GetAtomicNum() for atom in self.atoms) + self.charge

    @property
    def elements(self) -> List[str]:
        """
        Get the elements set.

        Returns:
            A list of element symbols dropped duplicates.
        """
        return list(set(self.atoms))

    @property
    def formula(self) -> str:
        """
        Get the formula.
        """
        return "".join(
            [
                f"{pt.GetElementSymbol(i)}{self.atoms.count(pt.GetElementSymbol(i))}"
                for i in range(1, 119)
                if self.atoms.count(pt.GetElementSymbol(i)) != 0
            ]
        )

    def __get_topology(self):
        if self.rdmol is None:
            return
        self.bonds = get_bond_pairs(self.rdmol)
        self.formal_charges = get_formal_charges(self.rdmol)
        self.formal_num_radicals = get_formal_num_radicals(self.rdmol)

    @property
    def omol(self) -> pybel.Molecule:
        """
        Get the openbabel molecule object.

        Returns:
            The openbabel molecule object.
        """
        if len(self.bonds):
            omol = pybel.readstring("xyz", self.to_XYZ_block())
            kekulized_rdmol = Chem.RWMol(self.rdmol)
            Chem.KekulizeIfPossible(kekulized_rdmol, clearAromaticFlags=True)
            for bond in ob.OBMolBondIter(omol.OBMol):
                omol.OBMol.DeleteBond(bond)
            for bond in kekulized_rdmol.GetBonds():
                bond_order = bond_list.index(bond.GetBondType())
                if bond_order > 3:
                    bond_order = 1
                omol.OBMol.AddBond(
                    bond.GetBeginAtomIdx() + 1, bond.GetEndAtomIdx() + 1, bond_order
                )
                omol.OBMol.GetBond(
                    omol.OBMol.GetAtomById(bond.GetBeginAtomIdx()),
                    omol.OBMol.GetAtomById(bond.GetEndAtomIdx()),
                ).SetBondOrder(bond_order)
            for atom, charge, radical in zip(
                omol.atoms, self.formal_charges, self.formal_num_radicals
            ):
                atom.OBAtom.SetFormalCharge(charge)
            return omol
        try:
            return xyz_block_to_omol(
                self.to_XYZ_block(), self.charge, int(self.multiplicity) - 1
            )
        except Exception as e:
            raise RuntimeError(f"{self.file_path}: {e}")

    @property
    def rdmol(self) -> Chem.rdchem.Mol:
        """
        Get the rdkit molecule object.

        Returns:
            The rdkit molecule object.
        """
        if self._rdmol is None:
            if len(self.bonds) == 0:
                try:
                    # try rdkit determinebonds first
                    # Issues known:
                    # - Can not recognize radicals
                    # - Can not recognize Metals
                    if len(self.atoms) <= 1:
                        raise ValueError("Single atom is not allowed in this function.")
                    raw_mol = Chem.MolFromXYZBlock(self.to_XYZ_block())
                    conn_mol = Chem.Mol(raw_mol)
                    rdDetermineBonds.DetermineBonds(conn_mol, charge=self.charge)
                    self._rdmol = conn_mol
                    self.__get_topology()
                except:
                    logger.debug(
                        f"{self.file_path}: rdkit determinebonds failed. Use MolOP structure recovery instead."
                    )
                    # If failed, use MolOP implementation
                    try:
                        omol = self.omol
                    except Exception as e:
                        logger.error(
                            f"{self.file_path}: MolOP structure recovery failed. {e}"
                        )
                        return None

                    try:
                        self._rdmol = omol_to_rdmol(
                            omol, self.charge, self.multiplicity - 1
                        )
                    except Exception as e:
                        logger.error(
                            f"{self.file_path}: RDKit structure recovery failed. {e}"
                        )
                    finally:
                        self.__get_topology()
                        return self._rdmol
            else:
                assert (
                    self.formal_charges and self.formal_num_radicals
                ), "If bonds given, formal charges and spins must be provided."
                rwmol = Chem.RWMol(Chem.MolFromXYZBlock(self.to_XYZ_block()))
                for bond in self.bonds:
                    rwmol.AddBond(bond[0], bond[1], bond_list[bond[2]])
                for atom, charge, spin in zip(
                    rwmol.GetAtoms(), self.formal_charges, self.formal_num_radicals
                ):
                    atom.SetFormalCharge(charge)
                    atom.SetNumRadicalElectrons(spin)
                self._rdmol = rwmol.GetMol()
        return self._rdmol

    @property
    def rdmol_no_conformer(self):
        """
        Get the rdkit molecule object without conformer.

        Returns:
            The rdkit molecule object without conformer.
        """
        rdmol = Chem.RWMol(self.rdmol)
        rdmol.RemoveAllConformers()
        return rdmol

    def to_XYZ_block(self) -> str:
        """
        Get the XYZ block.

        Returns:
            The XYZ block.
        """
        if self.coords is None:
            raise ValueError("No coordinates found!")
        return (
            f"{len(self.atoms)}\n"
            + f"charge {self.charge} multiplicity {self.multiplicity}\n"
            + "\n".join(
                [
                    f"{atom:10s}{x:10.5f}{y:10.5f}{z:10.5f}"
                    for atom, (x, y, z) in zip(self.atom_symbols, self.coords.m)
                ]
            )
        )

    def to_SDF_block(self) -> str:
        """
        Get the SDF block.

        Returns:
            The SDF block.
        """
        return Chem.MolToMolBlock(self.rdmol)

    def to_CML_block(self) -> str:
        """
        Get the CML block.

        Returns:
            The CML block.
        """
        return Chem.MolToMrvBlock(self.rdmol)

    def to_SMILES(self) -> str:
        """
        Get the SMILES with explicit hydrogens.

        Returns:
            The SMILES.
        """
        if self.rdmol is None:
            return ""
        smiles = Chem.MolToSmiles(self.rdmol)
        return smiles

    def to_standard_SMILES(self) -> str:
        """
        Get the SMILES with standardization.

        Returns:
            The SMILES.
        """
        # return rdMolStandardize.StandardizeSmiles(self.to_SMILES())
        return Chem.CanonSmiles(self.to_SMILES(), useChiral=True)

    def to_InChI(self) -> str:
        """
        Get the InChI.

        Returns:
            The InChI.
        """
        return Chem.MolToInchi(self.rdmol)

    def calc_rdkit_descs(self, desc_names: List[str] = None) -> Dict[str, float]:
        """
        Calculate the RDKit descriptors.

        Parameters:
            desc_names List[str]:
                The names of the descriptors. Must be a subset of the RDKit descriptors.

        Returns:
            The dictionary of the descriptors.
        """
        from molop.descriptor.descriptor import calc_rdkit_descs

        return calc_rdkit_descs(self.rdmol, desc_names=desc_names)

    def calc_dscribe_descs(
        self,
        desc_names: List[
            Literal[
                "SOAP",
                "ACSF",
                "MBTR",
                "LMBTR",
            ]
        ] = None,
    ) -> Dict[str, np.ndarray]:
        """
        Calculate the dscribe descriptors.

        Require additional dependencies:
            - [dscribe](https://singroup.github.io/dscribe/latest/)
            - [ase](https://wiki.fysik.dtu.dk/ase/index.html)

        Use the following command to install extra dependencies:
        ```bash
        pip install --extra-index-url http://10.72.201.58:13000/api/packages/tmj/pypi/simple/ --trusted-host 10.72.201.58 molop[full] --upgrade
        ```

        Parameters:
            desc_names List[str]:
                The names of the descriptors. Must be a subset of "SOAP", "ACSF", "MBTR", "LMBTR"

        Returns:
            The dictionary of the descriptors.
        """
        from molop.descriptor.descriptor import calc_dscribe_descs

        self.to_SDF_file(".temp.sdf")
        ans = calc_dscribe_descs(
            species=self.elements,
            atoms_path=".temp.sdf",
            desc_names=desc_names,
        )
        os.remove(".temp.sdf")
        return ans

    def calc_mordred_descs(self, **kwargs) -> Dict[str, float]:
        """
        Calculate the Mordred descriptors.

        Require additional dependencies:
            - [Mordred](https://github.com/mordred-descriptor/mordred)

        Use the following command to install extra dependencies:
        ```bash
        pip install --extra-index-url http://10.72.201.58:13000/api/packages/tmj/pypi/simple/ --trusted-host 10.72.201.58 molop[full] --upgrade
        ```

        Returns:
            The dictionary of the descriptors.
        """

        from molop.descriptor.descriptor import calc_mordred_descs

        return calc_mordred_descs(self.rdmol, **kwargs)

    def __hash__(self) -> int:
        return hash(str(self))

    def __len__(self) -> int:
        return len(self.atoms)

    def geometry_analysis(self, atom_idxs: Sequence[int], one_start=False) -> float:
        """
        Get the geometry infos among the atoms

        Parameters:
            atom_idxs Sequence[int]:
                A Sequence of index of the atoms, starts from 0
            one_start bool:
                If true, consider atom index starts from 1, so let index value subtracts 1 for all the atoms

        Returns:
            - If the length of atom_idxs is 2, the bond length with unit Angstrom between the two atoms will be returned.
            - If the length of atom_idxs is 3, the angle with unit degree between  the three atoms will be returned.
            - If the length of atom_idxs is 4, the dihedral angle with unit degree between the four atoms will be returned.
        """
        if one_start:
            atom_idxs = [atom_idx - 1 for atom_idx in atom_idxs]
        return get_geometry_info(self.rdmol, atom_idxs)
