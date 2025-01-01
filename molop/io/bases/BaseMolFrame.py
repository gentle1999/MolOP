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

from molop.descriptor.spms import SPMSCalculator
from molop.io.bases.DataClasses import BaseDataClassWithUnit
from molop.logger.logger import moloplogger
from molop.structure.geometry import calculate_rmsd, get_geometry_info
from molop.structure.structure import (
    bond_list,
    canonical_smiles,
    get_bond_pairs,
    get_formal_charges,
    get_formal_num_radicals,
)
from molop.structure.structure_recovery import rdmol_to_omol, xyz2rdmol
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
        return self.to_SMILES()

    @computed_field()
    @property
    def canonical_smiles(self) -> str:
        return self.to_canonical_SMILES()

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

        Returns:
            List[str]: A list of atom symbols.
        """
        return [Chem.Atom(atom).GetSymbol() for atom in self.atoms]

    @property
    def total_electrons(self) -> int:
        """
        Get the total electrons.

        Returns:
            int: The total electrons.
        """
        return sum(Chem.Atom(atom).GetAtomicNum() for atom in self.atoms) + self.charge

    @property
    def elements(self) -> List[str]:
        """
        Get the elements set.

        Returns:
            List[str]: A list of element symbols dropped duplicates.
        """
        return list(set(self.atoms))

    @property
    def formula(self) -> str:
        """
        Get the formula.

        Returns:
            str: The formula.
        """
        return "".join(
            [
                f"{pt.GetElementSymbol(i)}{self.atoms.count(pt.GetElementSymbol(i))}"
                for i in range(1, 119)
                if self.atoms.count(pt.GetElementSymbol(i)) != 0
            ]
        )

    def __get_topology(self):
        if self._rdmol is None:
            return
        self.bonds = get_bond_pairs(self.rdmol)
        self.formal_charges = get_formal_charges(self.rdmol)
        self.formal_num_radicals = get_formal_num_radicals(self.rdmol)

    @property
    def omol(self) -> pybel.Molecule:
        """
        Get the openbabel molecule object.

        Returns:
            pybel.Molecule: The openbabel molecule object.
        """
        return rdmol_to_omol(self.rdmol)

    @property
    def rdmol(self) -> Union[Chem.rdchem.Mol, None]:
        """
        Get the rdkit molecule object.

        If reconstruction failed, return None.

        Returns:
            Union[Chem.rdchem.Mol,None]: 
                The rdkit molecule object. If reconstruction failed, return None.
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
                    moloplogger.debug(
                        f"{self.file_path}: rdkit determinebonds failed. Use MolOP structure recovery instead."
                    )
                    # If failed, use MolOP implementation
                    try:
                        self._rdmol = xyz2rdmol(
                            self.to_XYZ_block(), self.charge, self.multiplicity - 1
                        )
                    except Exception as e:
                        moloplogger.error(
                            f"{self.file_path}: MolOP structure recovery failed. {e}"
                        )
                    finally:
                        self.__get_topology()
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
        if self._rdmol is not None:
            try:
                Chem.SanitizeMol(self._rdmol)
            except:
                moloplogger.error(
                    f"{self.file_path}: MolOP structure recovery failed. RDKit sanitization failed."
                )
        return self._rdmol

    @property
    def rdmol_no_conformer(self) -> Chem.rdchem.Mol:
        """
        Get the rdkit molecule object without conformer.

        Returns:
            Chem.rdchem.Mol: The rdkit molecule object without conformer.
        """
        rdmol = Chem.RWMol(self.rdmol)
        rdmol.RemoveAllConformers()
        return rdmol

    def to_XYZ_block(self) -> str:
        """
        Get the XYZ block.

        Returns:
            str: The XYZ block.
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

    def to_SDF_block(self, engine: Literal["rdkit", "openbabel"] = "rdkit") -> str:
        """
        Get the SDF block.

        Parameters:
            engine (Literal["rdkit", "openbabel"]):
                The engine to generate the SDF block. default is "rdkit".

        Returns:
            str: The SDF block.
        """
        if self.rdmol is None:
            moloplogger.warning(
                f"{self.file_path}: SDF building failed. No RDKit molecule found."
            )
            return ""
        if engine == "rdkit":
            return Chem.MolToMolBlock(self.rdmol)
        elif engine == "openbabel":
            return self.omol.write("sdf")
        else:
            raise ValueError(f"Unsupported engine: {engine}")

    def to_CML_block(self) -> str:
        """
        Get the CML block.

        Returns:
            str: The CML block.
        """
        if self.rdmol is None:
            moloplogger.warning(
                f"{self.file_path}: CML building failed. No RDKit molecule found."
            )
            return ""
        return Chem.MolToMrvBlock(self.rdmol)

    def to_SMILES(self) -> str:
        """
        Get the SMILES with explicit hydrogens.

        Returns:
            str: The SMILES.
        """
        if self.rdmol is None:
            return ""
        smi = Chem.MolToSmiles(self.rdmol)
        if Chem.MolFromSmiles(smi) is None:
            moloplogger.error(f"{self.file_path}: SMILES building failed.")
            return ""
        return smi

    def to_canonical_SMILES(self) -> str:
        """
        Get the SMILES with standardization.

        Returns:
            str: The standard SMILES.
        """
        return canonical_smiles(self.to_SMILES())

    def to_InChI(self) -> str:
        """
        Get the InChI.

        Returns:
            str: The InChI.
        """
        return Chem.MolToInchi(self.rdmol)

    def calc_rdkit_descs(self, desc_names: List[str] = None) -> Dict[str, float]:
        """
        Calculate the RDKit descriptors.

        Parameters:
            desc_names (List[str]):
                The names of the descriptors. Must be a subset of the RDKit descriptors.

        Returns:
            Dict[str,float]: The dictionary of the descriptors.
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
            desc_names (List[str]):
                The names of the descriptors. Must be a subset of "SOAP", "ACSF", "MBTR", "LMBTR"

        Returns:
            Dict[str,np.ndarray]: The dictionary of the descriptors.
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
            Dict[str,float]: The dictionary of the descriptors.
        """

        from molop.descriptor.descriptor import calc_mordred_descs

        return calc_mordred_descs(self.rdmol, **kwargs)

    def calc_spms_descriptor(
        self,
        anchor_list: Union[Sequence[int], None] = None,
        sphere_radius: Union[float, None] = None,
        atom_radius: Literal["vdw", "covalent"] = "vdw",
        latitudinal_resolution: int = 40,
        longitude_resolution: int = 40,
        precision: int = 8,
        *,
        custom_first_anchors: Union[Sequence[int], None] = None,
        custom_second_anchors: Union[Sequence[int], None] = None,
        custom_third_anchors: Union[Sequence[int], None] = None,
    ) -> SPMSCalculator:
        """
        Re-implementation of SPMS descriptor [A Molecular Stereostructure Descriptor based on
        Spherical Projection](https://www.thieme-connect.de/products/ejournals/abstract/10.1055/s-0040-1705977).
        GitHub repository: https://github.com/licheng-xu-echo/SPMS.git

        Parameters:
            anchor_list (Sequence[int] | None):
                The anchor atoms for SPMS calculation. If None, the default anchor atoms will be used.
            sphere_radius (float | None):
                The radius of the sphere for SPMS calculation. If None, the default sphere radius will be used.
            atom_radius:
                Atom radius type. Default is 'vdw', which means using Van der Waals radii as atom radii.
                If you want to use covalent radii, set this parameter to 'covalent'
            latitudinal_resolution (int):
                The number of latitudinal divisions for SPMS calculation.
            longitude_resolution (int):
                The number of longitudinal divisions for SPMS calculation.
            precision (int):
                The precision of the SPMS calculation.
            custom_first_anchors (Sequence[int] | None):
                The custom anchor atoms for SPMS calculation. If None, the default anchor atoms will be used.
            custom_second_anchors (Sequence[int] | None):
                The custom anchor atoms for SPMS calculation. If None, the default anchor atoms will be used.
            custom_third_anchors (Sequence[int] | None):
                The custom anchor atoms for SPMS calculation. If None, the default anchor atoms will be used.

        Returns:
            SPMSCalculator: The SPMS calculator object.
        """
        return SPMSCalculator(
            rdmol=self.rdmol,
            anchor_list=anchor_list,
            sphere_radius=sphere_radius,
            atom_radius=atom_radius,
            latitudinal_resolution=latitudinal_resolution,
            longitude_resolution=longitude_resolution,
            precision=precision,
            custom_first_anchors=custom_first_anchors,
            custom_second_anchors=custom_second_anchors,
            custom_third_anchors=custom_third_anchors,
        )

    def __hash__(self) -> int:
        return hash(str(self))

    def __len__(self) -> int:
        return len(self.atoms)

    def geometry_analysis(self, atom_idxs: Sequence[int], one_start=False) -> float:
        """
        Get the geometry infos among the atoms

        Parameters:
            atom_idxs (Sequence[int]):
                A Sequence of index of the atoms, starts from 0
            one_start (bool):
                If true, consider atom index starts from 1, so let index value subtracts 1 for all the atoms

        Returns:
            float:
                - If the length of atom_idxs is 2, the bond length with unit Angstrom between the two atoms will be returned.
                - If the length of atom_idxs is 3, the angle with unit degree between  the three atoms will be returned.
                - If the length of atom_idxs is 4, the dihedral angle with unit degree between the four atoms will be returned.
        """
        if one_start:
            atom_idxs = [atom_idx - 1 for atom_idx in atom_idxs]
        return get_geometry_info(self.rdmol, atom_idxs)

    def compare_rmsd(
        self,
        other: "BaseMolFrame",
        *,
        ignore_H: bool = False,
        centroid_align: bool = True,
        rotate_align: Literal["None", "kabsch", "quaternion"] = "kabsch",
        atom_idxs: Union[Sequence[int], None] = None,
    ) -> float:
        """
        Calculate the RMSD between two molecules.
        Based on [rmsd](https://github.com/charnley/rmsd) library.

        Parameters:
            other (BaseMolFrame):
                The other molecule to compare.
            ignore_H (bool):
                Whether to ignore the H atoms.
            centroid_align (bool):
                Whether to align the molecules by their centroids.
            rotate_align (Literal["None", "kabsch", "quaternion"]):
                Whether to align the molecules by their rotation matrices.
            atom_idxs (Union[Sequence[int], None]):
                A sequence of atom indices to calculate the RMSD.
                If None, all atoms will be used.

        Returns:
            float: The RMSD value.
        """
        assert isinstance(
            other, BaseMolFrame
        ), "The other object is not a BaseMolFrame object."
        return calculate_rmsd(
            self.rdmol,
            other.rdmol,
            ignore_H=ignore_H,
            centroid_align=centroid_align,
            rotate_align=rotate_align,
            atom_idxs=atom_idxs,
        )
