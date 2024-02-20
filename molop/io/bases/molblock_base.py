"""
Author: TMJ
Date: 2024-01-11 21:02:36
LastEditors: TMJ
LastEditTime: 2024-02-09 20:48:19
Description: 请填写简介
"""
import os
import re
from abc import ABC, abstractmethod
from typing import Any, Dict, List, Literal, Tuple, Union, Generator

import numpy as np
from openbabel import pybel
from pint.facets.plain import PlainQuantity
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds
from rdkit.Chem.MolStandardize import rdMolStandardize

from molop.logger.logger import logger
from molop.structure.geometry import get_geometry_info
from molop.structure.structure import (
    attempt_replacement,
    check_crowding,
    get_bond_pairs,
    get_formal_charges,
    get_formal_spins,
    reset_atom_index,
)
from molop.structure.structure_recovery import xyz_block_to_omol
from molop.unit import atom_ureg


class MolBlock(ABC):
    """
    Abstract class for a molecule block. The information is enough to recover the original molecule structure.

    Attributes:
        _atoms Union[List[str], List[int]]:
            **Need to be given**: Save the atoms. Consistent with the atomic order in the original file. Supports both element numbers and element symbols.
        _coords List[Tuple[PlainQuantity, PlainQuantity, PlainQuantity]]:
            **Need to be given**: Save the coordinates described by [pint](https://pint.readthedocs.io/en/stable/). Consistent with the atomic order in the original file. Defualt unit is angstrom.
        _charge int:
            **Need to be given**: Save the total charge. Allowed charge range is -3 to 3.
        _multiplicity int:
            **Need to be given**: Save the multiplicity. Allowed multiplicity range is 1 to 5.
        _bonds List[Tuple[int, int, int]]:
            **Predicted**: Save the bonds. The elements are: atom1 index, atom2 index, bond order(follow rdkit bond order `Chem.rdchem.BondType`).
        _formal_charges List[int]:
            **Predicted**: Save the formal charges. Consistent with the atomic order.
        _formal_spins List[int]:
            **Predicted**: Save the formal spins. Consistent with the atomic order.
        _omol pybel.Molecule:
            **Predicted**: Save the openbabel molecule object.
        _rdmol Union[Chem.rdchem.RWMol, Chem.rdchem.Mol]:
            **Predicted**: Save the rdkit molecule object.
    """

    _atoms: Union[List[str], List[int]]
    _coords: List[Tuple[PlainQuantity, PlainQuantity, PlainQuantity]]
    _charge: int
    _multiplicity: int
    _bonds: List[Tuple[int, int, int]]
    _formal_charges: List[int]
    _formal_spins: List[int]
    _omol: pybel.Molecule
    _rdmol: Union[Chem.rdchem.RWMol, Chem.rdchem.Mol]

    def __init__(self):
        self._atoms: Union[List[str], List[int]] = []
        self._coords: List[Tuple[float]] = []
        self._charge: int = 0
        self._multiplicity: int = 1
        self._bonds: List[Tuple[int, int, int]] = None
        self._formal_charges: List[int] = None
        self._formal_spins: List[int] = None
        self._omol: pybel.Molecule = None
        self._rdmol: Union[Chem.rdchem.RWMol, Chem.rdchem.Mol] = None

    @property
    def atoms(self) -> List[str]:
        """
        Get the atoms.

        Returns:
            A list of element symbols of the atoms.
        """
        return [Chem.Atom(atom).GetSymbol() for atom in self._atoms]

    @property
    def total_electrons(self) -> int:
        """
        Get the total electrons.

        Returns:
            The total electrons.
        """
        return sum(Chem.Atom(atom).GetAtomicNum() for atom in self._atoms) + self.charge

    @property
    def elements(self) -> List[str]:
        """
        Get the elements set.

        Returns:
            A list of element symbols dropped duplicates.
        """
        return list(set(self.atoms))

    @property
    def coords(self) -> List[Tuple[PlainQuantity, PlainQuantity, PlainQuantity]]:
        """
        Get the coordinates.

        Returns:
            A list of coordinates. Each coordinate is a tuple of three plain quantities. Default unit is angstrom.
        """
        return self._coords

    @property
    def dimensionless_coords(self) -> List[Tuple[float, float, float]]:
        """
        Get the dimensionless coordinates.

        Units will be tranformed:
            - `angstrom`

        Returns:
            A list of dimensionless coordinates
        """
        return [
            tuple([coord.to("angstrom").m for coord in atom]) for atom in self.coords
        ]

    @property
    def charge(self) -> int:
        """
        Get the charge.

        Returns:
            The total charge. Allowed charge range is -3 to 3.
        """
        return self._charge

    @property
    def multiplicity(self) -> int:
        """
        Get the multiplicity.

        Returns:
            The multiplicity. Allowed multiplicity range is 1 to 5.
        """
        return self._multiplicity

    @property
    def bonds(self) -> List[Tuple[int, int, int]]:
        """
        Get the bonds.

        Returns:
            A list of bonds. Each bond is a tuple of three integers. The first two integers are the atom indices of the bond. The third integer is the bond order(follow rdkit bond order `Chem.rdchem.BondType`).
        """
        if self._bonds is None:
            self._bonds = get_bond_pairs(self.rdmol)
            self._formal_charges = get_formal_charges(self.rdmol)
            self._formal_spins = get_formal_spins(self.rdmol)
        return self._bonds

    @property
    def formal_charges(self) -> List[int]:
        """
        Get the formal charges.

        Returns:
            A list of formal charges. Consistent with the atomic order.
        """
        if self._formal_charges is None:
            self._bonds = get_bond_pairs(self.rdmol)
            self._formal_charges = get_formal_charges(self.rdmol)
            self._formal_spins = get_formal_spins(self.rdmol)
        return self._formal_charges

    @property
    def formal_spins(self) -> List[int]:
        """
        Get the formal spins.

        Returns:
            A list of formal spins. Consistent with the atomic order.
        """
        if self._formal_spins is None:
            self._bonds = get_bond_pairs(self.rdmol)
            self._formal_charges = get_formal_charges(self.rdmol)
            self._formal_spins = get_formal_spins(self.rdmol)
        return self._formal_spins

    @property
    def omol(self) -> pybel.Molecule:
        """
        Get the openbabel molecule object.

        Returns:
            The openbabel molecule object.
        """
        if self._omol is None:
            if self._bonds is None:
                try:
                    self._omol = xyz_block_to_omol(
                        self.to_XYZ_block(),
                        self._charge,
                    )
                except Exception as e:
                    raise RuntimeError(f"{self._file_path}: {e}")
            else:
                assert (
                    self._formal_charges and self._formal_spins
                ), "If bonds given, formal charges and spins must be provided."
                omol = pybel.readstring("sdf", self.to_SDF_block())
                self._omol = omol
        return self._omol

    @property
    def rdmol(self) -> Union[Chem.rdchem.Mol, Chem.rdchem.RWMol]:
        """
        Get the rdkit molecule object.

        Returns:
            The rdkit molecule object.
        """
        if self._rdmol is None:
            if self._bonds is None:
                try:
                    # try rdkit determinebonds first
                    # Issues known:
                    # - Can not recognize radicals
                    # - Can not recognize Metals
                    raw_mol = Chem.MolFromXYZBlock(self.to_XYZ_block())
                    conn_mol = Chem.Mol(raw_mol)
                    rdDetermineBonds.DetermineBonds(conn_mol, charge=self._charge)
                    self._rdmol = conn_mol
                    self._bonds = get_bond_pairs(self._rdmol)
                    self._formal_charges = get_formal_charges(self._rdmol)
                    self._formal_spins = get_formal_spins(self._rdmol)
                except:
                    logger.debug(
                        f"{self._file_path}: rdkit determinebonds failed. Use MolOP structure recovery instead."
                    )
                    # If failed, use MolOP implementation
                    self._rdmol = Chem.MolFromMolBlock(
                        self.omol.write("sdf"), removeHs=False
                    )
                    self._bonds = get_bond_pairs(self._rdmol)
                    self._formal_charges = get_formal_charges(self._rdmol)
                    self._formal_spins = get_formal_spins(self._rdmol)
            else:
                assert (
                    self._formal_charges and self._formal_spins
                ), "If bonds given, formal charges and spins must be provided."
                rwmol = Chem.RWMol(Chem.MolFromXYZBlock(self.to_XYZ_block()))
                for bond in self._bonds:
                    rwmol.AddBond(bond[0], bond[1], bond[2])
                for atom, charge, spin in zip(
                    rwmol.GetAtoms(), self._formal_charges, self._formal_spins
                ):
                    atom.SetFormalCharge(charge)
                    atom.SetNumRadicalElectrons(spin)
                self._rdmol = rwmol

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

    def basic_check(self) -> bool:
        assert self._charge in [-3, -2, -1, 0, 1, 2, 3], "The charge must be -3 ~ +3."
        assert self._multiplicity in [1, 2, 3, 4, 5], "The multiplicity must be 1 ~ 5."
        assert len(self.atoms) == len(
            self.coords
        ), "The number of atoms and coordinates are not equal."
        for atom, coord in zip(self.atoms, self.coords):
            assert isinstance(
                atom, Union[str, int]
            ), "The atom must be a string or an integer."
            assert len(coord) == 3, "The coordinate must be a list of three floats."
        if self.bonds is not None:
            for bond in self.bonds:
                assert len(bond) == 3, "The bond must be a list of three integers."
                assert all(
                    isinstance(atom, Union[str, int]) for atom in bond
                ), "The atom must be a string or an integer."

    def to_XYZ_block(self) -> str:
        """
        Get the XYZ block.

        Returns:
            The XYZ block.
        """
        return (
            f"{len(self.atoms)}\n"
            + f"charge {self.charge} multiplicity {self.multiplicity}\n"
            + "\n".join(
                [
                    f"{atom:10s}{x.m:10.5f}{y.m:10.5f}{z.m:10.5f}"
                    for atom, x, y, z in zip(self.atoms, *zip(*self.coords))
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

    def to_SMILES(self) -> str:
        """
        Get the SMILES with explicit hydrogens.

        Returns:
            The SMILES.
        """

        smiles = Chem.MolToSmiles(self.rdmol)
        return smiles

    def to_standard_SMILES(self) -> str:
        """
        Get the SMILES with standardization.

        Returns:
            The SMILES.
        """
        return rdMolStandardize.StandardizeSmiles(self.to_SMILES())

    def to_InChI(self) -> str:
        """
        Get the InChI.

        Returns:
            The InChI.
        """
        return Chem.MolToInchi(self.rdmol)

    def to_XYZ_file(self, file_path: str = None) -> str:
        """
        Write the XYZ file.

        Parameters:
            file_path str:
                The file path. If not specified, will be generated in situ.
        Returns:
            The absolute path of the XYZ file.
        """
        if file_path is None:
            file_path = self._file_path
        if os.path.isdir(file_path):
            raise IsADirectoryError(f"{file_path} is a directory.")
        file_path = os.path.splitext(file_path)[0] + ".xyz"
        with open(file_path, "w") as f:
            f.write(self.to_XYZ_block())
        f.close()
        return os.path.abspath(file_path)

    def to_SDF_file(self, file_path: str = None) -> str:
        """
        Write the SDF file.

        Parameters:
            file_path str:
                The file path. If not specified, will be generated in situ.
        Returns:
            The absolute path of the SDF file.
        """
        if file_path is None:
            file_path = self._file_path
        if os.path.isdir(file_path):
            raise IsADirectoryError(f"{file_path} is a directory.")
        file_path = os.path.splitext(file_path)[0] + ".sdf"
        with open(file_path, "w") as f:
            f.write(self.to_SDF_block())
        f.close()
        return os.path.abspath(file_path)

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

    @abstractmethod
    def __str__(self) -> str:
        raise NotImplementedError

    def __hash__(self) -> int:
        return hash(str(self))

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({str(self)})"

    def __len__(self) -> int:
        return len(self.atoms)

    def to_GJF_block(
        self,
        charge: int = None,
        multiplicity: int = None,
        template: str = None,
        prefix: str = "# g16 gjf",
        suffix="",
    ) -> str:
        """
        Get the GJF block.

        Parameters:
            charge int:
                The forced charge. If specified, will be used to overwrite the charge in the gjf file.
            multiplicity int:
                The forced multiplicity. If specified, will be used to overwrite the multiplicity in the gjf file.
            template str:
                path to read a gjf file as a template.
            prefix str:
                prefix to add to the beginning of the gjf file, priority is lower than template.
            suffix str:
                suffix to add to the end of the gjf file, priority is lower than template.
        Returns:
            A modified GJF block.
        """
        if template is not None:
            if not os.path.isfile(template):
                raise FileNotFoundError(f"{template} is not found.")
            with open(template, "r") as f:
                lines = f.readlines()
            f.close()
            for i, line in enumerate(lines):
                if re.match(r"^\s*[\+\-\d]+\s+\d+$", line):
                    prefix = "".join(lines[: i - 2])
                if re.match(r"^\s*[A-Z][a-z]?(\s+\-?\d+(\.\d+)?){3}$", line):
                    suffix = "".join(lines[i + 1 :])
        else:
            prefix = prefix if prefix.endswith("\n") else prefix + "\n"
            prefix = prefix + "\n"
        return (
            prefix
            + f" Title: Generated by MolOP\n\n"
            + f"{charge if charge else self.charge} {multiplicity if multiplicity else self.multiplicity}\n"
            + "\n".join(
                [
                    f"{atom:10s}{x.m:10.5f}{y.m:10.5f}{z.m:10.5f}"
                    for atom, x, y, z in zip(self.atoms, *zip(*self.coords))
                ]
            )
            + "\n\n"
            + suffix
            + "\n\n"
        )

    def to_GJF_file(
        self,
        file_path: str = None,
        charge: int = None,
        multiplicity: int = None,
        template: str = None,
        prefix: str = "# g16 gjf \n",
        suffix="",
    ) -> str:
        """
        Write the GJF file.

        Parameters:
            file_path str:
                The path to write the GJF file. If not specified, will be generated in situ.
            charge int:
                The forced charge. If specified, will be used to overwrite the charge in the gjf file.
            multiplicity int:
                The forced multiplicity. If specified, will be used to overwrite the multiplicity in the gjf file.
            template str:
                path to read a gjf file as a template.
            prefix str:
                prefix to add to the beginning of the gjf file, priority is lower than template.
            suffix str:
                suffix to add to the end of the gjf file, priority is lower than template.
        Returns:
            The path to the GJF file.
        """
        if file_path is None:
            file_path = self._file_path
        if os.path.isdir(file_path):
            raise IsADirectoryError(f"{file_path} is a directory.")
        file_path = os.path.splitext(file_path)[0] + ".gjf"
        with open(file_path, "w") as f:
            f.write(
                self.to_GJF_block(
                    charge=charge,
                    multiplicity=multiplicity,
                    template=template,
                    prefix=prefix,
                    suffix=suffix,
                )
            )
        f.close()
        return os.path.abspath(file_path)

    def to_chemdraw(self, file_path: str = None, keep3D=True):
        """
        Write the ChemDraw file.

        Parameters:
            file_path str:
                The path to write the ChemDraw file. If not specified, will be generated in situ.
            keep3D bool:
                Whether to keep the 3D information.
        Returns:
            The path to the ChemDraw file.
        """
        if file_path is None:
            file_path = self._file_path
        if os.path.isdir(file_path):
            raise IsADirectoryError(f"{file_path} is a directory.")
        file_path = os.path.splitext(file_path)[0] + ".cdxml"
        if not keep3D:
            temp_rdmol = Chem.RWMol(self.rdmol)
            temp_rdmol.RemoveAllConformers()
            pybel.readstring("sdf", Chem.MolToMolBlock(temp_rdmol)).write(
                "cdxml", file_path, overwrite=True
            )
        else:
            self.omol.write("cdxml", file_path, overwrite=True)
        return os.path.abspath(file_path)

    def geometry_analysis(self, atom_idxs: Tuple[int], one_start=False) -> float:
        """
        Get the geometry infos among the atoms

        Parameters:
            atom_idxs Tuple[int]:
                A list of index of the atoms, starts from 0
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


class BaseBlockParser(MolBlock):
    """
    Base class for block parsers.
    Only need basic information.

    Attributes:
        _block_type str:
            The type of the block. eg. ".xyz", ".log", et. al.
        _file_path str:
            The absolute path to the file.
        _frameID int:
            The frameID of the block.
        _block str:
            The block content.
        _next_block BaseBlockParser:
            The next block. For calling as link list.
    """

    _block_type: str
    _file_path: str
    _frameID: int
    _block: str
    _next_block: "BaseBlockParser"

    def __init__(self, block: str) -> None:
        super().__init__()
        self._file_path = None
        self._frameID: int = 0
        self._block = block
        self._next_block = None

    @property
    def block(self) -> str:
        """
        Get the block content.

        Returns:
            The block content.
        """
        return self._block

    def __str__(self) -> str:
        return (
            f"{os.path.basename(self._file_path)}[{self._frameID}]\n"
            + self.to_XYZ_block()
        )

    def next(self) -> "BaseBlockParser":
        """
        Get the next block.

        Returns:
            The next block.
        """
        return self._next_block

    def to_dict(self) -> Dict[str, Any]:
        """
        Get the information of the block as a dict.

        Contents:
            - type: str
            - file_path: str
            - block: str
            - SMILES: str
            - frameID: int
            - atoms: List[str]
            - coords: List[Tuple[float, float, float]]
            - bonds: List[Tuple[int, int, int]]
            - total charge: int
            - total multiplicity: int
            - formal charges: List[int]

        Returns:
            The information of the block as a dict.
        """
        return {
            "type": self._block_type,
            "file_path": self._file_path,
            "block": self.block,
            "SMILES": self.to_SMILES(),
            "frameID": self._frameID,
            "atoms": self.atoms,
            "coords": self.dimensionless_coords,
            "bonds": self.bonds,
            "total charge": self.charge,
            "total multiplicity": self.multiplicity,
            "formal charges": self.formal_charges,
            "formal spins": self.formal_spins,
        }

    def summary(self):
        """
        Print the information of the block.
        """
        print(
            f"type: {self._block_type}\n"
            + f"file path: {self._file_path}\n"
            + f"frameID: {self._frameID}\n"
            + f"SMILES: {self.to_SMILES()}\n"
            + f"atom number: {len(self)}\n"
            + f"total charge: {self.charge}\n"
        )

    def rebuild_parser(
        self,
        new_mol: Union[Chem.rdchem.Mol, Chem.rdchem.RWMol],
        rebuild_type=Literal["mod", "reindex"],
    ) -> "BaseBlockParser":
        """
        Create a new parser with the given molecule.

        Parameters:
            new_mol Union[Chem.rdchem.Mol, Chem.rdchem.RWMol]:
                The new molecule.
            rebuild_type Literal["mod", "reindex"]:
                The tag of the new parser.
        Returns:
            The new parser.
        """
        new_parser = BaseBlockParser(block=Chem.MolToXYZBlock(new_mol))
        new_parser._atoms = [atom.GetSymbol() for atom in new_mol.GetAtoms()]
        new_parser._coords = [
            (x * atom_ureg.angstrom, y * atom_ureg.angstrom, z * atom_ureg.angstrom)
            for x, y, z in new_mol.GetConformer().GetPositions()
        ]
        new_parser._file_path = (
            os.path.splitext(self._file_path)[0] + f"_{rebuild_type}.xyz"
        )
        new_parser._rdmol = new_mol
        new_parser._charge = self.charge
        new_parser._multiplicity = self.multiplicity
        new_parser._bonds = get_bond_pairs(new_mol)
        new_parser._formal_charges = get_formal_charges(new_mol)
        new_parser._formal_spins = get_formal_spins(new_mol)
        return new_parser

    def replace_substituent(
        self,
        query_smi: str,
        replacement_smi: str,
        bind_idx: int = None,
        replace_all=False,
        attempt_num: int = 10,
    ) -> "BaseBlockParser":
        """
        Replace the substituent with the given SMARTS. The substituent is defined by the query_smi, and the new substituent is defined by the replacement_smi.

        Parameters:
            query_smi str:
                The SMARTS to query the substituent in the original molecule.
            replacement_smi str:
                The SMARTS of new substituent. The bind atom is the first atom of the replacement_smi.
            bind_idx int:
                The index of the atom to bind the new substituent. The default is None, which means to replace the first legal atom in original molecule.
                If specified, try to replace the atom. User should meke sure the atom is legal.
                Detail example in (Repalce Substituent)[Repalce Substituent]
            replace_all bool:
                If True, replace all the substituent queried in the original molecule.
            attempt_num int:
                Max attempt times to replace the substituent. Each time a new substituent conformation will be used for substitution.
        Returns:
            The new parser.
        """
        new_mol = attempt_replacement(
            self.rdmol,
            query_smi=query_smi,
            replacement_smi=replacement_smi,
            bind_idx=bind_idx,
            replace_all=replace_all,
            attempt_num=attempt_num,
        )
        return self.rebuild_parser(new_mol, rebuild_type="mod")

    def reset_atom_index(self, mapping_smarts: str) -> "BaseBlockParser":
        """
        Reset the atom index of the molecule according to the mapping SMARTS.

        Parameters:
            mapping_smarts str:
                The SMARTS to query the molecule substructure.
                The queried atoms will be renumbered and placed at the beginning of all atoms according to the order of the atoms in SMARTS. The relative order of the remaining atoms remains unchanged.

        Returns:
            The new parser.
        """
        smarts = Chem.MolFromSmarts(mapping_smarts)
        if not self.rdmol.HasSubstructMatch(smarts):
            logger.error(f"Failed to match {self._file_path} with {mapping_smarts}")
            raise ValueError("Failed to match")
        mapping = self.rdmol.GetSubstructMatches(smarts)
        if len(mapping) > 1:
            logger.warning(
                f"Multiple matches found in {self._file_path} with {mapping_smarts}"
            )
        mapping = mapping[0]
        rdmol = reset_atom_index(self.rdmol, mapping)
        return self.rebuild_parser(
            rdmol,
            rebuild_type="reindex",
        )


class QMBaseBlockParser(BaseBlockParser):
    """
    Base class for QM block parsers.
    Should fill some the QM information pre-defined.

    Attributes:
        _version str:
            The version of the QM software.
        _parameter_comment str:
            The comment of the QM parameters. Like calculation level, basis set, etc.
        _energy PlainQuantity:
            The total energy of the molecule.
        _partial_charges List[float]:
            The partial charges of the molecule.
        _spin_densities List[float]:
            The spin densities of the molecule.
        _spin_multiplicity float:
            The spin multiplicity of the molecule.
        _spin_eigenvalue float:
            The spin eigenvalue of the molecule. spin egigenvalue = sqrt(spin multiplicity * (spin multiplicity + 1)).
        _gradients List[Tuple[PlainQuantity, PlainQuantity, PlainQuantity]]:
            The gradients (forces) of the molecule.
        _frequencies List[Dict[str, Union[bool, PlainQuantity, List[Tuple[PlainQuantity]]]]]:
            The frequencies of the molecule. The frequency dict has the following keys:

                - is imaginary: bool
                - freq: PlainQuantity
                - reduced masses: PlainQuantity
                - IR intensities: PlainQuantity
                - force constants: PlainQuantity
                - normal coordinates: List[Tuple[PlainQuantity, PlainQuantity, PlainQuantity]]
        _alpha_FMO_orbits List[PlainQuantity]:
            The alpha FMO orbitals of the molecule.
        _beta_FMO_orbits List[PlainQuantity]:
            The beta FMO orbitals of the molecule.
        _alpha_energy Dict[str, PlainQuantity]:
            The alpha energy of the molecule. The energy dict has the following keys:

                - gap: PlainQuantity
                - homo: PlainQuantity
                - lumo: PlainQuantity
        _beta_energy Dict[str, PlainQuantity]:
            The beta energy of the molecule. The energy dict has the following keys:

                - gap: PlainQuantity
                - homo: PlainQuantity
                - lumo: PlainQuantity
        _sum_energy Dict[str, PlainQuantity]:
            The sum energy of the molecule. The energy dict has the following keys:

                - zero-point gas: PlainQuantity
                - E gas: PlainQuantity
                - H gas: PlainQuantity
                - G gas: PlainQuantity
                - zero-point correction: PlainQuantity
                - TCE: PlainQuantity
                - TCH: PlainQuantity
                - TCG: PlainQuantity
        _state Dict[str, bool]:
            The state of the molecule.
        _only_extract_structure bool:
            A flag to indicate whether to only extract the structure.
    """

    _version: str
    _parameter_comment: str
    _energy: PlainQuantity
    _partial_charges: List[float]
    _spin_densities: List[float]
    _gradients: List[Tuple[PlainQuantity, PlainQuantity, PlainQuantity]]
    _spin_multiplicity: float
    _spin_eigenvalue: float
    # _hessian: List[List[float]]
    _frequencies: List[
        Dict[
            Literal[
                "is imaginary",
                "freq",
                "reduced masses",
                "IR intensities",
                "force constants",
                "normal coordinates",
            ],
            Union[
                bool,
                PlainQuantity,
                List[Tuple[PlainQuantity, PlainQuantity, PlainQuantity]],
            ],
        ],
    ]
    _alpha_FMO_orbits: List[PlainQuantity]
    _beta_FMO_orbits: List[PlainQuantity]
    _alpha_energy: Dict[Literal["gap", "homo", "lumo"], PlainQuantity]
    _beta_energy: Dict[Literal["gap", "homo", "lumo"], PlainQuantity]
    _sum_energy: Dict[
        Literal[
            "zero-point gas",
            "E gas",
            "H gas",
            "G gas",
            "zero-point correction",
            "TCE",
            "TCH",
            "TCG",
        ],
        PlainQuantity,
    ]
    # _nbo_analysis: List[Dict[str, Dict[str, PlainQuantity]]]

    _state: Dict[str, bool]
    _only_extract_structure: bool

    def __init__(self, block: str, only_extract_structure=False) -> None:
        super().__init__(block)
        self._only_extract_structure = only_extract_structure

        self._version: str = None
        self._parameter_comment: str = None
        self._energy: PlainQuantity = None
        self._partial_charges: List[float] = []
        self._spin_densities: List[float] = []
        self._spin_multiplicity: float = None
        self._spin_eigenvalue: float = None
        self._gradients: List[Tuple[PlainQuantity, PlainQuantity, PlainQuantity]] = []
        # self._hessian: List[List[float]] = []
        self._alpha_FMO_orbits: List[PlainQuantity] = []
        self._beta_FMO_orbits: List[PlainQuantity] = []
        self._alpha_energy = {
            "gap": None,
            "homo": None,
            "lumo": None,
        }
        self._beta_energy = {
            "gap": None,
            "homo": None,
            "lumo": None,
        }
        self._sum_energy = {
            "zero-point gas": None,
            "E gas": None,
            "H gas": None,
            "G gas": None,
            "zero-point correction": None,
            "TCE": None,
            "TCH": None,
            "TCG": None,
        }
        self._frequencies = []
        # self._nbo_analysis = []
        self._state = {}

    @property
    def energy(self) -> PlainQuantity:
        """
        Get the total energy of the molecule.

        Returns:
            The total energy of the molecule.
        """
        return self._energy

    @property
    def dimensionless_energy(self) -> float:
        """
        Get the dimensionless total energy of the molecule.

        Unit will be transformed:
            - `hartree/particle`

        Returns:
            The dimensionless total energy of the molecule.
        """
        if self.energy is None:
            return None
        return self.energy.to("hartree/particle").m

    @property
    def partial_charges(self) -> List[float]:
        """
        Get the partial charges of the molecule.

        Returns:
            The partial charges of the molecule.
        """
        return self._partial_charges

    @property
    def spin_densities(self) -> List[float]:
        """
        Get the spin densities of the molecule.

        Returns:
            The spin densities of the molecule.
        """
        return self._spin_densities

    @property
    def spin_multiplicity(self) -> float:
        """
        Get the spin multiplicity of the molecule.

        Returns:
            The spin multiplicity of the molecule.
        """
        return self._spin_multiplicity

    @property
    def spin_eigenvalue(self) -> float:
        """
        Get the spin eigenvalue of the molecule.

        Returns:
            The spin eigenvalue of the molecule.
        """
        return self._spin_eigenvalue

    @property
    def gradients(self) -> List[Tuple[PlainQuantity, PlainQuantity, PlainQuantity]]:
        """
        Get the gradients of the molecule.

        Returns:
            The gradients of the molecule.
        """
        return self._gradients

    @property
    def dimensionless_gradients(
        self,
    ) -> Generator[Tuple[float, float, float], None, None]:
        """
        Get the dimensionless gradients of the molecule.

        Unit will be transformed:
            - `hartree/bohr`

        Returns:
            The dimensionless gradients of the molecule.
        """
        return (
            tuple([gradient.to("hartree/bohr").m for gradient in atom])
            for atom in self.gradients
        )

    # @property
    # def hessian(self) -> List[float]:
    # return self._hessian

    @property
    def is_TS(self) -> bool:
        """
        Check if the molecule is a TS.

        Returns:
            True if the molecule is a TS, False otherwise.
        """
        return len(self.imaginary_frequencies) == 1

    @property
    def imaginary_frequencies(
        self,
    ) -> List[
        Dict[
            Literal[
                "is imaginary",
                "freq",
                "reduced masses",
                "IR intensities",
                "force constants",
                "normal coordinates",
            ],
            Union[
                bool,
                PlainQuantity,
                List[Tuple[PlainQuantity, PlainQuantity, PlainQuantity]],
            ],
        ]
    ]:
        """
        Get the imaginary frequencies of the molecule.

        Returns:
            The imaginary frequencies of the molecule.
        """
        return [freq for freq in self._frequencies if freq["is imaginary"]]

    @property
    def dimensionless_imaginary_frequencies(
        self,
    ) -> Generator[
        Dict[
            Literal[
                "is imaginary",
                "freq",
                "reduced masses",
                "IR intensities",
                "force constants",
                "normal coordinates",
            ],
            Union[
                bool,
                float,
                List[Tuple[float, float, float]],
            ],
        ],
        None,
        None,
    ]:
        """
        Get the dimensionless imaginary frequencies of the molecule.

        Unit will be transformed:
            - frequency in `cm^-1`
            - reduced masses in `amu`
            - IR intensities in `kmol/mol`
            - force constants in `mdyne/angstrom`
            - normal coordinates in `angstrom`

        Returns:
            The dimensionless imaginary frequencies of the molecule.
        """
        return (freq for freq in self.dimensionless_frequencies if freq["is imaginary"])

    @property
    def frequencies(
        self,
    ) -> List[
        Dict[
            Literal[
                "is imaginary",
                "freq",
                "reduced masses",
                "IR intensities",
                "force constants",
                "normal coordinates",
            ],
            Union[
                bool,
                PlainQuantity,
                List[Tuple[PlainQuantity, PlainQuantity, PlainQuantity]],
            ],
        ]
    ]:
        """
        Get the frequencies of the molecule.

        Returns:
            The frequencies of the molecule.
        """
        return self._frequencies

    @property
    def dimensionless_frequencies(
        self,
    ) -> Generator[
        Dict[
            Literal[
                "is imaginary",
                "freq",
                "reduced masses",
                "IR intensities",
                "force constants",
                "normal coordinates",
            ],
            Union[
                bool,
                float,
                List[Tuple[float, float, float]],
            ],
        ],
        None,
        None,
    ]:
        """
        Get the dimensionless frequencies of the molecule.

        Unit will be transformed:
            - frequency in `cm^-1`
            - reduced masses in `amu`
            - IR intensities in `kmol/mol`
            - force constants in `mdyne/angstrom`
            - normal coordinates in `angstrom`

        Returns:
            The dimensionless frequencies of the molecule.
        """
        return (
            {
                "freq": freq["freq"].to("cm^-1").m,
                "is imaginary": freq["is imaginary"],
                "reduced masses": freq["reduced masses"].to("amu").m,
                "IR intensities": freq["IR intensities"].to("kmol/mol").m,
                "force constants": freq["force constants"].to("mdyne/angstrom").m,
                "normal coordinates": [
                    tuple([coord.to("angstrom").m for coord in atom])
                    for atom in freq["normal coordinates"]
                ],
            }
            for freq in self._frequencies
        )

    @property
    def alpha_FMO_orbits(self) -> List[PlainQuantity]:
        """
        Get the alpha FMO orbitals of the molecule.

        Returns:
            The alpha FMO orbitals of the molecule.
        """
        return self._alpha_FMO_orbits

    @property
    def dimensionless_alpha_FMO_orbits(self) -> Generator[float, None, None]:
        """
        Get the dimensionless alpha FMO orbitals of the molecule.

        Unit will be transformed:
            - `hartree/particle`

        Returns:
            The dimensionless alpha FMO orbitals of the molecule.
        """
        return (orbit.to("hartree/particle").m for orbit in self.alpha_FMO_orbits)

    @property
    def alpha_energy(self) -> Dict[Literal["gap", "homo", "lumo"], PlainQuantity]:
        """
        Get the alpha energy of the molecule.

        Returns:
            The alpha energy of the molecule.
        """
        return self._alpha_energy

    @property
    def dimensionless_alpha_energy(self) -> Dict[Literal["gap", "homo", "lumo"], float]:
        """
        Get the dimensionless alpha energy of the molecule.

        Unit will be transformed:
            - `hartree/particle`

        Returns:
            The dimensionless alpha energy of the molecule.
        """
        return {
            key: value.to("hartree/particle").m if value is not None else None
            for key, value in self.alpha_energy.items()
        }

    @property
    def beta_FMO_orbits(self) -> List[PlainQuantity]:
        """
        Get the beta FMO orbitals of the molecule.

        Returns:
            The beta FMO orbitals of the molecule.
        """
        return self._beta_FMO_orbits

    @property
    def dimensionless_beta_FMO_orbits(self) -> Generator[float, None, None]:
        """
        Get the dimensionless beta FMO orbitals of the molecule.

        Unit will be transformed:
            - `hartree/particle`

        Returns:
            The dimensionless beta FMO orbitals of the molecule.
        """
        return (orbit.to("hartree/particle").m for orbit in self.beta_FMO_orbits)

    @property
    def beta_energy(self) -> Dict[Literal["gap", "homo", "lumo"], PlainQuantity]:
        """
        Get the beta energy of the molecule.

        Returns:
            The beta energy of the molecule.
        """
        return self._beta_energy

    @property
    def dimensionless_beta_energy(self) -> Dict[Literal["gap", "homo", "lumo"], float]:
        """
        Get the dimensionless beta energy of the molecule.

        Unit will be transformed:
            - `hartree/particle`

        Returns:
            The dimensionless beta energy of the molecule.
        """
        return {
            key: value.to("hartree/particle").m if value is not None else None
            for key, value in self.beta_energy.items()
        }

    @property
    def sum_energy(
        self,
    ) -> Dict[
        Literal[
            "zero-point gas",
            "E gas",
            "H gas",
            "G gas",
            "zero-point correction",
            "TCE",
            "TCH",
            "TCG",
        ],
        PlainQuantity,
    ]:
        """
        Get the sum energy of the molecule.

        Returns:
            The sum energy of the molecule.
        """
        return self._sum_energy

    @property
    def dimensionless_sum_energy(
        self,
    ) -> Dict[
        Literal[
            "zero-point gas",
            "E gas",
            "H gas",
            "G gas",
            "zero-point correction",
            "TCE",
            "TCH",
            "TCG",
        ],
        PlainQuantity,
    ]:
        """
        Get the dimensionless sum energy of the molecule.

        Unit will be transformed:
            - `hartree/particle`

        Returns:
            The dimensionless sum energy of the molecule.
        """
        return {
            key: value.to("hartree/particle").m if value is not None else None
            for key, value in self.sum_energy.items()
        }

    @property
    def state(self) -> Dict[str, bool]:
        """
        Get the state of the calculation.

        Returns:
            The state of the calculation.
        """
        return self._state

    @property
    def parameter_comment(self) -> str:
        """
        Get the parameter comment of the calculation.

        Returns:
            The parameter comment of the calculation.
        """
        return self._parameter_comment

    @property
    def version(self) -> str:
        """
        Get the version of the QM software.

        Returns:
            The version of the QM software.
        """
        return self._version

    # @property
    # def nbo_analysis(self):
    # return self._nbo_analysis

    def to_dict(self):
        """
        Get the dictionary representation of the molecule.

        Contents:
            - type: str
            - file_path: str
            - block: str
            - SMILES: str
            - frameID: int
            - atoms: List[str]
            - coords: List[Tuple[float, float, float]]
            - bonds: List[Tuple[int, int, int]]
            - total charge: int
            - total multiplicity: int
            - formal charges: List[int]
            - version: str
            - parameter_comment: str
            - state: Dict[str, bool]
            - energy: float
            - sum_energy: Dict[str, float]
            - gradients: List[Tuple[float, float, float]]
            - frequencies: List[Dict[Literal['is imaginary', 'freq', 'reduced masses', 'IR intensities', 'force constants', 'normal coordinates'], bool | float | List[Tuple[float, float, float]]]]
            - imaginary_frequencies: List[Dict[Literal['is imaginary', 'freq', 'reduced masses', 'IR intensities', 'force constants', 'normal coordinates'], bool | float | List[Tuple[float, float, float]]]]
            - partial_charges: List[float]
            - spin_densities: List[float]
            - alpha_FMO_orbits: List[float]
            - alpha_energy: Dict[str, float]
            - beta_FMO_orbits: List[float]
            - beta_energy: Dict[str, float]
            - spin_eginvalue: float
            - spin_multiplicity: float

        Returns:
            The dictionary representation of the molecule.
        """
        return {
            **super().to_dict(),
            **{
                "version": self.version,
                "parameter_comment": self.parameter_comment,
                "state": self.state,
                "energy": self.dimensionless_energy,
                "sum_energy": self.dimensionless_sum_energy,
                "gradients": list(self.dimensionless_gradients),
                "frequencies": list(self.dimensionless_frequencies),
                "imaginary_frequencies": list(self.dimensionless_imaginary_frequencies),
                "partial_charges": self.partial_charges,
                "spin_densities": self.spin_densities,
                "alpha_FMO_orbits": list(self.dimensionless_alpha_FMO_orbits),
                "alpha_energy": self.dimensionless_alpha_energy,
                "beta_FMO_orbits": list(self.dimensionless_beta_FMO_orbits),
                "beta_energy": self.dimensionless_beta_energy,
                "spin_eginvalue": self.spin_eigenvalue,
                "spin_multiplicity": self.spin_multiplicity,
                # "nbo_analysis": self.nbo_analysis,
                # "hessian": self.hessian,
            },
        }

    def summary(self):
        """
        Print the summary of the molecule.
        """
        print(
            f"type: {self._block_type}\n"
            + f"file path: {self._file_path}\n"
            + f"frameID: {self._frameID}\n"
            + f"SMILES: {self.to_SMILES()}\n"
            + f"atom number: {len(self)}\n"
            + f"total charge: {self.charge}\n"
            + f"version: {self.version}\n"
            + f"parameter comment: \n{self.parameter_comment}\n"
        )
        if not self._only_extract_structure:
            print(
                f"state: {self.state}\n"
                + f"energy: {self.energy}\n"
                + f"sum energy: {self.sum_energy}\n"
                + f"gradients number: {len(self.gradients)}\n"
                + f"frequencies number: {len(self.frequencies)}\n"
                + f"imaginary frequencies number: {len(self.imaginary_frequencies)}\n"
                + f"partial charges number: {len(self.partial_charges)}\n"
                + f"spin densities number: {len(self.spin_densities)}\n"
                + f"alpha FMO orbits number: {len(self.alpha_FMO_orbits)}\n"
                + f"alpha energy: {self.alpha_energy}\n"
                + f"beta FMO orbits number: {len(self.beta_FMO_orbits)}\n"
                + f"beta energy: {self.beta_energy}\n"
                + f"spin eigenvalue: {self.spin_eigenvalue}\n"
                + f"spin multiplicity: {self.spin_multiplicity}\n"
                # + f"nbo analysis number: {len(self.nbo_analysis)}\n"
                # + f"hessian number: {len(self.hessian)}\n"
            )

    def ts_vibration(self) -> List[BaseBlockParser]:
        """
        Generate a list of base block parsers for transition state vibration calculations.

        Returns:
            A list of base block parsers for transition state vibration calculations.
        """
        if not self.is_TS:
            raise RuntimeError("This is not a TS")

        block_parsers = []  # Initialize a list of base block parsers

        # Iterate over a list of ratios
        for ratio in (-2.5, -2.25, -2, -1.5, -1, 1, 1.5, 2, 2.25, 2.5):
            # Calculate extreme coordinates based on current ratio
            extreme_coords = (
                self.rdmol.GetConformer().GetPositions()
                - np.array(
                    self.dimensionless_imaginary_frequencies.__next__()[
                        "normal coordinates"
                    ]
                )
                * ratio
            )

            # Convert extreme coordinates to openbabel molecule object
            omol = xyz_block_to_omol(
                f"{len(self.atoms)}\n"
                + f"charge {self.charge} multiplicity {self.multiplicity}\n"
                + "\n".join(
                    [
                        f"{atom:10s}{x:10.5f}{y:10.5f}{z:10.5f}"
                        for atom, x, y, z in zip(self.atoms, *zip(*extreme_coords))
                    ]
                ),
                given_charge=self.charge,
            )

            try:
                # Rebuild parser using the openbabel molecule object
                block_parser = self.rebuild_parser(
                    Chem.MolFromMolBlock(omol.write("sdf"), removeHs=False),
                    rebuild_type="reindex",
                )
            except:
                continue

            # Check if the molecule satisfies crowding conditions and append it to the list
            if check_crowding(block_parser.rdmol):
                block_parsers.append(block_parser)

        return block_parsers

    def ts_vibration_to_SDF_file(self, file_path: str = None) -> str:
        """
        Write the TS vibration calculations to an SDF file.

        Parameters:
            file_path str:
                The file path. If not specified, will be generated in situ.
        Returns:
            The absolute path of the SDF file.
        """
        if file_path is None:
            file_path = self._file_path
        if os.path.isdir(file_path):
            raise IsADirectoryError(f"{file_path} is a directory.")
        file_path = os.path.splitext(file_path)[0] + ".sdf"
        block_parsers = self.ts_vibration()
        with open(file_path, "w") as f:
            f.write("$$$$\n".join([frame.to_SDF_block() for frame in block_parsers]))
        f.close()
        return file_path

    def possible_pre_post_ts(self) -> Tuple[Chem.rdchem.Mol, Chem.rdchem.Mol]:
        """
        This method returns the possible pre- and post-transition state molecules.

        Returns:
            A tuple containing the possible pre- and post-transition state molecules.
        """
        block_parsers = self.ts_vibration()
        block_parsers[0].rdmol.RemoveAllConformers()
        block_parsers[-1].rdmol.RemoveAllConformers()
        return block_parsers[0].rdmol, block_parsers[-1].rdmol
