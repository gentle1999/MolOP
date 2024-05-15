import os
from abc import ABC, abstractmethod
from typing import Dict, List, Literal, Tuple, TypeVar, Union

import numpy as np
import numpy.typing as npt
from openbabel import pybel
from pint.facets.plain import PlainQuantity
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds
from rdkit.Chem.MolStandardize import rdMolStandardize

from molop.logger.logger import logger
from molop.structure.geometry import get_geometry_info
from molop.structure.structure import (
    get_bond_pairs,
    get_formal_charges,
    get_formal_num_radicals,
    bond_list,
)
from molop.utils.types import arrayNx3, RdMol
from molop.structure.structure_recovery import xyz_block_to_omol, omol_to_rdmol
from molop.unit import atom_ureg

pt = Chem.GetPeriodicTable()


class MolBlock(ABC):
    """
    Abstract class for a molecule block. The information is enough to recover the original molecule structure.

    Attributes:
        _atoms Union[List[str], List[int]]:
            **Need to be given**: Save the atoms. Consistent with the atomic order in the original file. Supports both element numbers and element symbols.
        _coords PlainQuantity:
            **Need to be given**: Save the coordinates described by [pint](https://pint.readthedocs.io/en/stable/). Consistent with the atomic order in the original file. Defualt unit is angstrom.
        _charge int:
            **Need to be given**: Save the total charge. Allowed charge range is -3 to 3.
        _multiplicity int:
            **Need to be given**: Save the multiplicity. Allowed multiplicity range is 1 to 5.
        _bonds List[Tuple[int, int, int]]:
            **Predicted**: Save the bonds. The elements are: atom1 index, atom2 index, bond order(follow rdkit bond order `Chem.rdchem.BondType`).
        _formal_charges List[int]:
            **Predicted**: Save the formal charges. Consistent with the atomic order.
        _formal_num_radicals List[int]:
            **Predicted**: Save the formal spins. Consistent with the atomic order.
        _omol str:
            **Predicted**: Save the sdf from openbabel.
        _rdmol Chem.rdchem.Mol:
            **Predicted**: Save the rdkit molecule object.
    """

    _atoms: Union[List[str], List[int]]
    _coords: PlainQuantity
    _standard_coords: PlainQuantity
    _charge: int
    _multiplicity: int
    _bonds: List[Tuple[int, int, int]]
    _formal_charges: List[int]
    _formal_num_radicals: List[int]
    _omol: str
    _rdmol: Chem.rdchem.Mol

    def __init__(self):
        self._atoms: Union[List[str], List[int]] = []
        self._coords: PlainQuantity = np.array([[]]) * atom_ureg.angstrom
        self._standard_coords: PlainQuantity = np.array([[]]) * atom_ureg.angstrom
        self._charge: int = 0
        self._multiplicity: int = 1
        self._bonds: List[Tuple[int, int, int]] = None
        self._formal_charges: List[int] = None
        self._formal_num_radicals: List[int] = None
        self._omol: pybel.Molecule = None
        self._rdmol: RdMol = None

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

    @property
    def coords(self) -> PlainQuantity:
        """
        Get the coordinates.

        Returns:
            A matrix of coordinates with unit. Default unit is angstrom.
        """
        if self._coords.m.size != 0:
            return self._coords
        else:
            return self._standard_coords

    @property
    def dimensionless_coords(self) -> arrayNx3[np.float32]:
        """
        Get the dimensionless coordinates.

        Units will be tranformed:
            - `angstrom`

        Returns:
            A `np.ndarray` of dimensionless coordinates
        """
        return self.coords.to("angstrom").m.astype(np.float32)

    @property
    def standard_coords(self) -> PlainQuantity:
        """
        Get the standard coordinates.
        """
        return self._standard_coords

    @property
    def dimensionless_standard_coords(self) -> np.ndarray:
        """
        Get the standard coordinates.

        Units will be tranformed:
            - `angstrom`

        Returns:
            A `np.ndarray` of dimensionless coordinates
        """
        return self.standard_coords.to("angstrom").m.astype(np.float32)

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

    def __get_topology(self):
        if self.rdmol is None:
            return
        self._bonds = get_bond_pairs(self.rdmol)
        self._formal_charges = get_formal_charges(self.rdmol)
        self._formal_num_radicals = get_formal_num_radicals(self.rdmol)

    @property
    def bonds(self) -> List[Tuple[int, int, int]]:
        """
        Get the bonds.

        Returns:
            A list of bonds. Each bond is a tuple of three integers. The first two integers are the atom indices of the bond. The third integer is the bond order(follow rdkit bond order `Chem.rdchem.BondType`).
        """
        if self._bonds is None:
            self.__get_topology()
        return self._bonds

    @property
    def formal_charges(self) -> List[int]:
        """
        Get the formal charges.

        Returns:
            A list of formal charges. Consistent with the atomic order.
        """
        if self._formal_charges is None:
            self.__get_topology()
        return self._formal_charges

    @property
    def formal_num_radicals(self) -> List[int]:
        """
        Get the formal spins.

        Returns:
            A list of formal spins. Consistent with the atomic order.
        """
        if self._formal_num_radicals is None:
            self.__get_topology()
        return self._formal_num_radicals

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
                        self.to_XYZ_block(), self._charge, int(self._multiplicity) - 1
                    ).write("cml")
                except Exception as e:
                    raise RuntimeError(f"{self._file_path}: {e}")
            else:
                assert (
                    self._formal_charges and self._formal_num_radicals
                ), "If bonds given, formal charges and spins must be provided."
                self._omol = self.to_CML_block()
        return pybel.readstring("cml", self._omol)

    @property
    def rdmol(self) -> Chem.rdchem.Mol:
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
                    if len(self.atoms) <= 1:
                        raise ValueError("Single atom is not allowed in this function.")
                    raw_mol = Chem.MolFromXYZBlock(self.to_XYZ_block())
                    conn_mol = Chem.Mol(raw_mol)
                    rdDetermineBonds.DetermineBonds(conn_mol, charge=self._charge)
                    self._rdmol = conn_mol
                    self.__get_topology()
                except:
                    logger.debug(
                        f"{self._file_path}: rdkit determinebonds failed. Use MolOP structure recovery instead."
                    )
                    # If failed, use MolOP implementation
                    try:
                        omol = self.omol
                    except Exception as e:
                        logger.error(
                            f"{self._file_path}: MolOP structure recovery failed. {e}"
                        )
                        return None

                    try:
                        self._rdmol = omol_to_rdmol(
                            omol, self._charge, self._multiplicity - 1
                        )
                    except Exception as e:
                        logger.error(
                            f"{self._file_path}: RDKit structure recovery failed. {e}"
                        )
                    finally:
                        return self._rdmol
            else:
                assert (
                    self._formal_charges and self._formal_num_radicals
                ), "If bonds given, formal charges and spins must be provided."
                rwmol = Chem.RWMol(Chem.MolFromXYZBlock(self.to_XYZ_block()))
                for bond in self._bonds:
                    rwmol.AddBond(bond[0], bond[1], bond_list[bond[2]])
                for atom, charge, spin in zip(
                    rwmol.GetAtoms(), self._formal_charges, self._formal_num_radicals
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
                    f"{atom:10s}{x:10.5f}{y:10.5f}{z:10.5f}"
                    for atom, (x, y, z) in zip(self.atoms, self.dimensionless_coords)
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

    @abstractmethod
    def __str__(self) -> str:
        raise NotImplementedError

    def __hash__(self) -> int:
        return hash(str(self))

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({str(self)})"

    def __len__(self) -> int:
        return len(self.atoms)

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
