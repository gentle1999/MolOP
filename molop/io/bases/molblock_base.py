"""
Author: TMJ
Date: 2024-01-07 13:47:18
LastEditors: TMJ
LastEditTime: 2024-01-07 13:52:31
Description: 请填写简介
"""
import os
from abc import ABC, abstractmethod
from typing import Any, List, Tuple, Union, str

from rdkit import Chem
from openbabel import pybel


class MolBlock(ABC):
    """
    Abstract class for a molecule block.
    """

    _atoms: List[Union[str, int]]
    _bonds: List[Tuple[int, int, int]] = None
    _partial_charges: List[int] = None
    _partial_spins: List[int] = None
    _charge: int = 0 # Only allow -2 ~ +2 
    _multiplicity: int = 0 # Only allow 1 ~ 5 
    _coords: List[float]

    def __init__(self):
        """
        Constructor.
        """

    @property
    def atoms(self) -> List[str]:
        """
        Get the atoms.
        """
        return [Chem.Atom(atom).GetSymbol() for atom in self._atoms]

    @property
    def bonds(self) -> List[Tuple[int, int, int]]:
        """
        Get the bonds.
        """
        return self._bonds

    @property
    def charge(self) -> int:
        """
        Get the charge.
        """
        return self._charge

    @property
    def multiplicity(self) -> int:
        """
        Get the multiplicity.
        """
        return self._multiplicity

    @property
    def coords(self) -> List[float]:
        """
        Get the coordinates.
        """
        return self._coords

    def basic_check(self) -> bool:
        """
        Check the basic information.
        """
        assert self._charge in [-2, -1, 0, 1, 2], "The charge must be -2 ~ +2."
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
        return (
            f"{len(self.atoms)}\n"
            + f"charge {self.charge} multiplicity {self.multiplicity}\n"
            + "\n".join(
                [
                    f"{atom} {x:.6f} {y:.6f} {z:.6f}"
                    for atom, x, y, z in zip(self.atoms, *zip(*self.coords))
                ]
            )
        )
    def to_SDF_block(self) -> str:
        if self._bonds is None:
            omol = pybel.readstring(format='xyz', string=self.to_XYZ_block())
            rdmol = Chem.RWmol(Chem.MolFromMolBlock(omol.write('sdf')))
            return omol.write('sdf')
        
        


    @abstractmethod
    def to_SMILES(self) -> str:
        raise NotImplementedError

    @abstractmethod
    def to_InChI(self) -> str:
        raise NotImplementedError

    @abstractmethod
    def to_RdMol(self) -> Chem.rdchem.Mol:
        if self._bonds:
            rwmol = Chem.MolFromXYZBlock(self.to_XYZ_block())
            for bond in self._bonds:
                rwmol.AddBond(bond[0], bond[1], bond[2])
        else:
            omol = pybel.readstring(format='xyz', string=self.to_XYZ_block())
            rwmol = Chem.RWmol(Chem.MolFromMolBlock(omol.write('sdf')))
        

    def __hash__(self) -> int:
        return hash(str(self))

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({str(self)})"


class BaseBlockParser(MolBlock):
    """
    Base class for block parsers.
    Only need basic information.
    """

    _filename: str
    _frameID: int = 0

    def __init__(self, block: str) -> None:
        self._block = block

    @property
    def block(self) -> str:
        return self._block

    def __str__(self) -> str:
        return f"{self.__class__.__name__}({os.path.basename(self._filename)})[{self._frameID}]"


class QMBaseBlockParser(BaseBlockParser):
    """
    Base class for QM block parsers.
    Should contain all the QM information pre-defined.
    """

    _energy: float
    _partial_charges: List[float] = []
    _gradient: List[float] = []
    _hessian: List[List[float]] = []
    _normal_modes: List[float] = []
    _normal_mode_frequencies: List[float] = []
    _FMO_orbits: List[List[float]] = []
    _HOMO_energy: float
    _LUMO_energy: float

    def __init__(self, block: str) -> None:
        super().__init__(block)

    @property
    def energy(self) -> float:
        return self._energy

    @property
    def partial_charges(self) -> List[float]:
        return self._partial_charges

    @property
    def gradient(self) -> List[float]:
        return self._gradient

    @property
    def hessian(self) -> List[List[float]]:
        return self._hessian

    @property
    def normal_modes(self) -> List[float]:
        return self._normal_modes

    @property
    def normal_mode_frequencies(self) -> List[float]:
        return self._normal_mode_frequencies

    @property
    def FMO_orbits(self) -> List[List[float]]:
        return self._FMO_orbits

    @property
    def HOMO_energy(self) -> float:
        return self._HOMO_energy

    @property
    def LUMO_energy(self) -> float:
        return self._LUMO_energy
