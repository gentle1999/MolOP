"""
Author: TMJ
Date: 2024-01-07 13:47:18
LastEditors: TMJ
LastEditTime: 2024-01-09 16:01:09
Description: 请填写简介
"""
import os
from abc import ABC, abstractmethod
from collections.abc import Iterable
from copy import deepcopy
from typing import Any, Dict, List, Literal, Tuple, Union

from openbabel import pybel
from pint.facets.plain import PlainQuantity
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

from molop.structure.structure import (
    get_bond_pairs,
    get_formal_charges,
    get_formal_spins,
)
from molop.structure.structure_recovery import xyz_block_to_omol


class MolBlock(ABC):
    """
    Abstract class for a molecule block.
    """

    _atoms: Union[List[str], List[int]]
    _coords: List[Tuple[float]]
    _charge: int  # Only allow -2 ~ +2
    _multiplicity: int  # Only allow 1 ~ 5
    _bonds: List[Tuple[int, int, int]]
    _formal_charges: List[int]
    _formal_spins: List[int]
    _omol: pybel.Molecule
    _rdmol: Union[Chem.rdchem.RWMol, Chem.rdchem.Mol]

    def __init__(self):
        """
        Constructor.
        """
        self._atoms: Union[List[str], List[int]] = []
        self._coords: List[Tuple[float]] = []
        self._charge: int = 0  # Only allow -2 ~ +2
        self._multiplicity: int = 1  # Only allow 1 ~ 5
        self._bonds: List[Tuple[int, int, int]] = None
        self._formal_charges: List[int] = None
        self._formal_spins: List[int] = None
        self._omol: pybel.Molecule = None
        self._rdmol: Union[Chem.rdchem.RWMol, Chem.rdchem.Mol] = None
        self._check_spin = False

    @property
    def atoms(self) -> List[str]:
        """
        Get the atoms.
        """
        return [Chem.Atom(atom).GetSymbol() for atom in self._atoms]

    @property
    def elements(self) -> List[str]:
        return list(set(self.atoms))

    @property
    def coords(self) -> List[Tuple[PlainQuantity]]:
        """
        Get the coordinates.
        """
        return self._coords

    @property
    def dimensionless_coords(self):
        return (
            str(self.coords[0][0].units),
            [[coord.m for coord in atom] for atom in (self.coords)],
        )

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
    def bonds(self) -> List[Tuple[int, int, int]]:
        """
        Get the bonds.
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
        """
        if self._formal_spins is None:
            self._bonds = get_bond_pairs(self.rdmol)
            self._formal_charges = get_formal_charges(self.rdmol)
            self._formal_spins = get_formal_spins(self.rdmol)
        return self._formal_spins

    @property
    def omol(self):
        if self._omol is None:
            if self._bonds is None:
                self._omol = xyz_block_to_omol(
                    self.to_XYZ_block(),
                    self._charge,
                    (self._multiplicity - 1) // 2,
                    self._check_spin,
                )
            else:
                omol = pybel.readstring("sdf", self.to_SDF_block())
                self._omol = omol
        return self._omol

    @property
    def rdmol(self) -> Chem.rdchem.Mol:
        if self._rdmol is None:
            if self._bonds is None:
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
                    f"{atom:10s}{x.m:10.5f}{y.m:10.5f}{z.m:10.5f}"
                    for atom, x, y, z in zip(self.atoms, *zip(*self.coords))
                ]
            )
        )

    def to_SDF_block(self) -> str:
        return Chem.MolToMolBlock(self.rdmol)

    def to_SMILES(self) -> str:
        return rdMolStandardize.StandardizeSmiles(Chem.MolToSmiles(self.rdmol))

    def to_InChI(self) -> str:
        return Chem.MolToInchi(self.rdmol)

    def to_XYZ_file(self, file_path: str):
        with open(file_path, "w") as f:
            f.write(self.to_XYZ_block())
        f.close()

    def to_SDF_file(self, file_path: str):
        with open(file_path, "w") as f:
            f.write(self.to_SDF_block())
        f.close()

    def calc_rdkit_descs(self, desc_names: List[str] = None) -> Dict[str, float]:
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
    ):
        """
        Require additional dependencies:
            - [dscribe](https://singroup.github.io/dscribe/latest/)
            - [ase](https://wiki.fysik.dtu.dk/ase/index.html)
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

    def calc_mordred_descs(self):
        from molop.descriptor.descriptor import calc_mordred_descs

        return calc_mordred_descs(self.rdmol)

    @abstractmethod
    def __str__(self) -> str:
        raise NotImplementedError

    def __hash__(self) -> int:
        return hash(str(self))

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({str(self)})"

    def __len__(self) -> int:
        return len(self.atoms)


class BaseBlockParser(MolBlock):
    """
    Base class for block parsers.
    Only need basic information.
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
        return self._block

    def __str__(self) -> str:
        return (
            f"{os.path.basename(self._file_path)}[{self._frameID}]\n"
            + self.to_XYZ_block()
        )

    def next(self):
        return self._next_block

    def to_dict(self):
        return {
            "type": self._block_type,
            "file_path": self._file_path,
            "SMILES": self.to_SMILES(),
            "frameID": self._frameID,
            "atoms": self.atoms,
            "coords": self.dimensionless_coords,
            "bonds": self.bonds,
            "total charge": self.charge,
            "total multiplicity": self.multiplicity,
        }

    def summary(self):
        print(
            f"type: {self._block_type}\n"
            + f"file path: {self._file_path}\n"
            + f"frameID: {self._frameID}\n"
            + f"SMILES: {self.to_SMILES()}\n"
            + f"atom number: {len(self)}\n"
            + f"total charge: {self.charge}\n"
        )


class QMBaseBlockParser(BaseBlockParser):
    """
    Base class for QM block parsers.
    Should contain all the QM information pre-defined.
    """

    _version: str
    _parameter_comment: str
    _energy: float
    _partial_charges: List[float]
    _spin_densities: List[float]
    _gradients: List[Tuple[PlainQuantity]]
    _hessian: List[List[float]]
    _frequencies: List[
        Dict[
            Literal[
                "is imaginary",
                "freq",
                "Reduced masses",
                "IR intensities",
                "force constants",
                "normal coordinates",
            ],
            Union[
                bool,
                PlainQuantity,
                List[Tuple[PlainQuantity]],
            ],
        ],
    ]
    _alpha_FMO_orbits: List[PlainQuantity]
    _beta_FMO_orbits: List[PlainQuantity]
    _alpha_energy: Dict[Literal["gap", "homo", "lumo"], PlainQuantity]
    _beta_energy: Dict[Literal["gap", "homo", "lumo"], PlainQuantity]
    _sum_energy: Dict[
        Literal[
            "zero-point",
            "thermal energy",
            "thermal enthalpy",
            "thermal gibbs free energy",
            "zero-point correction",
            "thermal energy correction",
            "thermal enthalpy correction",
            "thermal gibbs free energy correction",
        ],
        PlainQuantity,
    ]
    _nbo_analysis: List[Dict[str, Dict[str, PlainQuantity]]]

    _state: Dict[str, bool]
    _only_extract_structure: bool

    def __init__(self, block: str, only_extract_structure=False) -> None:
        super().__init__(block)
        self._check_spin = False
        self._only_extract_structure = only_extract_structure

        self._version: str = None
        self._parameter_comment: str = None
        self._energy: PlainQuantity = None
        self._partial_charges: List[float] = []
        self._spin_densities: List[float] = []
        self._gradients: List[Tuple[PlainQuantity]] = []
        self._hessian: List[List[float]] = []
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
        self._sum_energy = {}
        self._frequencies = []
        self._nbo_analysis = []
        self._state = {}

    @property
    def energy(self) -> PlainQuantity:
        return self._energy

    @property
    def dimensionless_energy(self):
        return (
            str(self.energy.units),
            self.energy.m,
        )

    @property
    def partial_charges(self) -> List[float]:
        return self._partial_charges

    @property
    def spin_densities(self) -> List[float]:
        return self._spin_densities

    @property
    def gradients(self):
        return self._gradients

    @property
    def dimensionless_gradients(self):
        return (
            str(self.gradients[0][0].units),
            [[gradient.m for gradient in atom] for atom in self.gradients],
        )

    @property
    def hessian(self) -> List[float]:
        return self._hessian

    @property
    def imaginary_frequencies(self):
        return [freq for freq in self._frequencies if freq["is imaginary"]]

    @property
    def dimensionless_imaginary_frequencies(self):
        return [freq for freq in self.dimensionless_frequencies if freq["is imaginary"]]

    @property
    def frequencies(self):
        return self._frequencies

    @property
    def dimensionless_frequencies(self):
        return [
            {
                "is imaginary": freq["is imaginary"],
                "freq": (str(freq["freq"].units), freq["freq"].m),
                "Reduced masses": (
                    str(freq["Reduced masses"].units),
                    freq["Reduced masses"].m,
                ),
                "IR intensities": (
                    str(freq["IR intensities"].units),
                    freq["IR intensities"].m,
                ),
                "force constants": (
                    str(freq["force constants"].units),
                    freq["force constants"].m,
                ),
                "normal coordinates": (
                    str(freq["normal coordinates"][0][0].units),
                    [
                        [coord.m for coord in atom]
                        for atom in freq["normal coordinates"]
                    ],
                ),
            }
            for freq in self._frequencies
        ]

    @property
    def alpha_FMO_orbits(self):
        return self._alpha_FMO_orbits

    @property
    def dimensionless_alpha_FMO_orbits(self):
        return (
            str(self.alpha_FMO_orbits[0].units),
            [orbit.m for orbit in self.alpha_FMO_orbits],
        )

    @property
    def alpha_energy(self):
        return self._alpha_energy

    @property
    def dimensionless_alpha_energy(self):
        return {
            key: (str(value.units), value.m) for key, value in self.alpha_energy.items()
        }

    @property
    def beta_FMO_orbits(self):
        return self._beta_FMO_orbits

    @property
    def dimensionless_beta_FMO_orbits(self):
        return (
            str(self.beta_FMO_orbits[0].units),
            [orbit.m for orbit in self.beta_FMO_orbits],
        )

    @property
    def beta_energy(self):
        return self._beta_energy

    @property
    def dimensionless_beta_energy(self):
        return {
            key: (str(value.units), value.m) for key, value in self.beta_energy.items()
        }

    @property
    def sum_energy(self):
        return self._sum_energy

    @property
    def dimensionless_sum_energy(self):
        return {
            key: (str(value.units), value.m) for key, value in self.sum_energy.items()
        }

    @property
    def state(self):
        return self._state

    @property
    def parameter_comment(self) -> str:
        return self._parameter_comment

    @property
    def version(self) -> str:
        return self._version

    @property
    def nbo_analysis(self):
        return self._nbo_analysis

    def to_dict(self):
        return {
            **super().to_dict(),
            **{
                "version": self.version,
                "parameter_comment": self.parameter_comment,
                "state": self.state,
                "energy": self.dimensionless_energy,
                "sum_energy": self.dimensionless_sum_energy,
                "gradients": self.dimensionless_gradients,
                "frequencies": self.dimensionless_frequencies,
                "imaginary_frequencies": self.dimensionless_imaginary_frequencies,
                "partial_charges": self.partial_charges,
                "spin_densities": self.spin_densities,
                "alpha_FMO_orbits": self.dimensionless_alpha_FMO_orbits,
                "alpha_energy": self.dimensionless_alpha_energy,
                "beta_FMO_orbits": self.dimensionless_beta_FMO_orbits,
                "beta_energy": self.dimensionless_beta_energy,
                "nbo_analysis": self.nbo_analysis,
                "hessian": self.hessian,
            },
        }
    
    def summary(self):
        print(
            f"type: {self._block_type}\n"
            + f"file path: {self._file_path}\n"
            + f"frameID: {self._frameID}\n"
            + f"SMILES: {self.to_SMILES()}\n"
            + f"atom number: {len(self)}\n"
            + f"total charge: {self.charge}\n"
            + f"version: {self.version}\n"
            + f"parameter comment: \n{self.parameter_comment}\n"
            + f"state: {self.state}\n"
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
            + f"nbo analysis number: {len(self.nbo_analysis)}\n"
            + f"hessian number: {len(self.hessian)}\n"
        )
