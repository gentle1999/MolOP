"""
Author: TMJ
Date: 2024-01-07 13:47:18
LastEditors: TMJ
LastEditTime: 2024-01-25 22:49:10
Description: 请填写简介
"""
import os
from abc import ABC, abstractmethod
from collections.abc import Iterable
import re
from typing import Any, Dict, List, Literal, Tuple, Union

from openbabel import pybel
from pint.facets.plain import PlainQuantity
from rdkit import Chem
from rdkit.Chem import AllChem, rdDetermineBonds
from rdkit.Chem.MolStandardize import rdMolStandardize
from molop.unit import atom_ureg

from molop.logger.logger import logger
from molop.structure.structure import (
    get_bond_pairs,
    get_formal_charges,
    get_formal_spins,
    get_resonance_structures,
    structure_score,
    attempt_replacement,
)
from molop.structure.structure_recovery import xyz_block_to_omol
from molop.structure.geometry import get_geometry_info


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
    _iter_resonance: bool

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
        self._iter_resonance = False

    @property
    def atoms(self) -> List[str]:
        """
        Get the atoms.
        """
        return [Chem.Atom(atom).GetSymbol() for atom in self._atoms]

    @property
    def total_electrons(self) -> int:
        return sum(Chem.Atom(atom).GetAtomicNum() for atom in self._atoms) + self.charge

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
                try:
                    self._omol = xyz_block_to_omol(
                        self.to_XYZ_block(),
                        self._charge,
                        (self._multiplicity - 1) // 2,
                    )
                except Exception as e:
                    raise RuntimeError(f"{self._file_path}: {e}")
            else:
                omol = pybel.readstring("sdf", self.to_SDF_block())
                self._omol = omol
        return self._omol

    @property
    def rdmol(self) -> Chem.rdchem.Mol:
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
                    # If failed, use our implementation
                    self._rdmol = Chem.MolFromMolBlock(
                        self.omol.write("sdf"), removeHs=False
                    )
                    if self._rdmol is None:
                        raise ValueError(
                            f"{self._file_path}: rdkit determinebonds failed. MolOP structure recovery failed."
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
            if self._iter_resonance:
                self._resonance()

        return self._rdmol

    def _resonance(self):
        resmols = get_resonance_structures(
            self._rdmol, flags=Chem.ResonanceFlags.ALLOW_INCOMPLETE_OCTETS
        )
        resmols.sort(key=structure_score)
        self._rdmol = resmols[0]

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
        smiles = Chem.MolToSmiles(self.rdmol)
        return smiles

    def to_standard_SMILES(self) -> str:
        return rdMolStandardize.StandardizeSmiles(self.to_SMILES())

    def to_InChI(self) -> str:
        return Chem.MolToInchi(self.rdmol)

    def to_XYZ_file(self, file_path: str = None):
        if file_path is None:
            file_path = self._file_path
        if os.path.isdir(file_path):
            raise IsADirectoryError(f"{file_path} is a directory.")
        file_path = os.path.splitext(file_path)[0] + ".xyz"
        with open(file_path, "w") as f:
            f.write(self.to_XYZ_block())
        f.close()
        return file_path

    def to_SDF_file(self, file_path: str = None):
        if file_path is None:
            file_path = self._file_path
        if os.path.isdir(file_path):
            raise IsADirectoryError(f"{file_path} is a directory.")
        file_path = os.path.splitext(file_path)[0] + ".sdf"
        with open(file_path, "w") as f:
            f.write(self.to_SDF_block())
        f.close()
        return file_path

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

    def calc_mordred_descs(self, **kwargs):
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
        prefix: str = f"# g16 gjf \n",
        suffix="",
    ) -> str:
        if template is not None:
            if not os.path.isfile(template):
                raise FileNotFoundError(f"{template} is not found.")
            with open(template, "r") as f:
                lines = f.readlines()
            f.close()
            for i, line in enumerate(lines):
                if re.match(r"^\s*[\+\-\d]+\s+\d+$", line):
                    prefix = "".join(lines[:i-2])
                if re.match(r"^\s*[A-Z][a-z]?(\s+\-?\d+(\.\d+)?){3}$", line):
                    suffix = "".join(lines[i+1:])
        else:
            prefix = prefix if prefix.endswith("\n") else prefix + "\n"
            prefix = prefix + "\n"
        return (
            prefix
            + f" Title\n\n"
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
        prefix: str = f"# g16 gjf \n",
        suffix="",
    ):
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
        return file_path

    def to_chemdraw(self, file_path: str = None, keep3D=True):
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
        return file_path

    def geometry_analysis(self, atom_idxs: Tuple[int], one_start=False):
        """
        Get the geometry infos among the atoms

        Parameters:
            atom_idxs Tuple[int]:
                A list of index of the atoms, starts from 0
            one_start bool:
                If true, consider atom index starts from 1, so let index value subtracts 1 for all the atoms

        Returns:
            A float value:
                If the length of atom_idxs is 2, the bond length with unit Angstrom between the two atoms will be returned.

                If the length of atom_idxs is 3, the angle with unit degree between  the three atoms will be returned.

                If the length of atom_idxs is 4, the dihedral angle with unit degree between the four atoms will be returned.
        """
        if one_start:
            atom_idxs = [atom_idx - 1 for atom_idx in atom_idxs]
        return get_geometry_info(self.rdmol, atom_idxs)


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
            "formal charges": self.formal_charges,
            "formal spins": self.formal_spins,
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

    def rebuild_parser(self, new_mol, rebuild_type=Literal["mod", "reindex"]):
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
    ):
        new_mol = attempt_replacement(
            self.rdmol,
            query_smi=query_smi,
            replacement_smi=replacement_smi,
            bind_idx=bind_idx,
            replace_all=replace_all,
            attempt_num=attempt_num,
        )
        return self.rebuild_parser(new_mol, rebuild_type="mod")

    def reset_atom_index(self, mapping_smarts: str):
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
        coords = [self.dimensionless_coords[1][idx] for idx in mapping] + [
            (x, y, z)
            for idx, (x, y, z) in enumerate(self.dimensionless_coords[1])
            if idx not in mapping
        ]
        atoms = [self.atoms[idx] for idx in mapping] + [
            atom for idx, atom in enumerate(self.atoms) if idx not in mapping
        ]
        omol = xyz_block_to_omol(
            f"{len(atoms)}\n"
            + f"charge {self.charge} multiplicity {self.multiplicity}\n"
            + "\n".join(
                [
                    f"{atom:10s}{x:10.5f}{y:10.5f}{z:10.5f}"
                    for atom, x, y, z in zip(atoms, *zip(*coords))
                ]
            ),
            given_charge=self.charge,
            given_spin=self.multiplicity,
        )
        return self.rebuild_parser(
            Chem.MolFromMolBlock(omol.write("sdf"), removeHs=False),
            rebuild_type="reindex",
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
        self._gradients: List[Tuple[PlainQuantity]] = []
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
        return self._energy

    @property
    def dimensionless_energy(self):
        """
        Unit will be transformed to hartree/particle
        """
        if self.energy is None:
            return None
        return self.energy.to("hartree/particle").m

    @property
    def partial_charges(self) -> List[float]:
        return self._partial_charges

    @property
    def spin_densities(self) -> List[float]:
        return self._spin_densities

    @property
    def spin_multiplicity(self) -> float:
        return self._spin_multiplicity

    @property
    def spin_eigenvalue(self) -> float:
        return self._spin_eigenvalue

    @property
    def gradients(self):
        return self._gradients

    @property
    def dimensionless_gradients(self):
        """
        Unit will be transformed to hartree/bohr
        """
        return (
            [gradient.to("hartree/bohr").m for gradient in atom]
            for atom in self.gradients
        )

    # @property
    # def hessian(self) -> List[float]:
    # return self._hessian

    @property
    def is_TS(self) -> bool:
        return len(list(self.imaginary_frequencies)) == 1

    @property
    def imaginary_frequencies(self):
        return (freq for freq in self._frequencies if freq["is imaginary"])

    @property
    def dimensionless_imaginary_frequencies(self):
        return (freq for freq in self.dimensionless_frequencies if freq["is imaginary"])

    @property
    def frequencies(self):
        return self._frequencies

    @property
    def dimensionless_frequencies(
        self,
    ):
        """
        freq: frequency in cm^-1
        Reduced masses: amu
        IR intensities: kmol/mol
        force constants: mdyne/angstrom
        normal coordinates: angstrom
        """
        return (
            {
                "freq": freq["freq"].to("cm^-1").m,
                "is imaginary": freq["is imaginary"],
                "reduced masses": freq["reduced masses"].to("amu").m,
                "IR intensities": freq["IR intensities"].to("kmol/mol").m,
                "force constants": freq["force constants"].to("mdyne/angstrom").m,
                "normal coordinates": [
                    [coord.to("angstrom").m for coord in atom]
                    for atom in freq["normal coordinates"]
                ],
            }
            for freq in self._frequencies
        )

    @property
    def alpha_FMO_orbits(self):
        return self._alpha_FMO_orbits

    @property
    def dimensionless_alpha_FMO_orbits(self):
        """
        Unit will be transformed to hartree/particle
        """
        return (orbit.to("hartree/particle").m for orbit in self.alpha_FMO_orbits)

    @property
    def alpha_energy(self):
        return self._alpha_energy

    @property
    def dimensionless_alpha_energy(self):
        """
        Unit will be transformed to hartree/particle
        """
        return {
            key: value.to("hartree/particle").m if value is not None else None
            for key, value in self.alpha_energy.items()
        }

    @property
    def beta_FMO_orbits(self):
        return self._beta_FMO_orbits

    @property
    def dimensionless_beta_FMO_orbits(self):
        """
        Unit will be transformed to hartree/particle
        """
        return (orbit.to("hartree/particle").m for orbit in self.beta_FMO_orbits)

    @property
    def beta_energy(self):
        return self._beta_energy

    @property
    def dimensionless_beta_energy(self):
        """
        Unit will be transformed to hartree/particle
        """
        return {
            key: value.to("hartree/particle").m if value is not None else None
            for key, value in self.beta_energy.items()
        }

    @property
    def sum_energy(self):
        return self._sum_energy

    @property
    def dimensionless_sum_energy(self):
        """
        Unit will be transformed to hartree/particle
        """
        return {
            key: value.to("hartree/particle").m if value is not None else None
            for key, value in self.sum_energy.items()
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

    # @property
    # def nbo_analysis(self):
    # return self._nbo_analysis

    def to_dict(self):
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
                # "nbo_analysis": self.nbo_analysis,
                # "hessian": self.hessian,
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
        )
        if not self._only_extract_structure:
            print(
                f"state: {self.state}\n"
                + f"energy: {self.energy}\n"
                + f"sum energy: {self.sum_energy}\n"
                + f"gradients number: {len(self.gradients)}\n"
                + f"frequencies number: {len(self.frequencies)}\n"
                + f"imaginary frequencies number: {len(list(self.imaginary_frequencies))}\n"
                + f"partial charges number: {len(self.partial_charges)}\n"
                + f"spin densities number: {len(self.spin_densities)}\n"
                + f"alpha FMO orbits number: {len(self.alpha_FMO_orbits)}\n"
                + f"alpha energy: {self.alpha_energy}\n"
                + f"beta FMO orbits number: {len(self.beta_FMO_orbits)}\n"
                + f"beta energy: {self.beta_energy}\n"
                # + f"nbo analysis number: {len(self.nbo_analysis)}\n"
                # + f"hessian number: {len(self.hessian)}\n"
            )
