"""
Author: TMJ
Date: 2024-01-11 21:02:36
LastEditors: TMJ
LastEditTime: 2024-02-09 20:48:19
Description: 请填写简介
"""

import os
import re
from typing import Any, Dict, List, Literal, Tuple, Union

import numpy as np
import pandas as pd
from openbabel import pybel
from pint.facets.plain import PlainQuantity
from rdkit import Chem

from molop.config import molopconfig
from molop.io.bases.MolBlock import MolBlock
from molop.logger.logger import logger
from molop.structure.geometry import standard_orient
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
from molop.utils.g16patterns import link0_parser
from molop.utils.types import arrayN, arrayNx3


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
    _prev_block: "BaseBlockParser"

    def __init__(self, block: str) -> None:
        super().__init__()
        self._file_path = None
        self._frameID: int = 0
        self._block = block
        self._next_block: "BaseBlockParser" = None
        self._prev_block: "BaseBlockParser" = None

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

    @property
    def file_name(self):
        return os.path.basename(self._file_path)

    @property
    def file_path(self):
        return self._file_path

    @property
    def pure_filename(self) -> str:
        return os.path.splitext(self.file_name)[0]

    @property
    def file_format(self) -> str:
        """
        Get the file format of the object.
        Returns:
            str: The file format of the object.
        """
        return os.path.splitext(self.file_path)[-1]

    def _check_path(self, file_path: str = None, format: str = ".xyz"):
        if file_path is None:
            return os.path.splitext(self._file_path)[0] + format
        if os.path.isdir(file_path):
            return os.path.join(
                file_path,
                self.pure_filename + format,
            )
        return file_path

    def to_XYZ_file(self, file_path: str = None) -> str:
        """
        Write the XYZ file.

        Parameters:
            file_path str:
                The file path. If not specified, will be generated in situ.
        Returns:
            The absolute path of the XYZ file.
        """
        _file_path = self._check_path(file_path, ".xyz")
        with open(_file_path, "w") as f:
            f.write(self.to_XYZ_block())
        f.close()
        return os.path.abspath(_file_path)

    def to_SDF_file(self, file_path: str = None) -> str:
        """
        Write the SDF file.

        Parameters:
            file_path str:
                The file path. If not specified, will be generated in situ.
        Returns:
            The absolute path of the SDF file.
        """
        _file_path = self._check_path(file_path, ".sdf")
        with open(_file_path, "w") as f:
            f.write(self.to_SDF_block())
        f.close()
        return os.path.abspath(_file_path)

    def to_GJF_block(
        self,
        file_path: str = None,
        charge: int = None,
        multiplicity: int = None,
        template: str = None,
        prefix: str = "#p opt b3lyp def2svp freq EmpiricalDispersion=GD3BJ NoSymm\n",
        suffix="",
        chk: bool = True,
        oldchk: bool = False,
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
            chk bool:
                If true, add the chk keyword to the link0 section. Will use the file name as the chk file name.
            oldchk bool:
                If true, add the oldchk keyword to the link0 section. Will use the file name as the chk file name.
        Returns:
            A modified GJF block.
        """
        _file_path = self._check_path(file_path, ".gjf")
        if template is not None:
            if not os.path.isfile(template):
                raise FileNotFoundError(f"{template} is not found.")
            with open(template, "r") as f:
                lines = f.readlines()
            f.close()
            for i, line in enumerate(lines):
                if re.match(r"^\s*[\+\-\d]+\s+\d+$", line):
                    prefix = "".join(lines[: i - 2])
                if re.search(
                    r"[a-zA-z0-9]+\s+([\s\-]\d+\.\d+)\s+([\s\-]\d+\.\d+)\s+([\s\-]\d+\.\d+)",
                    line,
                ):
                    suffix = lines[i + 1 :]
        suffix = "".join(suffix)
        link, route = prefix.split("#")
        link = link.replace("\n", " ")
        route = "#" + route.replace("\n", " ")

        _link0 = link0_parser(link)

        if chk:
            _link0[r"%chk"] = f"{os.path.splitext(os.path.basename(_file_path))[0]}.chk"
        if oldchk:
            _link0[r"%oldchk"] = (
                f"{os.path.splitext(os.path.basename(_file_path))[0]}.chk"
            )
        return (
            "\n".join([f"{key}={val}" for key, val in _link0.items()])
            + "\n"
            + route
            + "\n\n"
            + f" Title: {self.pure_filename}\n\n"
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
        prefix: str = "#p opt b3lyp def2svp freq EmpiricalDispersion=GD3BJ NoSymm\n",
        suffix="",
        chk: bool = True,
        oldchk: bool = False,
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
            chk bool:
                If true, add the chk keyword to the link0 section. Will use the file name as the chk file name.
            oldchk bool:
                If true, add the oldchk keyword to the link0 section. Will use the file name as the chk file name.
        Returns:
            The path to the GJF file.
        """
        _file_path = self._check_path(file_path, ".gjf")
        with open(_file_path, "w") as f:
            f.write(
                self.to_GJF_block(
                    file_path=file_path,
                    charge=charge,
                    multiplicity=multiplicity,
                    template=template,
                    prefix=prefix,
                    suffix=suffix,
                    chk=chk,
                    oldchk=oldchk,
                )
            )
        f.close()
        return os.path.abspath(_file_path)

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
        _file_path = self._check_path(file_path, ".cdxml")
        if not keep3D:
            temp_rdmol = Chem.RWMol(self.rdmol)
            temp_rdmol.RemoveAllConformers()
            pybel.readstring("sdf", Chem.MolToMolBlock(temp_rdmol)).write(
                "cdxml", _file_path, overwrite=True
            )
        else:
            self.omol.write("cdxml", _file_path, overwrite=True)
        return os.path.abspath(_file_path)

    def next(self) -> "BaseBlockParser":
        """
        Get the next block.

        Returns:
            The next block.
        """
        return self._next_block

    def to_dict(self) -> Dict[str, Any]:
        """
        Get the necessary information for molecule recovery as a dict.

        Contents:
            - smiles: str
            - atom_number: int
            - total_charge: int
            - total_multiplicity: int
            - atoms: List[str]
            - coords: List[Tuple[float, float, float]]
            - bonds: List[Tuple[int, int, int]]
            - formal_charges: List[int]
            - formal_spins: List[int]

        Returns:
            The information of the block as a dict.
        """
        return {
            "smiles": self.to_standard_SMILES(),
            "atom_number": len(self),
            "total_charge": self.charge,
            "total_multiplicity": self.multiplicity,
            "atoms": [Chem.Atom(atom).GetAtomicNum() for atom in self.atoms],
            "coords": self.dimensionless_coords.tolist(),
            "bonds": self.bonds,
            "formal_charges": self.formal_charges,
            "formal_spins": self.formal_spins,
        }

    @staticmethod
    def _flatten(sequence):
        if sequence is None:
            return np.array([])
        else:
            return np.array(sequence).flatten()

    def rebuild_parser(
        self,
        new_mol: Union[Chem.rdchem.Mol, Chem.rdchem.RWMol],
        rebuild_type=Literal["mod", "reindex", "transform"],
    ) -> "BaseBlockParser":
        """
        Create a new parser with the given molecule.

        Parameters:
            new_mol Union[Chem.rdchem.Mol, Chem.rdchem.RWMol]:
                The new molecule.
            rebuild_type Literal["mod", "reindex", "transform"]:
                The tag of the new parser.
        Returns:
            The new parser.
        """
        new_parser = BaseBlockParser(block=Chem.MolToXYZBlock(new_mol))
        new_parser._atoms = [atom.GetSymbol() for atom in new_mol.GetAtoms()]
        new_parser._coords = new_mol.GetConformer().GetPositions() * atom_ureg.angstrom
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

    def standard_orient(
        self,
        anchor_list: List[int],
    ):
        """
        Depending on the input `idx_list`, `translate_anchor`, `rotate_anchor_to_X`, and `rotate_anchor_to_XY` are executed in order to obtain the normalized oriented molecule.

        Sub-functions:
            - `translate_anchor`: Translate the entire molecule so that the specified atom reaches the origin.
            - `rotate_anchor_to_X`: Rotate the specified second atom along the axis passing through the origin so that it reaches the positive half-axis of the X-axis.
            - `rotate_anchor_to_XY`: Rotate along the axis passing through the origin so that the specified third atom reaches quadrant 1 or 2 of the XY plane.

        Parameters:
            anchor_list List[int]:
                A list of indices of the atoms to be translated to origin, rotated to X axis, and rotated again to XY face:

                - If length is 1, execute `translate_anchor`
                - If length is 2, execute `translate_anchor` and `rotate_anchor_to_X`
                - If length is 3, execute `translate_anchor`, `rotate_anchor_to_X` and `rotate_anchor_to_XY`
                - If the length of the input `idx_list` is greater than 3, subsequent atomic numbers are not considered.
        """
        mol = Chem.RWMol(self.rdmol)
        standard_orient(mol, anchor_list)
        return self.rebuild_parser(
            mol,
            rebuild_type="transform",
        )

    def is_error(self) -> bool:
        """
        Abstrcact method to check if the current frame is an error. The details are implemented in the derived classes.
        """
        return False

    def is_TS(self) -> bool:
        """
        Check if the molecule is a TS. Can not check if the molecule is a TS without frequency information. Thus this function returns False.

        Returns:
            False
        """
        return False

    def to_summary_series(self) -> pd.Series:
        """
        Generate a summary series for the current frame.
        """
        return pd.Series(
            {
                "parser": self.__class__.__name__,
                "file_path": self.file_path,
                "file_name": self.file_name,
                "file_format": self.file_format,
                "frame_index": self._frameID,
                "charge": self.charge,
                "multiplicity": self.multiplicity,
                "SMILES": self.to_standard_SMILES(),
            }
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
        _scf_energy PlainQuantity:
            The total energy of the molecule.
        _mulliken_charges List[float]:
            The mulliken_charges of the molecule.
        _spin_densities List[float]:
            The spin densities of the molecule.
        _spin_multiplicity float:
            The spin multiplicity of the molecule.
        _spin_eigenvalue float:
            The spin eigenvalue of the molecule. spin egigenvalue = sqrt(spin multiplicity * (spin multiplicity + 1)).
        _gradients PlainQuantity:
            The gradients (forces) of the molecule.
        _frequencies List[Dict[str, Union[bool, PlainQuantity]]]:
            The frequencies of the molecule. The frequency dict has the following keys:

                - is imaginary: bool
                - freq: PlainQuantity
                - reduced masses: PlainQuantity
                - IR intensities: PlainQuantity
                - force constants: PlainQuantity
                - normal coordinates: PlainQuantity
        _alpha_FMO_orbits PlainQuantity:
            The alpha FMO orbitals of the molecule.
        _beta_FMO_orbits PlainQuantity:
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

                - zero-point sum: PlainQuantity
                - E sum: PlainQuantity
                - H sum: PlainQuantity
                - G sum: PlainQuantity
                - zero-point correction: PlainQuantity
                - TCE: PlainQuantity
                - TCH: PlainQuantity
                - TCG: PlainQuantity
        _state Dict[str, bool]:
            The state of the molecule.
        _only_extract_structure bool:
            A flag to indicate whether to only extract the structure.
    """

    def __init__(self, block: str, only_extract_structure=False) -> None:
        super().__init__(block)
        self._only_extract_structure = only_extract_structure

        # qm software basic information
        self.qm_software: str = ""
        self.version: str = ""
        self.parameter_comment: str = ""
        self.basis: str = ""
        self.functional: str = ""
        self.temperature: float = 298.15

        # Solvent properties
        self.solvent_model: str = ""
        self.solvent: str = ""
        self.solvent_eps: float = None
        self.solvent_eps_inf: float = None

        # Molecule properties
        self._total_energy: PlainQuantity = None
        self.scf_energy: PlainQuantity = None
        self.mp2_energy: PlainQuantity = None
        self.mp3_energy: PlainQuantity = None
        self.mp4_energy: PlainQuantity = None
        self.ccsd_energy: PlainQuantity = None
        self.spin_multiplicity: float = None
        self.spin_eigenvalue: float = None
        self.alpha_FMO_orbits: PlainQuantity = (
            np.array([]) * atom_ureg.hartree / atom_ureg.particle
        )
        self.beta_FMO_orbits: PlainQuantity = (
            np.array([]) * atom_ureg.hartree / atom_ureg.particle
        )
        self.alpha_energy: Dict[Literal["gap", "homo", "lumo"], PlainQuantity] = {
            "gap": None,
            "homo": None,
            "lumo": None,
        }
        self.beta_energy: Dict[Literal["gap", "homo", "lumo"], PlainQuantity] = {
            "gap": None,
            "homo": None,
            "lumo": None,
        }
        self.mulliken_charges: List[float] = []
        self.mulliken_spin_densities: List[float] = []
        self.apt_charges: List[float] = []
        self.lowdin_charges: List[float] = []
        self.dipole: PlainQuantity = np.array([]) * atom_ureg.debye
        self.quadrupole: PlainQuantity = (
            np.array([]) * atom_ureg.debye * atom_ureg.angstrom
        )
        self.octapole: PlainQuantity = (
            np.array([]) * atom_ureg.debye * atom_ureg.angstrom * atom_ureg.angstrom
        )
        self.hexadecapole: PlainQuantity = (
            np.array([])
            * atom_ureg.debye
            * atom_ureg.angstrom
            * atom_ureg.angstrom
            * atom_ureg.angstrom
        )
        self.hirshfeld_charges: List[float] = []
        self.hirshfeld_spins: List[float] = []
        self.hirshfeld_q_cm5: List[float] = []
        self.sum_energy = {
            "zero-point sum": None,
            "E sum": None,
            "H sum": None,
            "G sum": None,
            "zero-point correction": None,
            "TCE": None,
            "TCH": None,
            "TCG": None,
        }
        self.frequencies: List[
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
                ],
            ],
        ] = []
        self._dimensionless_frequencies = []
        self.gradients: PlainQuantity = (
            np.array([[]]) * atom_ureg.hartree / atom_ureg.bohr
        )
        self.wiberg_bond_order: np.ndarray[np.float32] = np.array([[]])
        self.mo_bond_order: np.ndarray[np.float32] = np.array([[]])
        self.atom_atom_overlap_bond_order: np.ndarray[np.float32] = np.array([[]])
        self.nbo_bond_order: List[Tuple[int, int, int, PlainQuantity]] = []
        self.nbo_charges: List[float] = []
        self.status = {}

    @property
    def total_energy(self) -> PlainQuantity:
        """
        Get the total energy of the molecule.
        """
        if self._total_energy is None:
            keys = list(self.energy.keys())
            if len(keys) > 0:
                return self.energy[keys[0]]
            else:
                return None
        return self._total_energy

    @property
    def energy(self):
        return {
            energy_type: getattr(self, f"{energy_type}_energy")
            for energy_type in ("ccsd", "mp4", "mp3", "mp2", "scf")
            if getattr(self, f"{energy_type}_energy") is not None
        }

    @property
    def dimensionless_total_energy(self) -> float:
        """
        Get the dimensionless total energy of the molecule.

        Unit will be transformed:
            - `hartree/particle`

        Returns:
            The dimensionless total energy of the molecule.
        """
        if self.total_energy is None:
            return None
        return self.total_energy.to("hartree/particle").m

    @property
    def dimensionless_scf_energy(self) -> float:
        """
        Get the dimensionless total energy of the molecule.

        Unit will be transformed:
            - `hartree/particle`

        Returns:
            The dimensionless total energy of the molecule.
        """
        if self.scf_energy is None:
            return None
        return self.scf_energy.to("hartree/particle").m

    @property
    def dimensionless_gradients(
        self,
    ) -> arrayNx3[np.float32]:
        """
        Get the dimensionless gradients of the molecule.

        Unit will be transformed:
            - `hartree/bohr`

        Returns:
            The dimensionless gradients of the molecule.
        """
        return self.gradients.to("hartree/bohr").m

    # @property
    # def hessian(self) -> List[float]:
    # return self._hessian

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
            Union[bool, PlainQuantity],
        ]
    ]:
        """
        Get the imaginary frequencies of the molecule.

        Returns:
            The imaginary frequencies of the molecule.
        """
        return [freq for freq in self.frequencies if freq["is imaginary"]]

    @property
    def dimensionless_imaginary_frequencies(
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
            Union[bool, float, arrayNx3[np.float32]],
        ]
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
        return [freq for freq in self.dimensionless_frequencies if freq["is imaginary"]]

    @property
    def dimensionless_frequencies(
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
            Union[bool, float, arrayNx3[np.float32]],
        ]
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
        if len(self._dimensionless_frequencies) == 0:
            self._dimensionless_frequencies = [
                {
                    "freq": freq["freq"].to("cm^-1").m,
                    "is imaginary": freq["is imaginary"],
                    "reduced masses": freq["reduced masses"].to("amu").m,
                    "IR intensities": freq["IR intensities"].to("kmol/mol").m,
                    "force constants": freq["force constants"].to("mdyne/angstrom").m,
                    "normal coordinates": freq["normal coordinates"].to("angstrom").m,
                }
                for freq in self.frequencies
            ]
        return self._dimensionless_frequencies

    def first_freq(self, dimensionless=False):
        """
        Get the first frequency of the molecule.
        """
        if len(self.frequencies) < 1:
            return {
                "freq": None,
                "is imaginary": None,
                "reduced masses": None,
                "IR intensities": None,
                "force constants": None,
                "normal coordinates": None,
            }
        else:
            if dimensionless:
                return self.dimensionless_frequencies[0]
            else:
                return self.frequencies[0]

    def second_freq(self, dimensionless=False):
        """
        Get the second frequency of the molecule.
        """
        if len(self.frequencies) < 2:
            return {
                "freq": None,
                "is imaginary": None,
                "reduced masses": None,
                "IR intensities": None,
                "force constants": None,
                "normal coordinates": None,
            }
        else:
            if dimensionless:
                return self.dimensionless_frequencies[1]
            else:
                return self.frequencies[1]

    @property
    def dimensionless_alpha_FMO_orbits(self) -> arrayN[np.float32]:
        """
        Get the dimensionless alpha FMO orbitals of the molecule.

        Unit will be transformed:
            - `hartree/particle`

        Returns:
            The dimensionless alpha FMO orbitals of the molecule.
        """
        return self.alpha_FMO_orbits.to("hartree/particle").m.astype(np.float32)

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
    def dimensionless_beta_FMO_orbits(self) -> arrayN[np.float32]:
        """
        Get the dimensionless beta FMO orbitals of the molecule.

        Unit will be transformed:
            - `hartree/particle`

        Returns:
            The dimensionless beta FMO orbitals of the molecule.
        """
        return self.beta_FMO_orbits.to("hartree/particle").m

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
    def dimensionless_sum_energy(
        self,
    ) -> Dict[
        Literal[
            "zero-point sum",
            "E sum",
            "H sum",
            "G sum",
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
    def dimensionless_dipole(self) -> np.ndarray:
        """
        Get the dimensionless dipole of the molecule.

        Unit will be transformed:
            - `debye`
        """
        return self.dipole.to("debye").m

    @property
    def dimensionless_quadrupole(self) -> np.ndarray:
        """
        Get the dimensionless quadrupole of the molecule.
        """
        return self.quadrupole.to("debye*angstrom").m

    @property
    def dimensionless_octapole(self) -> np.ndarray:
        """
        Get the dimensionless octapole of the molecule.
        """
        return self.octapole.to("debye*angstrom*angstrom").m

    @property
    def dimensionless_hexadecapole(self) -> np.ndarray:
        """
        Get the dimensionless hexadecapole of the molecule.
        """
        return self.hexadecapole.to("debye*angstrom*angstrom*angstrom").m

    @property
    def dimensionless_nbo_bond_order(self) -> List[Tuple[int, int, int, float]]:
        """
        Get the NBO bond order and dimensionless energy of the molecule.

        Unit will be transformed:
            - `hartree/particle`
        """
        return [
            (a, b, c, d.to("hartree/particle").m) for a, b, c, d in self.nbo_bond_order
        ]

    def to_dict(self):
        return {
            **super().to_dict(),
            **{
                "qm_software": self.qm_software,
                "qm_software_version": self.version,
                "functional": self.functional,
                "basis": self.basis,
                "solvent_model": self.solvent_model,
                "solvent": self.solvent,
                "temperature": self.temperature,
                "mulliken_charge": self.mulliken_charges,
                "spin_densities": self.mulliken_spin_densities,
                "gradients": self.dimensionless_gradients.tolist(),
                "single_point_energy": self.dimensionless_total_energy,
                "zero_point_correction": self.dimensionless_sum_energy[
                    "zero-point correction"
                ],
                "energy_correction": self.dimensionless_sum_energy["TCE"],
                "enthalpy_correction": self.dimensionless_sum_energy["TCH"],
                "gibbs_free_energy_correction": self.dimensionless_sum_energy["TCG"],
                "zero_point_sum": self.dimensionless_sum_energy["zero-point sum"],
                "thermal_energy_sum": self.dimensionless_sum_energy["E sum"],
                "thermal_enthalpy_sum": self.dimensionless_sum_energy["H sum"],
                "thermal_free_energy_sum": self.dimensionless_sum_energy["G sum"],
                "alpha_orbital_energies": self._flatten(
                    self.dimensionless_alpha_FMO_orbits
                ).tolist(),
                "beta_orbital_energies": self._flatten(
                    self.dimensionless_beta_FMO_orbits
                ).tolist(),
                "alpha_homo": self.dimensionless_alpha_energy["homo"],
                "alpha_lumo": self.dimensionless_alpha_energy["homo"],
                "alpha_gap": self.dimensionless_alpha_energy["gap"],
                "beta_homo": self.dimensionless_beta_energy["homo"],
                "beta_lumo": self.dimensionless_beta_energy["lumo"],
                "beta_gap": self.dimensionless_beta_energy["gap"],
                "first_frequency": self.first_freq(dimensionless=True)["freq"],
                "first_vibration_mode": self._flatten(
                    self.first_freq(dimensionless=True)["normal coordinates"]
                ).tolist(),
                "second_frequency": self.second_freq(dimensionless=True)["freq"],
                "second_vibration_mode": self._flatten(
                    self.second_freq(dimensionless=True)["normal coordinates"]
                ).tolist(),
                "freqs": self.dimensionless_frequencies,
                "is_TS": self.is_TS(),
                "spin_eginvalue": self.spin_eigenvalue,
                "spin_multiplicity": self.spin_multiplicity,
                "dipole": self._flatten(self.dimensionless_dipole).tolist(),
                "quadrupole": self._flatten(self.dimensionless_quadrupole).tolist(),
                "octapole": self._flatten(self.dimensionless_octapole).tolist(),
                "hexadecapole": self._flatten(self.dimensionless_hexadecapole).tolist(),
                "nbo_bond_order": self.dimensionless_nbo_bond_order,
                "wiberg_bond_order": self.wiberg_bond_order.tolist(),
                "mo_bond_order": self.mo_bond_order.tolist(),
                "atom_atom_overlap_bond_order": self.atom_atom_overlap_bond_order.tolist(),
                "nbo_charges": self.nbo_charges,
                "lowdin_charges": self.lowdin_charges,
                "hirshfeld_charges": self.hirshfeld_charges,
            },
        }

    def ts_vibration(self) -> List[BaseBlockParser]:
        """
        Generate a list of base block parsers for transition state vibration calculations.

        Returns:
            A list of base block parsers for transition state vibration calculations.
        """
        if not self.is_TS():
            raise RuntimeError("This is not a TS")

        block_parsers = []  # Initialize a list of base block parsers

        # Iterate over a list of ratios
        for ratio in (-2.5, -2.25, -2, -1.5, -1, 1, 1.5, 2, 2.25, 2.5):
            # Calculate extreme coordinates based on current ratio
            extreme_coords = (
                self.rdmol.GetConformer().GetPositions()
                - self.dimensionless_imaginary_frequencies[0]["normal coordinates"]
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

    def __repr__(self) -> str:
        return (
            f"{self.__class__.__name__}({str(self)})\n"
            + f"functional: {self.functional}\n"
            + f"basis: {self.basis}\n"
            + f"solvent_model: {self.solvent_model}\n"
            + f"solvent: {self.solvent}"
        )

    def to_summary_series(self) -> pd.Series:
        """
        Generate a summary series for the current frame.
        """
        return pd.Series(
            {
                "parser": self.__class__.__name__,
                "file_path": self.file_path,
                "file_name": self.file_name,
                "file_format": self.file_format,
                "version": self.version,
                "frame_index": self._frameID,
                "charge": self.charge,
                "multiplicity": self.multiplicity,
                "SMILES": self.to_standard_SMILES(),
                "functional": self.functional,
                "basis": self.basis,
                "solvent_model": self.solvent_model,
                "solvent": self.solvent,
                "temperature": self.temperature,
                "status": str(self.status),
                "ZPE": self.dimensionless_sum_energy["zero-point correction"],
                "TCE": self.dimensionless_sum_energy["TCE"],
                "TCH": self.dimensionless_sum_energy["TCH"],
                "TCG": self.dimensionless_sum_energy["TCG"],
                "ZPE-Gas": self.dimensionless_sum_energy["zero-point sum"],
                "E-Gas": self.dimensionless_sum_energy["E sum"],
                "H-Gas": self.dimensionless_sum_energy["H sum"],
                "G-Gas": self.dimensionless_sum_energy["G sum"],
                "sp": self.dimensionless_total_energy,
                "HOMO": self.dimensionless_alpha_energy["homo"],
                "LUMO": self.dimensionless_alpha_energy["lumo"],
                "GAP": self.dimensionless_alpha_energy["gap"],
                "first freq": self.first_freq(dimensionless=True)["freq"],
                "first freq tag": self.first_freq(dimensionless=True)["is imaginary"],
                "second freq": self.second_freq(dimensionless=True)["freq"],
                "second freq tag": self.second_freq(dimensionless=True)["is imaginary"],
                "S**2": self.spin_eigenvalue,
                "S": self.spin_multiplicity,
            }
        )
