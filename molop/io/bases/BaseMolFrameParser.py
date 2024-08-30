"""
Author: TMJ
Date: 2024-06-18 13:04:26
LastEditors: TMJ
LastEditTime: 2024-06-21 21:19:20
Description: 请填写简介
"""

import os
import re
from typing import List, Sequence, Tuple, TypeVar, Union

import numpy as np
import pandas as pd
from openbabel import pybel
from pint.facets.plain import PlainQuantity
from pydantic import Field, PrivateAttr, computed_field, model_validator
from rdkit import Chem
from typing_extensions import Self

from molop.io.bases.BaseMolFrame import BaseMolFrame
from molop.io.bases.DataClasses import (
    BondOrders,
    ChargeSpinPopulations,
    Energies,
    MolecularOrbitals,
    Polarizability,
    SinglePointProperties,
    ThermalEnergies,
    TotalSpin,
    Vibrations,
    Status,
    GeometryOptimizationStatus,
)
from molop.logger.logger import moloplogger
from molop.structure.geometry import standard_orient
from molop.structure.structure import (
    attempt_replacement,
    get_bond_pairs,
    get_formal_charges,
    get_formal_num_radicals,
    get_total_charge,
    get_total_multiplicity,
    reset_atom_index,
)
from molop.structure.structure_recovery import xyz_block_to_rdmol
from molop.unit import atom_ureg
from molop.utils.g16patterns import link0_parser
from molop.utils.types import RdMol


class BaseMolFrameParser(BaseMolFrame):
    frame_id: int = Field(default=0, description="Frame ID")
    file_path: str = Field(default="", repr=False, exclude=True)
    frame_content: str = Field(default="", repr=False, exclude=True)
    _frame_type: str = PrivateAttr(default="")
    _next_frame: "BaseMolFrameParser" = PrivateAttr(default=None)
    _prev_frame: "BaseMolFrameParser" = PrivateAttr(default=None)
    status: Status = Field(default=Status(), description="Status of the frame")
    geometry_optimization_status: GeometryOptimizationStatus = Field(
        default=GeometryOptimizationStatus(),
        description="Geometry optimization status",
    )

    def _parse(self):
        raise NotImplementedError

    @model_validator(mode="after")
    def __parse_frame__(self) -> Self:
        if len(self.atoms) > 0:
            return self
        self._parse()
        return self

    @computed_field
    @property
    def filename(self) -> str:
        """
        Get the filename of the frame.
        Returns:
            str: The filename of the frame.
        """
        return os.path.basename(self.file_path)

    @property
    def pure_filename(self) -> str:
        """
        Get the pure filename of the frame.
        Returns:
            str: The pure filename of the frame.
        """
        return os.path.splitext(self.filename)[0]

    @property
    def file_dir_path(self) -> str:
        """
        Get the file directory path of the frame.
        Returns:
            str: The file directory path of the frame.
        """
        return os.path.dirname(self.file_path)

    @property
    def file_format(self) -> str:
        """
        Get the file format of the object.
        Returns:
            str: The file format of the object.
        """
        return os.path.splitext(self.file_path)[-1]

    def _check_path(self, file_path: str = None, format: str = ".xyz") -> str:
        """
        Check if the file path is valid or modify it.

        Parameters:
            file_path (str):
                The file path. If not specified, will be generated in situ.
            format (str):
                The file format. If not specified, will be generated in situ.
        Returns:
            str: The absolute path of the file.
        """
        if file_path is None:
            return os.path.splitext(self.file_path)[0] + format
        if os.path.isdir(file_path):
            return os.path.join(
                file_path,
                self.pure_filename + format,
            )
        if not file_path.endswith(format):
            file_path = file_path + format
        if not os.path.isabs(file_path):
            file_path = os.path.abspath(file_path)
        return file_path

    def to_XYZ_file(self, file_path: str = None) -> str:
        """
        Write the XYZ file.

        Parameters:
            file_path (str):
                The file path. If not specified, will be generated in situ.
        Returns:
            str: The absolute path of the XYZ file.
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
            file_path (str):
                The file path. If not specified, will be generated in situ.
        Returns:
            str: The absolute path of the SDF file.
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
        prefix: str = "",
        suffix="",
        chk: bool = True,
        oldchk: bool = False,
    ) -> str:
        """
        Get the GJF frame.

        Parameters:
            file_path (str):
                The path to write the GJF file. If not specified, will be generated in situ.
            charge (int):
                The forced charge. If specified, will be used to overwrite the charge in the gjf file.
            multiplicity (int):
                The forced multiplicity. If specified, will be used to overwrite the multiplicity in the gjf file.
            template (str):
                path to read a gjf file as a template.
            prefix (str):
                prefix to add to the beginning of the gjf file, priority is higher than template.
            suffix (str):
                suffix to add to the end of the gjf file, priority is higher than template.
            chk (bool):
                If true, add the chk keyword to the link0 section. Will use the file name as the chk file name.
            oldchk (bool):
                If true, add the oldchk keyword to the link0 section. Will use the file name as the chk file name.
        Returns:
            str: A modified GJF frame.
        """
        if not isinstance(prefix, str) or not isinstance(suffix, str):
            raise TypeError("prefix and suffix must be strings.")
        _prefix, _suffix = "", ""
        _file_path = self._check_path(file_path, ".gjf")
        if template is not None:
            if not os.path.isfile(template):
                raise FileNotFoundError(f"{template} is not found.")
            with open(template, "r") as f:
                lines = f.readlines()
            f.close()
            for i, line in enumerate(lines):
                if re.match(r"^\s*[\+\-\d]+\s+\d+$", line):
                    _prefix = "".join(lines[: i - 2])
                if re.search(
                    r"[a-zA-z0-9]+\s+([\s\-]\d+\.\d+)\s+([\s\-]\d+\.\d+)\s+([\s\-]\d+\.\d+)",
                    line,
                ):
                    try:
                        _suffix = lines[i + 2 :]
                    except:
                        raise ValueError(
                            "The template file is not a valid GJF file. Please add more blank lines after the section."
                        )
        if prefix != "":
            _prefix = prefix
        if suffix != "":
            _suffix = suffix
        _suffix = "".join(_suffix)
        link, route = _prefix.split("#")
        link = link.replace("\n", " ")
        route = "#" + route.replace("\n", " ")

        _link0 = link0_parser(link)

        if chk:
            _link0[r"%chk"] = f"{os.path.splitext(os.path.basename(_file_path))[0]}.chk"
        if oldchk:
            _link0[r"%oldchk"] = (
                f"{os.path.splitext(os.path.basename(_file_path))[0]}.chk"
            )
        link0_lines = (
            "\n".join([f"{key}={val}" for key, val in _link0.items()]) + "\n"
            if _link0
            else ""
        )
        return (
            link0_lines
            + route
            + "\n\n"
            + f" Title: {self.pure_filename}\n\n"
            + f"{charge if charge else self.charge} {multiplicity if multiplicity else self.multiplicity}\n"
            + "\n".join(
                [
                    f"{Chem.Atom(atom).GetSymbol():10s}{x:10.5f}{y:10.5f}{z:10.5f}"
                    for atom, (x, y, z) in zip(self.atoms, self.coords.m)
                ]
            )
            + "\n\n"
            + _suffix
            + "\n\n"
        )

    def to_GJF_file(
        self,
        file_path: str = None,
        charge: int = None,
        multiplicity: int = None,
        template: str = None,
        prefix: str = "",
        suffix="",
        chk: bool = True,
        oldchk: bool = False,
    ) -> str:
        """
        Write the GJF file.

        Parameters:
            file_path (str):
                The path to write the GJF file. If not specified, will be generated in situ.
            charge (int):
                The forced charge. If specified, will be used to overwrite the charge in the gjf file.
            multiplicity (int):
                The forced multiplicity. If specified, will be used to overwrite the multiplicity in the gjf file.
            template (str):
                path to read a gjf file as a template.
            prefix (str):
                prefix to add to the beginning of the gjf file, priority is higher than template.
            suffix (str):
                suffix to add to the end of the gjf file, priority is higher than template.
            chk (bool):
                If true, add the chk keyword to the link0 section. Will use the file name as the chk file name.
            oldchk (bool):
                If true, add the oldchk keyword to the link0 section. Will use the file name as the chk file name.
        Returns:
            str: The path to the GJF file.
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

    def to_chemdraw(self, file_path: str = None, keep3D=True) -> str:
        """
        Write the ChemDraw file.

        Parameters:
            file_path (str):
                The path to write the ChemDraw file. If not specified, will be generated in situ.
            keep3D (bool):
                Whether to keep the 3D information.
        Returns:
            str: The path to the ChemDraw file.
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

    @property
    def next(self) -> "BaseMolFrameParser":
        """
        Get the next frame.

        Returns:
            BaseMolFrameParser: The next frame.
        """
        if self._next_frame.file_path == self.file_path:
            return self._next_frame
        raise ValueError("No next frame in the file.")

    @property
    def prev(self) -> "BaseMolFrameParser":
        """
        Get the previous frame.
        Returns:
            BaseMolFrameParser: The previous frame.
        """
        if self._prev_frame.file_path == self.file_path:
            return self._prev_frame
        raise ValueError("No previous frame in the file.")

    @classmethod
    def rebuild_parser(
        cls,
        new_mol: Chem.rdchem.Mol,
        path: str = os.path.join(os.getcwd(), "temp.xyz"),
    ) -> "BaseMolFrameParser":
        """
        Create a new parser with the given molecule.

        Parameters:
            new_mol (Chem.rdchem.Mol):
                The new molecule.
            path (str):
                The path to save the new file.
        Returns:
            BaseMolFrameParser: The new parser.
        """
        return cls(
            atoms=[atom.GetAtomicNum() for atom in new_mol.GetAtoms()],
            coords=new_mol.GetConformer().GetPositions() * atom_ureg.angstrom,
            file_path=path,
            frame_content=Chem.MolToXYZBlock(new_mol),
            charge=get_total_charge(new_mol),
            multiplicity=get_total_multiplicity(new_mol),
            bonds=get_bond_pairs(new_mol),
            formal_charges=get_formal_charges(new_mol),
            formal_num_radicals=get_formal_num_radicals(new_mol),
        )

    def replace_substituent(
        self,
        query: Union[str, RdMol],
        replacement: Union[str, RdMol] = None,
        bind_idx: int = None,
        replace_all=False,
        attempt_num=10,
        crowding_threshold=0.75,
        angle_split=10,
        randomSeed=114514,
        start_idx: int = None,
        end_idx: int = None,
    ) -> "BaseMolFrameParser":
        """
        Replace the substituent with the given SMARTS. The substituent is defined by the query_smi,
        and the new substituent is defined by the replacement_smi.

        Parameters:
            query (str | RdMol):
                The SMARTS or Mol object to query the substituent in the original molecule.
            replacement (str | RdMol):
                The SMARTS or Mol object of new substituent.
            bind_idx (int):
                The index of the atom to bind the new substituent. The default is None, which means
                to replace the first legal atom in original molecule.
                If specified, try to replace the legal substruct where the atom in it. User should
                meke sure the atom is legal.
                Detail example in (Repalce Substituent)[Repalce Substituent]
            replace_all (bool):
                If True, replace all the substituent queried in the original molecule.
            attempt_num (int):
                Max attempt times to replace the substituent. Each time a new substituent conformation
                will be used for substitution.
            crowding_threshold (float):
                The threshold of crowding. If the new substituent is too crowded
                `d(a-b) < threshold * (R(a)+R(b))`, the substitution will be rejected.
            angle_split (int):
                Decide how many equal parts of 360° you want to divide. The larger the number the finer
                the rotation will be attempted but the slower the calculation will be.
            randomSeed (int):
                The random seed.
            start_idx (int):
                If both `start_idx` and `end_idx` are specified, simply ignore the `query`, break the
                key between `start_idx` and `end_idx` and replace the base group where `end_idx` is located
            end_idx (int):
                If both `start_idx` and `end_idx` are specified, simply ignore the `query`, break the
                key between `start_idx` and `end_idx` and replace the base group where `end_idx` is located
        Returns:
            BaseMolFrameParser: The new parser.
        """
        rebuild_type = "mod"
        new_mol = attempt_replacement(
            self.rdmol,
            query=query,
            replacement=replacement,
            bind_idx=bind_idx,
            replace_all=replace_all,
            attempt_num=attempt_num,
            crowding_threshold=crowding_threshold,
            angle_split=angle_split,
            randomSeed=randomSeed,
            start_idx=start_idx,
            end_idx=end_idx,
        )
        return self.rebuild_parser(
            new_mol, path=os.path.splitext(self.file_path)[0] + f"_{rebuild_type}.xyz"
        )

    def reset_atom_index(self, mapping_smarts: str) -> "BaseMolFrameParser":
        """
        Reset the atom index of the molecule according to the mapping SMARTS.

        Parameters:
            mapping_smarts (str):
                The SMARTS to query the molecule substructure.
                The queried atoms will be renumbered and placed at the beginning of all atoms according to the order of the atoms in SMARTS. The relative order of the remaining atoms remains unchanged.

        Returns:
            BaseMolFrameParser: The new parser.
        """
        rebuild_type = "reindex"
        smarts = Chem.MolFromSmarts(mapping_smarts)
        if not self.rdmol.HasSubstructMatch(smarts):
            moloplogger.error(f"Failed to match {self.file_path} with {mapping_smarts}")
            raise ValueError("Failed to match")
        mapping = self.rdmol.GetSubstructMatches(smarts)
        if len(mapping) > 1:
            moloplogger.warning(
                f"Multiple matches found in {self.file_path} with {mapping_smarts}"
            )
        mapping = mapping[0]
        rdmol = reset_atom_index(self.rdmol, mapping)
        return self.rebuild_parser(
            rdmol, path=os.path.splitext(self.file_path)[0] + f"_{rebuild_type}.xyz"
        )

    def standard_orient(
        self,
        anchor_list: Sequence[int],
    ) -> "BaseMolFrameParser":
        """
        Depending on the input `idx_list`, `translate_anchor`, `rotate_anchor_to_X`, and `rotate_anchor_to_XY` are executed in order to obtain the normalized oriented molecule.

        Sub-functions:
            - `translate_anchor`: Translate the entire molecule so that the specified atom reaches the origin.
            - `rotate_anchor_to_X`: Rotate the specified second atom along the axis passing through the origin so that it reaches the positive half-axis of the X-axis.
            - `rotate_anchor_to_XY`: Rotate along the axis passing through the origin so that the specified third atom reaches quadrant 1 or 2 of the XY plane.

        Parameters:
            anchor_list (Sequence[int]):
                A Sequence of indices of the atoms to be translated to origin, rotated to X axis, and rotated again to XY face:

                - If length is 1, execute `translate_anchor`
                - If length is 2, execute `translate_anchor` and `rotate_anchor_to_X`
                - If length is 3, execute `translate_anchor`, `rotate_anchor_to_X` and `rotate_anchor_to_XY`
                - If the length of the input `idx_list` is greater than 3, subsequent atomic numbers are not considered.
        Returns:
            BaseMolFrameParser: The new parser.
        """
        rebuild_type = "transform"
        mol = Chem.RWMol(self.rdmol)
        standard_orient(mol, anchor_list)
        return self.rebuild_parser(
            mol, path=os.path.splitext(self.file_path)[0] + f"_{rebuild_type}.xyz"
        )

    @computed_field(
        description="Abstrcact method to check if the current frame is an error. "
        "The details are implemented in the derived classes."
    )
    @property
    def is_error(self) -> bool:
        """
        Abstrcact method to check if the current frame is an error.
        The details are implemented in the derived classes.
        """
        return not self.status.normal_terminated

    @computed_field(description="Check if the molecule is normal.")
    @property
    def is_normal(self) -> bool:
        return not self.is_error

    @computed_field(
        description="Check if the molecule is a TS. Can not check if the molecule "
        "is a TS without frequency information. Thus this function returns False."
    )
    @property
    def is_TS(self) -> bool:
        """
        Check if the molecule is a TS. Can not check if the molecule
        is a TS without frequency information. Thus this function returns False.
        """
        return False

    @property
    def is_optimized(self) -> bool:
        """
        Check if the molecule is optimized.
        """
        return self.geometry_optimization_status.geometry_optimized

    def to_summary_series(self, full: bool = False) -> pd.Series:
        """
        Generate a summary series for the current frame.

        Parameters:
            full (bool):
                Not used but for API compatibility.
        Returns:
            pd.Series: The summary series.
        """
        return pd.Series(
            {
                "parser": self.__class__.__name__,
                "file_path": self.file_path,
                "file_name": self.filename,
                "file_format": self.file_format,
                "frame_index": self.frame_id,
                "charge": self.charge,
                "multiplicity": self.multiplicity,
                "SMILES": self.to_standard_SMILES(),
            }
        )


class BaseQMMolFrameParser(BaseMolFrameParser):
    """
    Base class for QM molecule frames.
    """

    # QM software
    qm_software: str = Field(
        default="",
        description="QM software used to perform the calculation",
    )
    qm_software_version: str = Field(
        default="",
        description="QM software version used to perform the calculation",
    )
    # QM parameters
    keywords: str = Field(
        default="",
        description="Keywords for the QM parameters",
    )
    method: str = Field(
        default="",
        description="QM method used to perform the calculation. e.g. DFT or GFN2-xTB",
    )
    basis: str = Field(
        default="",
        description="Basis set used in the QM calculation, only for DFT calculations",
    )
    functional: str = Field(
        default="",
        description="Functional used in the QM calculation, only for DFT calculations",
    )
    # solvation
    solvent: str = Field(
        default="",
        description="Solvent used in the QM calculation",
    )
    solvent_model: str = Field(
        default="",
        description="Solvent model used in the QM calculation",
    )
    solvent_epsilon: float = Field(
        default=0.0,
        description="Solvent dielectric constant used in the QM calculation",
    )
    solvent_epsilon_infinite: float = Field(
        default=0.0,
        description="Solvent epsilon infinite used in the QM calculation",
    )
    # physical settings
    temperature: PlainQuantity = Field(
        default=298.15 * atom_ureg.K,
        description="Temperature used in the QM calculation, unit is `K`",
    )
    electron_temperature: PlainQuantity = Field(
        default=298.15 * atom_ureg.K,
        description="Electron temperature used in the QM calculation, unit is `K`",
    )
    # QM properties
    forces: Union[PlainQuantity, None] = Field(
        default=np.array([[]]) * atom_ureg.hartree / atom_ureg.bohr,
        description="Forces of each atom, unit is `hartree/bohr`",
    )
    rotation_constants: Union[PlainQuantity, None] = Field(
        default=np.array([[]]) * atom_ureg.Ghertz,
        description="Rotational constants, unit is `GHz`",
    )
    energies: Energies = Field(
        default=Energies(),
        description="Energies",
    )
    thermal_energies: ThermalEnergies = Field(
        default=ThermalEnergies(),
        description="Thermal Energies",
    )
    molecular_orbitals: MolecularOrbitals = Field(
        default=MolecularOrbitals(),
        description="Molecular Orbitals",
    )
    vibrations: Vibrations = Field(default=Vibrations(), description="vibrations")
    charge_spin_populations: ChargeSpinPopulations = Field(
        default=ChargeSpinPopulations(), description="Charge and spin populations"
    )
    polarizability: Polarizability = Field(
        default=Polarizability(), description="Polarizability of the molecule"
    )
    bond_orders: BondOrders = Field(
        default=BondOrders(), description="Bond orders of the molecule"
    )
    total_spin: TotalSpin = Field(
        default=TotalSpin(), description="Total spin of the molecule"
    )
    single_point_properties: SinglePointProperties = Field(
        default=SinglePointProperties(),
        description="Single point properties of the molecule",
    )
    running_time: PlainQuantity = Field(
        default=0.0 * atom_ureg.second,
        description="Running time of the QM calculation, unit is `second`",
    )

    def ts_vibration(self) -> List[BaseMolFrameParser]:
        """
        Generate a list of base block parsers for transition state vibration calculations.

        Returns:
            List[BaseMolFrameParser]: A list of base block parsers for transition state vibration calculations.
        """
        if not self.is_TS:
            raise RuntimeError("This is not a TS")

        block_parsers = []  # Initialize a list of base block parsers

        # Iterate over a list of ratios
        for ratio in (-1.75, -1.5, -1.25, -1, 1, 1.25, 1.5, 1.75):
            # Calculate extreme coordinates based on current ratio
            extreme_coords = (
                self.coords.m
                - self.vibrations.imaginary_vibrations[0].vibration_mode.m * ratio
            )

            try:
                # Convert extreme coordinates to openbabel molecule object
                rdmol = xyz_block_to_rdmol(
                    f"{len(self.atoms)}\n"
                    + f"charge {self.charge} multiplicity {self.multiplicity}\n"
                    + "\n".join(
                        [
                            f"{Chem.Atom(atom).GetSymbol():10s}{x:10.5f}{y:10.5f}{z:10.5f}"
                            for atom, x, y, z in zip(self.atoms, *zip(*extreme_coords))
                        ]
                    ),
                    total_charge=self.charge,
                    total_radical=self.multiplicity - 1,
                )
                # Rebuild parser using the openbabel molecule object
                block_parser = self.rebuild_parser(rdmol)
            except Exception as e:
                moloplogger.warning(
                    f"Failed to rebuild for {self.to_SMILES()} with error {e}"
                )
                continue

            # Check if the molecule satisfies crowding conditions and append it to the list
            block_parsers.append(block_parser)

        return block_parsers

    def ts_vibration_to_SDF_file(self, file_path: str = None) -> str:
        """
        Write the TS vibration calculations to an SDF file.

        Parameters:
            file_path (str):
                The file path. If not specified, will be generated in situ.
        Returns:
            str: The absolute path of the SDF file.
        """
        if file_path is None:
            file_path = self.file_path
        if os.path.isdir(file_path):
            raise IsADirectoryError(f"{file_path} is a directory.")
        file_path = os.path.splitext(file_path)[0] + ".sdf"
        block_parsers = self.ts_vibration()
        with open(file_path, "w") as f:
            f.write("$$$$\n".join([frame.to_SDF_block() for frame in block_parsers]))
        f.close()
        return file_path

    def possible_pre_post_ts(
        self, show_3D=False
    ) -> Tuple[Chem.rdchem.Mol, Chem.rdchem.Mol]:
        """
        This method returns the possible pre- and post-transition state molecules.

        Parameters:
            show_3D (bool): Whether to show 3D coordinates. Defaults to False.

        Returns:
            Tuple[Chem.rdchem.Mol, Chem.rdchem.Mol]: A tuple containing the possible pre- and post-transition state molecules.
        """
        block_parsers = self.ts_vibration()
        if not show_3D:
            block_parsers[0].rdmol.RemoveAllConformers()
            block_parsers[-1].rdmol.RemoveAllConformers()
        return block_parsers[0].rdmol, block_parsers[-1].rdmol

    def to_summary_series(self, full: bool = False) -> pd.Series:
        """
        Generate a summary series for the current frame.

        Parameters:
            full (bool): Whether to include all properties. Defaults to False.

        Returns:
            pd.Series: A summary series for the current frame.
        """
        brief_dict = {
            "parser": self.__class__.__name__,
            "file_path": self.file_path,
            "file_name": self.filename,
            "file_format": self.file_format,
            "version": self.qm_software_version,
            "frame_index": self.frame_id,
            "charge": self.charge,
            "multiplicity": self.multiplicity,
            "SMILES": self.to_standard_SMILES(),
            "keywords": self.keywords,
            "functional": self.functional,
            "basis": self.basis,
            "solvent_model": self.solvent_model,
            "solvent": self.solvent,
            "temperature": self.temperature,
            **self.status.model_dump(),
            "is_TS": self.is_TS,
        }
        if full:
            brief_dict = {
                **brief_dict,
                "rotational_constants": self.rotation_constants,
                **self.energies.model_dump(),
                **self.thermal_energies.model_dump(),
                **self.total_spin.model_dump(),
                **self.molecular_orbitals.model_dump(),
                **self.vibrations.model_dump(),
                **self.polarizability.model_dump(),
                **self.bond_orders.model_dump(),
                **self.single_point_properties.model_dump(),
                **self.geometry_optimization_status.model_dump(),
                "running_time": self.running_time,
            }
        return pd.Series(brief_dict)

    @computed_field
    @property
    def is_error(self) -> bool:
        """
        Abstrcact method to check if the current frame is an error. The details are implemented in the derived classes.
        """
        if self.energies.total_energy is None:
            return True
        return not self.status.normal_terminated

    @computed_field
    @property
    def is_TS(self) -> bool:
        """
        Abstract method to check if the current frame is a transition state. The details are implemented in the derived classes.
        """
        if self.is_error:
            return False
        if not self.is_optimized:
            return False
        return len(self.vibrations.imaginary_idxs) == 1


MolFrameType = TypeVar("MolFrameType", bound=BaseMolFrameParser, covariant=True)
QMMolFrameType = TypeVar("QMMolFrameType", bound=BaseQMMolFrameParser, covariant=True)
