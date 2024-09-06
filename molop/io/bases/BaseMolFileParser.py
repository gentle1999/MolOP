"""
Author: TMJ
Date: 2024-06-20 16:20:28
LastEditors: TMJ
LastEditTime: 2024-06-20 16:47:32
Description: 请填写简介
"""

import os
from typing import Generic, List, Sequence, Tuple, Union

import pandas as pd
from pint.facets.plain import PlainQuantity
from pydantic import (
    Field,
    PrivateAttr,
    computed_field,
    field_validator,
    model_validator,
)
from typing_extensions import Self

from molop.io.bases.BaseMolFrameParser import MolFrameType, QMMolFrameType
from molop.io.bases.DataClasses import BaseDataClassWithUnit
from molop.unit import atom_ureg, unit_transform
from molop.utils.types import RdMol


class BaseMolFileParser(BaseDataClassWithUnit, Generic[MolFrameType]):
    file_path: str = Field(
        default="",
        description="File path",
    )

    charge: int = Field(default=0, description="charge")
    multiplicity: int = Field(default=1, description="multiplicity")
    only_extract_structure: bool = Field(default=False, exclude=True, repr=False)
    only_last_frame: bool = Field(default=False, exclude=True, repr=False)
    _allowed_formats: Tuple[str] = PrivateAttr(default=())
    __frames: List[MolFrameType] = PrivateAttr(default=[])
    __index: int = PrivateAttr(default=0)

    def _parse(self):
        raise NotImplementedError

    def _set_default_units(self):
        pass

    @model_validator(mode="after")
    def __parse_file__(self) -> Self:
        self._parse()
        return self

    @field_validator("file_path", mode="before")
    @classmethod
    def check_file_path(cls, file_path: str) -> str:
        """
        Check if the file exists and is a file, and if its format is allowed.

        Parameters:
            file_path (str): The path of the file to be checked.

        Raises:
            FileNotFoundError: If the file does not exist.
            IsADirectoryError: If the file is a directory.
            ValueError: If the file format is not allowed.
        """
        if file_path == "":
            return ""
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"{file_path} not found.")
        elif not os.path.isfile(file_path):
            raise IsADirectoryError(f"{file_path} is not a file.")
        else:
            _, file_format = os.path.splitext(file_path)
            if file_format not in cls._allowed_formats.default:
                raise ValueError(
                    f"File format {file_format} not in {cls._allowed_formats}"
                )
        return os.path.abspath(file_path)

    @computed_field
    @property
    def file_format(self) -> str:
        """
        return the file format of the file.
        """
        _, file_format = os.path.splitext(self.file_path)
        return file_format

    def _check_path(self, file_path: str = None, format: str = ".xyz"):
        """
        Check the file path and format. If the file path is a directory, will generate a file path with the same name as the
        input file and the specified format.

        Parameters:
            file_path (str):
                The file path. If not specified, will use the file path in the object.
            format (str):
                The file format.

        Returns:
            strThe absolute path of the file.
        """
        if file_path is None:
            return os.path.splitext(self.file_path)[0] + format
        if os.path.isdir(file_path):
            return os.path.join(
                file_path,
                self.pure_filename + format,
            )
        return file_path

    def __iter__(self):
        self.__index = 0
        return self

    def __next__(
        self,
    ) -> MolFrameType:
        if self.__index >= len(self):
            raise StopIteration
        else:
            self.__index += 1
            return self.__frames[self.__index - 1]

    def __getitem__(self, frameID: int) -> MolFrameType:
        return self.__frames[frameID]

    def __len__(self) -> int:
        return len(self.__frames)

    @property
    def file_name(self) -> str:
        """
        Get the file name from the file path.

        Returns:
            str: The file name.
        """
        return os.path.basename(self.file_path)

    @property
    def pure_filename(self) -> str:
        return os.path.splitext(self.file_name)[0]

    @property
    def frames(
        self,
    ) -> List[MolFrameType]:
        """
        Get the list of parsed frames.

        Returns:
            List[MolFrameType]: The list of parsed frames.
        """
        return self.__frames

    def __append__(
        self,
        frame: MolFrameType,
    ):
        """
        Append a frame to the list of frames.

        Parameters:
            frame (MolFrameType): The frame to be appended.
        """
        frame.frame_id = len(self.__frames)
        if frame.frame_id > 0:
            self.__frames[frame.frame_id - 1]._next_frame = frame
            frame._prev_frame = self.__frames[frame.frame_id - 1]
        self.__frames.append(frame)

    def to_XYZ_block(self) -> str:
        """
        Convert all frames to a string representation in XYZ format.

        Returns:
            str: The string representation of the object in XYZ format.
        """
        return "\n".join([frame.to_XYZ_block() for frame in self.__frames])

    def to_SDF_block(self) -> str:
        """
        Convert all frames to a string representation in SDF format.

        Returns:
           str: The string representation in SDF format.
        """
        return "$$$$\n".join([frame.to_SDF_block() for frame in self.__frames])

    def to_GJF_block(
        self,
        file_path: str = None,
        charge: int = None,
        multiplicity: int = None,
        prefix: str = "",
        suffix="",
        template: str = None,
        chk: bool = True,
        oldchk: bool = False,
        frameID=-1,
    ) -> str:
        """
        Convert the frame at the specified index to a Gaussian input block in the Gaussian 16 format (GJF).

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
            str: A modified GJF block.
        """
        return self.__frames[frameID].to_GJF_block(
            file_path=file_path,
            charge=charge,
            multiplicity=multiplicity,
            prefix=prefix,
            suffix=suffix,
            template=template,
            chk=chk,
            oldchk=oldchk,
        )

    def to_XYZ_file(self, file_path: str = None):
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

    def to_SDF_file(self, file_path: str = None):
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
        frameID=-1,
    ):
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
            frameID (int):
                The frame ID to write.
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
                    prefix=prefix,
                    suffix=suffix,
                    template=template,
                    chk=chk,
                    oldchk=oldchk,
                    frameID=frameID,
                )
            )
        f.close()
        return os.path.abspath(_file_path)

    def to_chemdraw(self, file_path: str = None, frameID=-1, keep3D=True):
        """
        Write the ChemDraw file.

        Parameters:
            file_path (str):
                The path to write the ChemDraw file. If not specified, will be generated in situ.
            frameID (int):
                The frame ID to write.
            keep3D (bool):
                Whether to keep the 3D information.
        Returns:
            str: The path to the ChemDraw file.
        """
        return self.__frames[frameID].to_chemdraw(file_path, keep3D=keep3D)

    def geometry_analysis_df(
        self, key_atoms: List[List[int]], precision: int = 1, one_start=False
    ) -> pd.DataFrame:
        """
        Get the geometry infos among the atoms with all frames in the file

        Parameters:
            key_atoms (List[List[int]]):
                A list of list of index of the atoms, starts from 0:

                    - If the length of atom_idxs is 2, the bond length with unit Angstrom between the two atoms will be returned.
                    - If the length of atom_idxs is 3, the angle with unit degree between  the three atoms will be returned.
                    - If the length of atom_idxs is 4, the dihedral angle with unit degree between the four atoms will be returned.
            precision (int):
                The precision of the geometry analysis. Default is 1. e.g. 1 means 1.0001 ==> 1.0
            one_start (bool):
                If true, consider atom index starts from 1, so let index value subtracts 1 for all the atoms

        Returns:
            pd.DataFrame:
                The columns of the DataFrame is the atom_idxs, and the index are the frames.
        """
        values = {
            "-".join([str(idx) for idx in atom_idxs]): [
                round(molblock.geometry_analysis(atom_idxs, one_start), precision)
                for molblock in self.__frames
            ]
            for atom_idxs in key_atoms
        }
        return pd.DataFrame(values)

    def geometry_analysis(
        self,
        key_atoms: List[List[int]],
        file_path: str = None,
        precision: int = 1,
        one_start=False,
    ):
        """
        Get the geometry infos among the atoms with all frames in the file

        Parameters:
            key_atoms (List[List[int]]):
                A list of list of index of the atoms, starts from 0

                    - If the length of atom_idxs is 2, the bond length with unit Angstrom between the two atoms will be returned.
                    - If the length of atom_idxs is 3, the angle with unit degree between  the three atoms will be returned.
                    - If the length of atom_idxs is 4, the dihedral angle with unit degree between the four atoms will be returned.
            file_path (str):
                The path of the csv file to be saved. If None, the file will be saved in the same directory of the file_path.
            precision (int):
                The precision of the geometry analysis. Default is 1. e.g. 1 means 1.0001 ==> 1.0
            one_start (bool):
                If true, consider atom index starts from 1, so let index value subtracts 1 for all the atoms

        Returns:
            str: file_path
        """
        df = self.geometry_analysis_df(
            key_atoms, precision=precision, one_start=one_start
        )
        if file_path is None:
            file_path = self.file_path
        if os.path.isdir(file_path):
            raise IsADirectoryError(f"{file_path} is a directory.")
        file_path = os.path.splitext(file_path)[0] + "_geometry_analysis.csv"
        df.to_csv(file_path)
        return file_path

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
        frameID=-1,
    ) -> MolFrameType:
        """
        Replace the substituent with the given SMARTS. The substituent is defined by the query_smi, and the new substituent is defined by the replacement_smi.

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
            frameID int:
                The frame ID to replace.
        Returns:
            The new parser.
        """
        return self.__frames[frameID].replace_substituent(
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

    def reset_atom_index(
        self,
        mapping_smarts: str,
        frameID: int = -1,
    ) -> MolFrameType:
        """
        Reset the atom index of the molecule according to the mapping SMARTS.

        Parameters:
            mapping_smarts (str):
                The SMARTS to query the molecule substructure.
                The queried atoms will be renumbered and placed at the beginning of all atoms according to the order of the atoms in SMARTS. The relative order of the remaining atoms remains unchanged.
            frameID (int):
                The frame ID to reset.
        Returns:
            MolFrameType: The new parser.
        """
        return self.__frames[frameID].reset_atom_index(mapping_smarts)

    def standard_orient(
        self,
        anchor_list: Sequence[int],
        frameID: int = -1,
    ) -> MolFrameType:
        """
        Depending on the input `idx_list`, `translate_anchor`, `rotate_anchor_to_axis`, and `rotate_anchor_to_plane` are executed in order to obtain the normalized oriented molecule.

        Sub-functions:
            - `translate_anchor`: Translate the entire molecule so that the specified atom reaches the origin.
            - `rotate_anchor_to_axis`: Rotate the specified second atom along the axis passing through the origin so that it reaches the positive half-axis of the X-axis.
            - `rotate_anchor_to_plane`: Rotate along the axis passing through the origin so that the specified third atom reaches quadrant 1 or 2 of the XY plane.

        Parameters:
            anchor_list (Sequence[int]):
                A list of indices of the atoms to be translated to origin, rotated to X axis, and rotated again to XY face:

                - If length is 1, execute `translate_anchor`
                - If length is 2, execute `translate_anchor` and `rotate_anchor_to_axis`
                - If length is 3, execute `translate_anchor`, `rotate_anchor_to_axis` and `rotate_anchor_to_plane`
                - If the length of the input `idx_list` is greater than 3, subsequent atomic numbers are not considered.
            frameID (int):
                The frame ID to standardize orientation.

        Returns:
            MolFrameType: The new parser.
        """
        return self.__frames[frameID].standard_orient(anchor_list)

    def to_summary_df(self, full: bool = False, with_units: bool = True) -> pd.DataFrame:
        """
        Get the summary information of the parser.

        Parameters:
            full (bool):
                If True, include the full information of the parser.

        Returns:
            pd.DataFrame:
                A pandas DataFrame containing the summary information of the parser.
        """
        self.recover_structures()
        return pd.concat(
            [frame.to_summary_series(full, with_units) for frame in self.__frames], axis=1
        ).T

    def recover_structures(self) -> List[str]:
        """
        Recover the structures of the frames.

        Returns:
            List[str]: A list of SMILES strings of the frames.
        """
        return [frame.to_canonical_SMILES() for frame in self.frames]


class BaseQMMolFileParser(BaseMolFileParser[QMMolFrameType]):
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
    running_time: PlainQuantity = Field(
        default=0.0 * atom_ureg.second,
        description="Running time of the QM calculation, unit is `second`",
    )

    def _set_default_units(self):
        self.temperature = unit_transform(self.temperature, "K")
        self.electron_temperature = unit_transform(self.electron_temperature, "K")
        self.running_time = unit_transform(self.running_time, "s")

    @property
    def sort_by_optimization(self) -> List[QMMolFrameType]:
        """
        Sort the frames by the optimization status. The closer the frame is to the optimized state, the higher the priority.

        Returns:
            List[QMMolFrameType]: A list of frames sorted by the optimization status.
        """
        frames = self.frames
        return sorted(frames, key=lambda frame: frame.geometry_optimization_status)

    @property
    def closest_optimized_frame(self) -> QMMolFrameType:
        """
        Get the frame with the closest optimization status to the optimized state.

        Returns:
            QMMolFrameType: The frame with the closest optimization status to the optimized state.
        """
        return self.sort_by_optimization[0]
