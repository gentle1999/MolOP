"""
Author: TMJ
Date: 2024-01-07 13:50:55
LastEditors: TMJ
LastEditTime: 2024-01-25 22:56:04
Description: 请填写简介
"""
import os
from typing import List, Tuple, TypeVar, Union
import pandas as pd

from molop.io.bases.molblock_base import BaseBlockParser, QMBaseBlockParser
from molop.io.coords_file.GJFBlockParser import GJFBlockParser
from molop.io.coords_file.SDFBlockParser import SDFBlockParser
from molop.io.coords_file.XYZBlockParser import XYZBlockParser
from molop.io.qm_file.G16FCHKBlockParser import G16FCHKBlockParser
from molop.io.qm_file.G16IRCBlockParser import G16IRCBlockParser
from molop.io.qm_file.G16LOGBlockParser import G16LOGBlockParser
from molop.io.qm_file.XTBOUTBlockParser import XTBOUTBlockParser

BlockType = Union[
    GJFBlockParser,
    SDFBlockParser,
    XYZBlockParser,
    G16LOGBlockParser,
    G16IRCBlockParser,
    G16FCHKBlockParser,
    XTBOUTBlockParser,
]


class BaseFileParser:
    """
    Base class for multi-frame parsers.
    """

    _file_path: str
    __frames: List[BlockType]
    __index: int
    _allowed_formats: Tuple[str]

    def __init__(self, file_path: str) -> None:
        self._file_path = os.path.abspath(file_path)
        self.__frames = []
        self.__index: int = 0

    def _check_formats(self, file_path: str) -> None:
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"{file_path} not found.")
        elif not os.path.isfile(file_path):
            raise IsADirectoryError(f"{file_path} is not a file.")
        else:
            _, file_format = os.path.splitext(file_path)
            if file_format not in self._allowed_formats:
                raise ValueError(
                    f"File format {file_format} not in {self._allowed_formats}"
                )
        self._file_format = file_format

    def __iter__(self):
        self.__index = 0
        return self

    def __next__(
        self,
    ) -> BlockType:
        if self.__index >= len(self):
            raise StopIteration
        else:
            self.__index += 1
            return self.__frames[self.__index - 1]

    def __getitem__(self, frameID: int) -> BlockType:
        return self.__frames[frameID]

    def __len__(self) -> int:
        return len(self.__frames)

    def __str__(self) -> str:
        return f"{self.__class__.__name__}({os.path.basename(self._file_path)})"

    def __hash__(self) -> int:
        return hash(str(self))

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({str(self)})"

    @property
    def file_path(self) -> str:
        return self._file_path

    @property
    def file_name(self) -> str:
        return os.path.basename(self._file_path)

    @property
    def frames(
        self,
    ) -> List[BlockType]:
        return self.__frames

    def append(
        self,
        frame: BlockType,
    ) -> None:
        if not issubclass(type(frame), BaseBlockParser):
            raise TypeError(f"{type(frame)} is not a subclass of {BaseBlockParser}")
        frame._frameID = len(self.__frames)
        if frame._frameID > 0:
            self.__frames[frame._frameID - 1]._next_block = frame
        self.__frames.append(frame)

    def to_XYZ_block(self) -> str:
        return "\n".join([frame.to_XYZ_block() for frame in self.__frames])

    def to_SDF_block(self) -> str:
        return "$$$$\n".join([frame.to_SDF_block() for frame in self.__frames])

    def to_GJF_block(
        self,
        charge: int = None,
        multiplicity: int = None,
        prefix: str = f"# g16 gjf \n",
        suffix="",
        template: str = None,
        frameID=-1,
    ) -> str:
        """Only extract one frame."""
        return self.__frames[frameID].to_GJF_block(
            charge=charge,
            multiplicity=multiplicity,
            prefix=prefix,
            suffix=suffix,
            template=template,
        )

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

    def to_GJF_file(
        self,
        file_path: str = None,
        charge: int = None,
        multiplicity: int = None,
        prefix: str = f"# g16 gjf \n",
        suffix="\n\n",
        template: str = None,
        frameID=-1,
    ):
        """Only extract one frame."""
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
                    prefix=prefix,
                    suffix=suffix,
                    template=template,
                    frameID=frameID,
                )
            )
        f.close()
        return file_path

    def to_chemdraw(self, file_path: str = None, frameID=-1, keep3D=False):
        return self.__frames[frameID].to_chemdraw(file_path, keep3D=keep3D)

    def summary(self):
        print(
            f"file path: {self._file_path}\n"
            + f"frame num: {len(self)}\n"
            + f"first SMILES: {self[0].to_SMILES()}\n"
            + f"last SMILES: {self[-1].to_SMILES()}\n"
        )

    def geometry_analysis_df(
        self, key_atoms: List[List[int]], precision: int = 1, one_start=False
    ):
        """
        Get the geometry infos among the atoms with all frames in the file

        Parameters:
            key_atoms List[List[int]]:
                A list of list of index of the atoms, starts from 0
                    If the length of atom_idxs is 2, the bond length with unit Angstrom between the two atoms will be returned.

                    If the length of atom_idxs is 3, the angle with unit degree between  the three atoms will be returned.

                    If the length of atom_idxs is 4, the dihedral angle with unit degree between the four atoms will be returned.
            precision int:
                The precision of the geometry analysis. Default is 1. e.g. 1 means 1.0001 ==> 1.0
            one_start bool:
                If true, consider atom index starts from 1, so let index value subtracts 1 for all the atoms

        Returns:
            A pandas DataFrame:
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
            key_atoms List[List[int]]:
                A list of list of index of the atoms, starts from 0
                    If the length of atom_idxs is 2, the bond length with unit Angstrom between the two atoms will be returned.

                    If the length of atom_idxs is 3, the angle with unit degree between  the three atoms will be returned.

                    If the length of atom_idxs is 4, the dihedral angle with unit degree between the four atoms will be returned.
            file_path str:
                The path of the csv file to be saved. If None, the file will be saved in the same directory of the file_path.
            precision int:
                The precision of the geometry analysis. Default is 1. e.g. 1 means 1.0001 ==> 1.0
            one_start bool:
                If true, consider atom index starts from 1, so let index value subtracts 1 for all the atoms

        Returns:
            file_path
        """
        df = self.geometry_analysis_df(
            key_atoms, precision=precision, one_start=one_start
        )
        if file_path is None:
            file_path = self._file_path
        if os.path.isdir(file_path):
            raise IsADirectoryError(f"{file_path} is a directory.")
        file_path = os.path.splitext(file_path)[0] + "_geometry_analysis.csv"
        df.to_csv(file_path)
        return file_path


class BaseQMFileParser(BaseFileParser):
    """
    Base class for QM multi-frame parsers.
    """

    _parameter_comment: str
    _show_progress: bool
    _only_extract_structure: bool
    _only_last_frame: bool
    _version: str

    def __init__(
        self,
        file_path: str,
        show_progress=False,
        only_extract_structure=False,
        only_last_frame=False,
    ) -> None:
        super().__init__(file_path)
        self._parameter_comment: str = None
        self._show_progress: bool = show_progress
        self._only_extract_structure: bool = only_extract_structure
        self._only_last_frame = only_last_frame
        self._version = None

    @property
    def parameter_comment(self) -> str:
        return self._parameter_comment

    @property
    def version(self) -> str:
        return self._version
