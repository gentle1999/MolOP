"""
Author: TMJ
Date: 2024-01-07 13:50:55
LastEditors: TMJ
LastEditTime: 2024-01-11 09:56:24
Description: 请填写简介
"""
import os
from typing import List, Union, Tuple, TypeVar

from molop.io.bases.molblock_base import BaseBlockParser, QMBaseBlockParser
from molop.io.coords_file.GJFBlockParser import GJFBlockParser
from molop.io.coords_file.SDFBlockParser import SDFBlockParser
from molop.io.coords_file.XYZBlockParser import XYZBlockParser
from molop.io.qm_file.G16IRCBlockParser import G16IRCBlockParser
from molop.io.qm_file.G16LOGBlockParser import G16LOGBlockParser
from molop.io.qm_file.XTBOUTBlockParser import XTBOUTBlockParser

BlockType = Union[
    GJFBlockParser,
    SDFBlockParser,
    XYZBlockParser,
    G16LOGBlockParser,
    G16IRCBlockParser,
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
        frame._file_path = self._file_path
        frame._frameID = len(self.__frames)
        if frame._frameID > 0:
            self.__frames[frame._frameID - 1]._next_block = frame
        self.__frames.append(frame)

    def to_XYZ_block(self) -> str:
        return "\n".join([frame.to_XYZ_block() for frame in self.__frames])

    def to_SDF_block(self) -> str:
        return "$$$$\n".join([frame.to_SDF_block() for frame in self.__frames])

    def to_XYZ_file(self, file_path: str):
        with open(file_path, "w") as f:
            f.write(self.to_XYZ_block())
        f.close()

    def to_SDF_file(self, file_path: str):
        with open(file_path, "w") as f:
            f.write(self.to_SDF_block())
        f.close()

    def summary(self):
        print(
            f"file path: {self._file_path}\n"
            + f"frame num: {len(self)}\n"
            + f"first SMILES: {self[0].to_SMILES()}\n"
            + f"last SMILES: {self[-1].to_SMILES()}\n"
        )


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
