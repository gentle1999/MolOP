"""
Author: TMJ
Date: 2024-01-07 13:50:55
LastEditors: TMJ
LastEditTime: 2024-01-11 09:56:24
Description: 请填写简介
"""
import os
from typing import List, Union

from molop.io.bases.molblock_base import BaseBlockParser, QMBaseBlockParser
from molop.io.coords_file.GJFBlockParser import GJFBlockParser
from molop.io.coords_file.SDFBlockParser import SDFBlockParser
from molop.io.coords_file.XYZBlockParser import XYZBlockParser
from molop.io.qm_file.G16LOGBlockParser import G16LOGBlockParser


class BaseFileParser:
    """
    Base class for multi-frame parsers.
    """

    _file_path: str
    __frames: List[
        Union[
            BaseBlockParser,
            GJFBlockParser,
            SDFBlockParser,
            XYZBlockParser,
            G16LOGBlockParser,
        ],
    ]
    __index: int

    def __init__(self, file_path: str) -> None:
        self._file_path = file_path
        self.__frames = []
        self.__index: int = 0

    def __iter__(self):
        self.__index = 0
        return self

    def __next__(
        self,
    ) -> Union[
        BaseBlockParser,
        GJFBlockParser,
        SDFBlockParser,
        XYZBlockParser,
        G16LOGBlockParser,
    ]:
        if self.__index >= len(self):
            raise StopIteration
        else:
            self.__index += 1
            return self.__frames[self.__index - 1]

    def __getitem__(
        self, frameID: int
    ) -> Union[
        BaseBlockParser,
        GJFBlockParser,
        SDFBlockParser,
        XYZBlockParser,
        G16LOGBlockParser,
    ]:
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
    ) -> List[
        Union[
            BaseBlockParser,
            GJFBlockParser,
            SDFBlockParser,
            XYZBlockParser,
            G16LOGBlockParser,
        ]
    ]:
        return self.__frames

    def append(
        self,
        frame: Union[
            BaseBlockParser,
            GJFBlockParser,
            SDFBlockParser,
            XYZBlockParser,
            G16LOGBlockParser,
        ],
    ) -> None:
        if not issubclass(type(frame), BaseBlockParser):
            raise TypeError(f"{type(frame)} is not a subclass of {BaseBlockParser}")
        frame._file_path = self._file_path
        frame._frameID = len(self.__frames)
        if frame._frameID > 0:
            self.__frames[frame._frameID - 1]._next_block = frame
        self.__frames.append(frame)


class BaseQMFileParser(BaseFileParser):
    """
    Base class for QM multi-frame parsers.
    """

    _parameter_comment: str
    _show_progress: bool
    _only_extract_structure: bool

    def __init__(
        self, file_path: str, show_progress=False, only_extract_structure=False
    ) -> None:
        super().__init__(file_path)
        self._parameter_comment: str = None
        self._show_progress: bool = show_progress
        self._only_extract_structure: bool = only_extract_structure

    @property
    def parameter_comment(self) -> str:
        return self._parameter_comment
