"""
Author: TMJ
Date: 2024-01-07 13:50:55
LastEditors: TMJ
LastEditTime: 2024-01-07 14:06:43
Description: 请填写简介
"""
import os
from typing import Any, List, Tuple, Union, str, Literal

from .molblock_base import MolBlock, BaseBlockParser, QMBaseBlockParser


class BaseFileParser:
    """
    Base class for multi-frame parsers.
    """

    _file_path: str
    _blocks: List[str] = []
    __frames: List[BaseBlockParser] = []
    __index: int = 0

    def __init__(self, file_path) -> None:
        self._file_path = file_path

    def __iter__(self):
        self.__index = 0
        return self

    def __next__(self) -> BaseBlockParser:
        if self.__index >= len(self):
            raise StopIteration
        else:
            self.__index += 1
            return self.__frames[self.__index - 1]

    def __getitem__(self, frameID: int) -> BaseBlockParser:
        return self.__frames[frameID]

    def __len__(self) -> int:
        return len(self._blocks)

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
    def blocks(self) -> List[str]:
        return self._blocks

    @property
    def frames(self) -> List[BaseBlockParser]:
        return self.__frames
    
class BaseQMFileParser(BaseFileParser):
    """
    Base class for QM multi-frame parsers.
    """

    _state: Literal["Normal Termination", "Convergence"]

    def __init__(self, file_path) -> None:
        super().__init__(file_path)
    
