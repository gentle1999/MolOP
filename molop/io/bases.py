import os
from abc import ABC, abstractmethod
from typing import Any, List, str, Tuple, Union

from rdkit import Chem


class MolBlock(ABC):
    """
    Abstract class for a molecule block.
    """

    _atoms: List[Union[str, int]]
    _bonds: List[Tuple[int, int, int]]
    _charge: int
    _multiplicity: int
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
        return self._atoms

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

    @abstractmethod
    def to_XYZ_block(self) -> str:
        raise NotImplementedError

    @abstractmethod
    def to_SDF_block(self) -> str:
        raise NotImplementedError

    @abstractmethod
    def to_SMILES(self) -> str:
        raise NotImplementedError

    @abstractmethod
    def to_InChI(self) -> str:
        raise NotImplementedError

    @abstractmethod
    def to_RdMol(self) -> Chem.rdchem.Mol:
        raise NotImplementedError


class BaseBlockParser(MolBlock):
    """
    Base class for block parsers.
    """
    _file_path : str
    _frameID : int = 0

    def __init__(self, block: str) -> None:
        self._block = block

    @property
    def block(self) -> str:
        return self._block

    def __str__(self) -> str:
        return f"{self.__class__.__name__}({os.path.basename(self._file_path)})[{self._frameID}]"

    def __hash__(self) -> int:
        return hash(str(self))

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({str(self)})"
    

class BaseFileParser:
    """
    Base class for multi-frame parsers.
    """
    _file_path : str
    _blocks : List[str] = []
    _frames : List[BaseBlockParser] = []
    __index : int = 0

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
        return self._frames[frameID]
        
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
        return self._frames



    
    
