"""
Author: TMJ
Date: 2023-10-30 14:40:42
LastEditors: TMJ
LastEditTime: 2023-11-01 13:22:48
Description: 请填写简介
"""
from abc import ABC
import os
from molop.mol import Molecule
from molop.utils.types import RdMol
from molop.utils import geometry, structure

class BaseParser(ABC):
    """
    Base class for all parsers.
    Only One frame.
    """

    def __init__(self, block):
        self._block = block
        self._FilePath_attach_info = None
        self._frameID_attach_info = None
        self._atoms_attach_info = []
        self._coords_attach_info = []
        self._charge_attach_info = None
        self._multi_attach_info = None
        self._bond_pairs_attach_info = None
        # TODO local spin
        # TODO local charge

    def _parse(self):
        """
        Parse the file.
        """
        raise NotImplementedError

    @property
    def coords(self):
        if len(self._coords_attach_info) == 0:
            return None
        return self._coords_attach_info

    @property
    def charge(self):
        return self._charge_attach_info

    @property
    def multi(self):
        return self._multi_attach_info

    @property
    def atoms(self):
        if len(self._atoms_attach_info) == 0:
            return None
        return self._atoms_attach_info

    def __str__(self) -> str:
        return f"{os.path.basename(self._FilePath_attach_info)}({self._frameID_attach_info})"

    def __hash__(self) -> int:
        return hash(str(self))

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({str(self)})"

    def _get_attach_infos(self):
        """
        Get all methods described by "property"
        """
        return [
            attr
            for attr in dir(self)
            if attr.endswith("_attach_info") and self.__getattribute__(attr) is not None
        ]

    def to_XYZ_block(self):
        attrs = [
            attr
            for attr in self._get_attach_infos()
            if attr not in ("_atoms_attach_info", "_coords_attach_info")
        ]
        comments = "; ".join(
            [f"{attr.split('_')[1]}=={self.__getattribute__(attr)}" for attr in attrs]
        )
        coords = [self.coords[i : i + 3] for i in range(0, len(self.coords), 3)]
        block = f"{len(self.atoms)}\nCreated by MolOP; {comments}\n"
        for atom, coord in zip(self.atoms, coords):
            block += f"{atom:3} {coord[0]:>15.8f} {coord[1]:>15.8f} {coord[2]:>15.8f}\n"
        return block

    def to_XYZ_file(self, file_path=None):
        if file_path is None:
            try:
                file_path = os.path.splitext(self._FilePath_attach_info)[0] + ".xyz"
            except:
                raise NotImplementedError("file_path not be defined")
        block = self.to_XYZ_block()
        with open(file_path, "w") as f:
            f.write(block)
        f.close()

    def to_RdMol(self) -> RdMol:
        if all(self._bond_pairs_attach_info, self._atoms_attach_info, self._coords_attach_info):
            mol = structure.bond_pairs_2_mol(self._atoms_attach_info, self.get_bond_pairs())

            return
        

    def to_MolOP_Molecule(self) -> Molecule:
        pass

    def get_bond_pairs(self):
        if self._bond_pairs_attach_info:
            return self._bond_pairs_attach_info
        else:
            mol = self.to_RdMol()



class MultiFrameBaseParser(ABC):
    """
    The base class for multiple frames files
    """

    def __init__(self, file_path) -> None:
        self._FilePath_attach_info = file_path
        self._blocks = []
        self.__frames = []
        self.__index = 0

    def __iter__(self):
        self.__index = 0
        return self

    def __next__(self) -> BaseParser:
        if self.__index >= len(self):
            raise StopIteration
        else:
            self.__index += 1
            return self.__frames[self.__index - 1]

    def __len__(self):
        return len(self.__frames)

    def __getitem__(self, index) -> BaseParser:
        return self.__frames[index]

    def __repr__(self):
        return (
            f"{self.__class__.__name__}({os.path.basename(self._FilePath_attach_info)})"
        )

    def _parse_meta(self):
        raise NotImplementedError

    def _split_frames(self):
        raise NotImplementedError

    def _parse(self):
        raise NotImplementedError

    @property
    def file_path(self):
        return self._FilePath_attach_info

    def _get_attach_infos(self):
        """
        Get all methods described by "property"
        """
        return [
            attr
            for attr in dir(self)
            if attr.endswith("_attach_info") and self.__getattribute__(attr) is not None
        ]

    def __set_frame_meta(self, frame: BaseParser):
        """
        设置属于整个文件的统一的元数据（例如基组等只会在文件头部出现一次的属性）
        """
        if not isinstance(frame, BaseParser):
            raise TypeError(f"{frame} is not a BaseParser")
        for attr in self._get_attach_infos():
            frame.__setattr__(attr, self.__getattribute__(attr))

    def _append(self, frame: BaseParser):
        if not isinstance(frame, BaseParser):
            raise TypeError(f"{frame} is not a BaseParser")
        frame._frameID_attach_info = len(self)
        self.__set_frame_meta(frame)
        self.__frames.append(frame)

    def to_XYZ_block(self):
        block = "".join([frame.to_XYZ_block() for frame in self])
        return block

    def to_XYZ_file(self, file_path):
        if file_path is None:
            try:
                file_path = os.path.splitext(self._file_path_attach_info)[0] + ".xyz"
            except:
                raise NotImplementedError("file_path not be defined")
        with open(file_path, "w") as f:
            f.write(self.to_XYZ_block())
        f.close()
