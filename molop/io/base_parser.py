"""
Author: TMJ
Date: 2023-10-30 14:40:42
LastEditors: TMJ
LastEditTime: 2023-10-30 17:39:15
Description: 请填写简介
"""
import os


class BaseParser:
    """
    Base class for all parsers.
    """

    def __init__(self, file_path):
        self.file_path = file_path
        self._atoms = []
        self._coords = []
        self._charge = None
        self._multi = None

    def parse(self):
        """
        Parse the file.
        """
        raise NotImplementedError

    @property
    def coords(self):
        if len(self._coords) == 0:
            return None
        return self._coords

    @property
    def charge(self):
        return self._charge

    @property
    def multi(self):
        return self._multi

    @property
    def atoms(self):
        if len(self._atoms) == 0:
            return None
        return self._atoms

    @staticmethod
    def get_property_methods(obj):
        """
        Get all methods described by "property"
        """
        # TODO 未确认是否会放入未被填充的元素
        property_methods = []
        for name in dir(type(obj)):
            member = getattr(type(obj), name)
            if (
                isinstance(member, property)
                and callable(member.fget)
                and member.fget.__get__(obj, type(obj)) is not None
            ):
                property_methods.append(member.fget)
        return property_methods

    def to_XYZ_file(self, file_path=None):
        if file_path is None:
            file_path = os.path.splitext(self.file_path)[0] + ".xyz"
        props = self.get_property_methods(self)
        props = [prop for prop in props if prop.__name__ not in ("atoms", "coords")]
        comments = "; ".join(
            [
                f"{prop.__name__}={self.__getattribute__(prop.__name__)}"
                for prop in props
            ]
        )
        coords = [self.coords[i : i + 3] for i in range(0, len(self.coords), 3)]
        with open(file_path, "w") as f:
            f.write(f"{len(self.atoms)}\n")
            f.write(f"Created by MolOP; {comments}\n")
            for atom, coord in zip(self.atoms, coords):
                f.write(
                    f"{atom:3} {coord[0]:>15.8f} {coord[1]:>15.8f} {coord[2]:>15.8f}\n"
                )
        f.close()
