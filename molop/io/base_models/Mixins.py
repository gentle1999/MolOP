"""
Author: TMJ
Date: 2025-07-29 12:36:28
LastEditors: TMJ
LastEditTime: 2025-12-10 23:58:55
Description: 请填写简介
"""

import os
from typing import TYPE_CHECKING, Any, Dict, Protocol, Sequence

from pint.facets.numpy.quantity import NumpyQuantity
from pydantic import Field, PrivateAttr, computed_field, field_validator
from rdkit import Chem
from typing_extensions import Self

from molop.io.base_models.Bases import BaseDataClassWithUnit
from molop.io.patterns.G16Patterns import options_parser


class MemoryStorageMixin(BaseDataClassWithUnit): ...


class DiskStorageMixin(BaseDataClassWithUnit):
    file_path: str = Field(default="", repr=False, exclude=True)
    _allowed_formats_: Sequence[str] = PrivateAttr(default_factory=tuple)

    @field_validator("file_path")
    @classmethod
    def _check_file_path(cls, v: str) -> str:
        """
        Validate the file path.

        Parameters:
            v (str): The file path to validate.

        Returns:
            str: The validated file path.

        Raises:
            ValueError: If the file path does not exist or is not a file.
        """
        if not os.path.exists(v):
            raise ValueError(f"File {v} does not exist.")
        if not os.path.isfile(v):
            raise ValueError(f"Path {v} is not a file.")
        for fmt in cls._allowed_formats_.default:  # type: ignore
            if v.endswith(fmt):
                return v
        raise ValueError(
            f"File {v} has an invalid format. "
            f"Allowed formats: {cls._allowed_formats_.default}."  # type: ignore
        )

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

    def __ge__(self, other: Self) -> bool:
        return self.file_path >= other.file_path

    def __gt__(self, other: Self) -> bool:
        return self.file_path > other.file_path

    def __le__(self, other: Self) -> bool:
        return self.file_path <= other.file_path

    def __lt__(self, other: Self) -> bool:
        return self.file_path < other.file_path

    def __eq__(self, other: Self) -> bool:
        return self.file_path == other.file_path

    def to_summary_dict(self, **kwargs) -> Dict[tuple, Any]:
        return {
            ("DiskStorage", "FilePath"): self.file_path,
            ("DiskStorage", "FileFormat"): self.file_format,
            **super().to_summary_dict(**kwargs),  # type: ignore
        }


class ChemFrameProtocol(Protocol):
    atoms: list[int]
    coords: NumpyQuantity
    charge: int
    multiplicity: int
    pure_filename: str


if TYPE_CHECKING:

    class _ChemFrameProtocol(ChemFrameProtocol, DiskStorageMixin): ...
else:

    class _ChemFrameProtocol(DiskStorageMixin): ...


class DiskStorageWithFrameMixin(_ChemFrameProtocol):
    def to_GJF_block(
        self,
        options: str = "",
        route: str = "#p",
        title_card: str | None = None,
        suffix: str = "",
        chk: bool = False,
        oldchk: bool | str = False,
        **kwargs,
    ) -> str:
        """
        Generate a GJF block for the frame.

        Parameters:
            options (str, optional): The options for the GJF block. Defaults to "".
            route (str, optional): The route for the GJF block. Defaults to "#p".
            title_card (str | None, optional): The title card for the GJF block. Defaults to None and use the pure filename.
            suffix (str, optional): The suffix for the GJF block. Defaults to "".
            chk (bool, optional): Whether to include a checkpoint file. Defaults to False.
            oldchk (bool | str, optional): Whether to include an old checkpoint file. Defaults to False.

        Returns:
            str: The GJF block for the frame.
        """
        _options = options_parser(options)
        if chk:
            _options[r"%chk"] = f"{self.pure_filename}.chk"
        if oldchk:
            _options[r"%oldchk"] = f"{self.pure_filename}.chk"
        options_lines = (
            "\n".join([f"{key}={val}" for key, val in _options.items()]) + "\n"
            if _options
            else ""
        )
        return (
            options_lines
            + route
            + "\n\n"
            + f"{title_card or self.pure_filename}\n\n"
            + f"{self.charge} {self.multiplicity}\n"
            + "\n".join(
                [
                    f"{Chem.Atom(atom).GetSymbol():10s}{x:18.10f}{y:18.10f}{z:18.10f}"
                    for atom, (x, y, z) in zip(self.atoms, self.coords.m, strict=True)
                ]
            )
            + "\n\n"
            + suffix
            + "\n\n"
        )
