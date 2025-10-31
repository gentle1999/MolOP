"""
Author: TMJ
Date: 2025-07-29 12:36:28
LastEditors: TMJ
LastEditTime: 2025-10-29 14:17:47
Description: 请填写简介
"""

import os
from typing import Any, Dict, Sequence

from pydantic import Field, PrivateAttr, computed_field, field_validator
from typing_extensions import Self

from molop.io.base_models.Bases import BaseDataClassWithUnit


class MemoryStorageMixin(BaseDataClassWithUnit): ...


class DiskStorageMixin(BaseDataClassWithUnit):
    file_path: str = Field(default="", repr=False, exclude=True)
    _allowed_formats_: Sequence[str] = PrivateAttr(default_factory=tuple)

    @field_validator("file_path")
    @classmethod
    def _check_file_path(cls, v: str) -> str:
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

    def to_summary_dict(self, **kwargs) -> Dict[str, Any]:
        return {
            "FilePath": self.file_path,
            "FileFormat": self.file_format,
            **super().to_summary_dict(**kwargs),
        }
