"""
Author: TMJ
Date: 2025-07-29 12:36:28
LastEditors: TMJ
LastEditTime: 2025-12-16 14:56:32
Description: 请填写简介
"""

import os
from abc import abstractmethod
from collections.abc import Sequence
from typing import TYPE_CHECKING, Any, Literal, Protocol, overload

from pydantic import Field, computed_field
from typing_extensions import Self

from molop.io.base_models.Bases import BaseDataClassWithUnit


class MemoryStorageMixin(BaseDataClassWithUnit): ...


class DiskStorageMixin(BaseDataClassWithUnit):
    file_path: str = Field(default="", repr=False)

    @computed_field()  # type: ignore[misc]
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

    def __eq__(self, other: Self) -> bool:  # type: ignore[override]
        return self.file_path == other.file_path

    def to_summary_dict(self, **kwargs) -> dict[tuple, Any]:
        return {
            ("DiskStorage", "FilePath"): self.file_path,
            ("DiskStorage", "FileFormat"): self.file_format,
            **super().to_summary_dict(**kwargs),  # type: ignore
        }


class FileProtocol(Protocol):
    frames: Sequence


if TYPE_CHECKING:

    class _FileProtocol(FileProtocol, BaseDataClassWithUnit): ...
else:

    class _FileProtocol(BaseDataClassWithUnit): ...


class FileMixin(_FileProtocol):
    @overload
    def _render(
        self,
        frameID: Sequence[int] | int | Literal["all"] = -1,
        embed_in_one_file: Literal[True] = True,
        **kwargs,
    ) -> str: ...
    @overload
    def _render(
        self,
        frameID: Sequence[int] | int | Literal["all"] = -1,
        embed_in_one_file: Literal[False] = False,
        **kwargs,
    ) -> list[str]: ...
    def _render(
        self,
        frameID: Sequence[int] | int | Literal["all"] = -1,
        embed_in_one_file: bool = True,
        **kwargs,
    ) -> str | list[str]:
        if isinstance(frameID, int):
            frameIDs = [frameID if frameID >= 0 else len(self.frames) + frameID]
        elif frameID == "all":
            frameIDs = list(range(len(self.frames)))
        elif isinstance(frameID, Sequence):
            frameIDs = []
            for i in frameID:
                assert isinstance(i, int), "frameID should be a sequence of integers"
                frameIDs.append(i)
        if embed_in_one_file:
            return self._render_frames_in_one_file(frameIDs, **kwargs)
        else:
            return self._render_frames(frameIDs, **kwargs)

    @abstractmethod
    def _render_frames_in_one_file(self, frameID: Sequence[int], **kwargs) -> str:
        raise NotImplementedError

    @abstractmethod
    def _render_frames(self, frameID: Sequence[int], **kwargs) -> list[str]:
        raise NotImplementedError
