"""
Author: TMJ
Date: 2025-07-29 12:36:28
LastEditors: TMJ
LastEditTime: 2026-02-04 16:04:37
Description: 请填写简介
"""

import os
from collections.abc import Sequence
from typing import Any, Literal, Protocol, cast, overload

from pydantic import Field, computed_field
from typing_extensions import Self


class _RenderableFrame(Protocol):
    frame_id: int

    def _render(self, **kwargs) -> str: ...


class _HasRenderableFrames(Protocol):
    frames: Sequence[_RenderableFrame]

    def _render_frames_in_one_file(self, frameID: Sequence[int], **kwargs) -> str: ...

    def _render_frames(self, frameID: Sequence[int], **kwargs) -> list[str]: ...


class MemoryStorageMixin: ...


class DiskStorageMixin:
    file_path: str = Field(default="", repr=False)

    @computed_field()  # type: ignore[prop-decorator] # type: ignore[prop-decorator]
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

    def __eq__(self, other: Any) -> bool:
        assert isinstance(other, DiskStorageMixin), (
            "other should be an instance of DiskStorageMixin"
        )
        return self.file_path == other.file_path

    def to_summary_dict(self, brief: bool = True, **kwargs) -> dict[tuple[str, str], Any]:
        return {
            ("DiskStorage", "FilePath"): self.file_path,
            ("DiskStorage", "FileFormat"): self.file_format,
            **super().to_summary_dict(brief=brief, **kwargs),  # type: ignore
        }


class FileMixin:
    @overload
    def _render(
        self,
        frameID: Sequence[int] | int | Literal["all"],
        embed_in_one_file: Literal[True],
        **kwargs,
    ) -> str: ...
    @overload
    def _render(
        self,
        frameID: Sequence[int] | int | Literal["all"],
        embed_in_one_file: Literal[False],
        **kwargs,
    ) -> list[str]: ...
    def _render(
        self,
        frameID: Sequence[int] | int | Literal["all"] = -1,
        embed_in_one_file: bool = True,
        **kwargs,
    ) -> str | list[str]:
        typed_self = cast(_HasRenderableFrames, self)
        if isinstance(frameID, int):
            frameIDs = [frameID if frameID >= 0 else len(typed_self.frames) + frameID]
        elif frameID == "all":
            frameIDs = list(range(len(typed_self.frames)))
        elif isinstance(frameID, Sequence):
            frameIDs = []
            for i in frameID:
                assert isinstance(i, int), "frameID should be a sequence of integers"
                frameIDs.append(i)
        else:
            raise ValueError("frameID should be an integer, a sequence of integers, or 'all'")
        if embed_in_one_file:
            return typed_self._render_frames_in_one_file(frameIDs, **kwargs)
        else:
            return typed_self._render_frames(frameIDs, **kwargs)

    def _render_frames_in_one_file(self, frameID: Sequence[int], **kwargs) -> str:
        raise NotImplementedError(
            f"{self.__class__.__name__} has not implemented _render_frames_in_one_file method yet"
        )

    def _render_frames(self, frameID: Sequence[int], **kwargs) -> list[str]:
        raise NotImplementedError(
            f"{self.__class__.__name__} has not implemented _render_frames method yet"
        )
