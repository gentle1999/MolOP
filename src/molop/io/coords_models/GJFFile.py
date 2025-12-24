"""
Author: TMJ
Date: 2025-07-29 17:00:29
LastEditors: TMJ
LastEditTime: 2025-12-14 21:28:19
Description: 请填写简介
"""

from typing import TYPE_CHECKING, Protocol

from molop.io.base_models.ChemFile import BaseCoordsFile
from molop.io.base_models.Mixins import DiskStorageMixin, FileMixin, MemoryStorageMixin
from molop.io.coords_models.GJFFileFrame import GJFFileFrameDisk, GJFFileFrameMemory


class GJFFileProtocol(Protocol):
    frames: list[GJFFileFrameDisk | GJFFileFrameMemory]


if TYPE_CHECKING:

    class _GJFFileProtocol(GJFFileProtocol, FileMixin): ...
else:

    class _GJFFileProtocol(FileMixin): ...


class GJFFileMixin(_GJFFileProtocol):
    def _render_frames_in_one_file(self, frameID: list[int], **kwargs) -> str:
        return "\n--Link1--\n".join(self._render_frames(frameID, **kwargs))

    def _render_frames(self, frameID: list[int], **kwargs) -> list[str]:
        return [frame._render(**kwargs) for frame in self.frames if frame.frame_id in frameID]


class GJFFileMemory(MemoryStorageMixin, GJFFileMixin, BaseCoordsFile[GJFFileFrameMemory]): ...


class GJFFileDisk(DiskStorageMixin, GJFFileMixin, BaseCoordsFile[GJFFileFrameDisk]): ...
