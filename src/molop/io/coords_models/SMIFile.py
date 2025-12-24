"""
Author: TMJ
Date: 2025-12-14 23:26:19
LastEditors: TMJ
LastEditTime: 2025-12-14 23:26:59
Description: 请填写简介
"""

from collections.abc import Sequence
from typing import TYPE_CHECKING, Protocol

from molop.io.base_models.ChemFile import BaseCoordsFile
from molop.io.base_models.Mixins import DiskStorageMixin, FileMixin, MemoryStorageMixin
from molop.io.coords_models.SMIFileFrame import SMIFileFrameDisk, SMIFileFrameMemory


class SMIFileProtocol(Protocol):
    frames: list[SMIFileFrameDisk | SMIFileFrameMemory]


if TYPE_CHECKING:

    class _SMIFileProtocol(SMIFileProtocol, FileMixin): ...
else:

    class _SMIFileProtocol(FileMixin): ...


class SMIFileMixin(_SMIFileProtocol):
    def _render_frames_in_one_file(self, frameID: Sequence[int], **kwargs) -> str:
        return "\n".join(self._render_frames(frameID, **kwargs))

    def _render_frames(self, frameID: Sequence[int], **kwargs) -> list[str]:
        return [frame._render(**kwargs) for frame in self.frames if frame.frame_id in frameID]


class SMIFileMemory(MemoryStorageMixin, SMIFileMixin, BaseCoordsFile[SMIFileFrameMemory]): ...


class SMIFileDisk(DiskStorageMixin, SMIFileMixin, BaseCoordsFile[SMIFileFrameDisk]): ...
