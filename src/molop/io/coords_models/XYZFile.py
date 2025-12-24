"""
Author: TMJ
Date: 2025-07-29 16:54:42
LastEditors: TMJ
LastEditTime: 2025-12-14 16:00:19
Description: 请填写简介
"""

from collections.abc import Sequence
from typing import TYPE_CHECKING, Protocol

from molop.io.base_models.ChemFile import BaseCoordsFile
from molop.io.base_models.Mixins import DiskStorageMixin, FileMixin, MemoryStorageMixin
from molop.io.coords_models.XYZFileFrame import XYZFileFrameDisk, XYZFileFrameMemory


class XYZFileProtocol(Protocol):
    frames: list[XYZFileFrameDisk | XYZFileFrameMemory]


if TYPE_CHECKING:

    class _XYZFileProtocol(XYZFileProtocol, FileMixin): ...
else:

    class _XYZFileProtocol(FileMixin): ...


class XYZFileMixin(_XYZFileProtocol):
    def _render_frames_in_one_file(self, frameID: Sequence[int], **kwargs) -> str:
        return "\n".join(self._render_frames(frameID, **kwargs))

    def _render_frames(self, frameID: Sequence[int], **kwargs) -> list[str]:
        return [frame._render(**kwargs) for frame in self.frames if frame.frame_id in frameID]


class XYZFileMemory(MemoryStorageMixin, XYZFileMixin, BaseCoordsFile[XYZFileFrameMemory]): ...


class XYZFileDisk(DiskStorageMixin, XYZFileMixin, BaseCoordsFile[XYZFileFrameDisk]): ...
