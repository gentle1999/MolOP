"""
Author: TMJ
Date: 2025-07-29 16:59:36
LastEditors: TMJ
LastEditTime: 2025-12-14 16:00:41
Description: 请填写简介
"""

from typing import TYPE_CHECKING, Protocol, Sequence

from molop.io.base_models.ChemFile import BaseCoordsFile
from molop.io.base_models.Mixins import DiskStorageMixin, FileMixin, MemoryStorageMixin
from molop.io.coords_models.SDFFileFrame import SDFFileFrameDisk, SDFFileFrameMemory


class SDFFileProtocol(Protocol):
    frames: list[SDFFileFrameDisk | SDFFileFrameMemory]


if TYPE_CHECKING:

    class _SDFFileProtocol(SDFFileProtocol, FileMixin): ...
else:

    class _SDFFileProtocol(FileMixin): ...


class SDFFileMixin(_SDFFileProtocol):
    def _render_frames_in_one_file(self, frameID: Sequence[int], **kwargs) -> str:
        return "$$$$\n".join(self._render_frames(frameID, **kwargs))

    def _render_frames(self, frameID: Sequence[int], **kwargs) -> list[str]:
        return [
            frame._render(**kwargs)
            for frame in self.frames
            if frame.frame_id in frameID
        ]


class SDFFileMemory(
    MemoryStorageMixin, SDFFileMixin, BaseCoordsFile[SDFFileFrameMemory]
): ...


class SDFFileDisk(DiskStorageMixin, SDFFileMixin, BaseCoordsFile[SDFFileFrameDisk]): ...
