"""
Author: TMJ
Date: 2025-12-14 23:26:19
LastEditors: TMJ
LastEditTime: 2026-02-04 11:40:05
Description: 请填写简介
"""

from collections.abc import Sequence
from typing import cast

from molop.io.base_models.ChemFile import BaseCoordsFile
from molop.io.base_models.Mixins import (
    DiskStorageMixin,
    FileMixin,
    MemoryStorageMixin,
    _HasRenderableFrames,
)
from molop.io.coords_models.SMIFileFrame import SMIFileFrameDisk, SMIFileFrameMemory


class SMIFileMixin(FileMixin):
    def _render_frames_in_one_file(self, frameID: Sequence[int], **kwargs) -> str:
        typed_self = cast(_HasRenderableFrames, self)
        return "\n".join(
            frame._render(**kwargs) for frame in typed_self.frames if frame.frame_id in frameID
        )

    def _render_frames(self, frameID: Sequence[int], **kwargs) -> list[str]:
        typed_self = cast(_HasRenderableFrames, self)
        return [frame._render(**kwargs) for frame in typed_self.frames if frame.frame_id in frameID]


class SMIFileMemory(MemoryStorageMixin, SMIFileMixin, BaseCoordsFile[SMIFileFrameMemory]): ...


class SMIFileDisk(DiskStorageMixin, SMIFileMixin, BaseCoordsFile[SMIFileFrameDisk]): ...
