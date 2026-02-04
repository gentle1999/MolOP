"""
Author: TMJ
Date: 2025-07-29 16:59:36
LastEditors: TMJ
LastEditTime: 2026-02-04 09:41:56
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
from molop.io.coords_models.SDFFileFrame import SDFFileFrameDisk, SDFFileFrameMemory


class SDFFileMixin(FileMixin):
    def _render_frames_in_one_file(self, frameID: Sequence[int], **kwargs) -> str:
        typed_self = cast(_HasRenderableFrames, self)
        return "$$$$\n".join(
            frame._render(**kwargs) for frame in typed_self.frames if frame.frame_id in frameID
        )

    def _render_frames(self, frameID: Sequence[int], **kwargs) -> list[str]:
        typed_self = cast(_HasRenderableFrames, self)
        return [frame._render(**kwargs) for frame in typed_self.frames if frame.frame_id in frameID]


class SDFFileMemory(MemoryStorageMixin, SDFFileMixin, BaseCoordsFile[SDFFileFrameMemory]): ...


class SDFFileDisk(DiskStorageMixin, SDFFileMixin, BaseCoordsFile[SDFFileFrameDisk]): ...
