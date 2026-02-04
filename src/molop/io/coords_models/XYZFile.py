"""
Author: TMJ
Date: 2025-07-29 16:54:42
LastEditors: TMJ
LastEditTime: 2026-02-04 11:42:13
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
from molop.io.coords_models.XYZFileFrame import XYZFileFrameDisk, XYZFileFrameMemory


class XYZFileMixin(FileMixin):
    def _render_frames_in_one_file(self, frameID: Sequence[int], **kwargs) -> str:
        typed_self = cast(_HasRenderableFrames, self)
        rendered_frames = [
            frame._render(**kwargs) for frame in typed_self.frames if frame.frame_id in frameID
        ]
        return "\n".join(rendered_frames)

    def _render_frames(self, frameID: Sequence[int], **kwargs) -> list[str]:
        typed_self = cast(_HasRenderableFrames, self)
        return [frame._render(**kwargs) for frame in typed_self.frames if frame.frame_id in frameID]


class XYZFileMemory(MemoryStorageMixin, XYZFileMixin, BaseCoordsFile[XYZFileFrameMemory]): ...


class XYZFileDisk(DiskStorageMixin, XYZFileMixin, BaseCoordsFile[XYZFileFrameDisk]): ...
