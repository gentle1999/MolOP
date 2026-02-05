"""
Author: TMJ
Date: 2025-12-14 23:26:19
LastEditors: TMJ
LastEditTime: 2026-02-05 19:53:16
Description: 请填写简介
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING, cast

from molop.io.base_models.ChemFile import BaseCoordsFile
from molop.io.base_models.Mixins import (
    DiskStorageMixin,
    FileMixin,
    MemoryStorageMixin,
    _HasRenderableFrames,
)
from molop.io.logic.coords_frame_models.SMIFileFrame import SMIFileFrameDisk, SMIFileFrameMemory


if TYPE_CHECKING:
    from molop.io.codec_registry import Registry


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


def register(registry: Registry) -> None:
    """Register this file model as a writer codec (renderer)."""

    from typing import cast

    from molop.io.codecs._shared.writer_helpers import (
        FileRendererWriter,
        StructureLevel,
        WriterCodec,
    )

    priority = 100

    @registry.writer_factory(
        format_id="smi",
        required_level=StructureLevel.GRAPH,
        priority=priority,
    )
    def _factory() -> WriterCodec:
        return cast(
            WriterCodec,
            FileRendererWriter(
                format_id="smi",
                required_level=StructureLevel.GRAPH,
                file_cls=SMIFileDisk,
                frame_cls=SMIFileFrameDisk,
                priority=priority,
            ),
        )
