"""
Author: TMJ
Date: 2026-02-10 00:00:00
LastEditors: TMJ
LastEditTime: 2026-02-10 00:00:00
Description: ORCA input file models
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING, cast

from molop.io.base_models.ChemFile import BaseQMInputFile
from molop.io.base_models.Mixins import (
    DiskStorageMixin,
    FileMixin,
    MemoryStorageMixin,
    _HasRenderableFrames,
)
from molop.io.logic.qminput_frame_models.ORCAInpFileFrame import (
    ORCAInpFileFrameDisk,
    ORCAInpFileFrameMemory,
)


if TYPE_CHECKING:
    from molop.io.codec_registry import Registry


class ORCAInpFileMixin(FileMixin):
    def _render_frames_in_one_file(self, frameID: Sequence[int], **kwargs) -> str:
        typed_self = cast(_HasRenderableFrames, self)
        rendered_frames = [
            frame._render(**kwargs) for frame in typed_self.frames if frame.frame_id in frameID
        ]
        return "\n$new_job\n".join(rendered_frames)

    def _render_frames(self, frameID: Sequence[int], **kwargs) -> list[str]:
        typed_self = cast(_HasRenderableFrames, self)
        return [frame._render(**kwargs) for frame in typed_self.frames if frame.frame_id in frameID]


class ORCAInpFileMemory(
    MemoryStorageMixin, ORCAInpFileMixin, BaseQMInputFile[ORCAInpFileFrameMemory]
): ...


class ORCAInpFileDisk(
    DiskStorageMixin, ORCAInpFileMixin, BaseQMInputFile[ORCAInpFileFrameDisk]
): ...


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
        format_id="orcainp",
        required_level=StructureLevel.COORDS,
        priority=priority,
    )
    def _factory() -> WriterCodec:
        return cast(
            WriterCodec,
            FileRendererWriter(
                format_id="orcainp",
                required_level=StructureLevel.COORDS,
                file_cls=ORCAInpFileDisk,
                frame_cls=ORCAInpFileFrameDisk,
                priority=priority,
            ),
        )
