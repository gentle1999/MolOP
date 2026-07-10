"""
Author: TMJ
Date: 2025-07-29 17:00:29
LastEditors: TMJ
LastEditTime: 2026-03-23 18:58:35
Description: 请填写简介
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from molop.io.base_models.ChemFile import BaseQMInputFile
from molop.io.base_models.Mixins import (
    DiskStorageMixin,
    FileMixin,
    MemoryStorageMixin,
)
from molop.io.logic.gaussian.input.frame_models.GJFFileFrame import (
    GJFFileFrameDisk,
    GJFFileFrameMemory,
)


if TYPE_CHECKING:
    from molop.io.codec_registry import Registry


class GJFFileMixin(FileMixin):
    file_frame_separator = "\n--Link1--\n"


class GJFFileMemory(MemoryStorageMixin, GJFFileMixin, BaseQMInputFile[GJFFileFrameMemory]): ...


class GJFFileDisk(DiskStorageMixin, GJFFileMixin, BaseQMInputFile[GJFFileFrameDisk]): ...


def register(registry: Registry) -> None:
    """Register this file model as a writer codec (renderer)."""

    from typing import cast

    from molop.io.codecs._shared.writer_helpers import (
        FileRendererWriter,
        FrameRendererWriter,
        StructureLevel,
        WriterCodec,
    )

    priority = 100

    @registry.writer_factory(
        format_id="gjf",
        required_level=StructureLevel.COORDS,
        domain="file",
        default_graph_policy="prefer",
        priority=priority,
    )
    def _factory() -> WriterCodec:
        return cast(
            WriterCodec,
            FileRendererWriter(
                format_id="gjf",
                required_level=StructureLevel.COORDS,
                file_cls=GJFFileDisk,
                frame_cls=GJFFileFrameDisk,
                priority=priority,
            ),
        )

    @registry.writer_factory(
        format_id="gjf",
        required_level=StructureLevel.COORDS,
        domain="frame",
        default_graph_policy="prefer",
        priority=priority,
    )
    def _frame_factory() -> WriterCodec:
        return cast(
            WriterCodec,
            FrameRendererWriter(
                format_id="gjf",
                required_level=StructureLevel.COORDS,
                frame_cls=GJFFileFrameDisk,
                priority=priority,
            ),
        )
