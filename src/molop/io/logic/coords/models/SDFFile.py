"""
Author: TMJ
Date: 2025-07-29 16:59:36
LastEditors: TMJ
LastEditTime: 2026-02-04 09:41:56
Description: 请填写简介
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from molop.io.base_models.ChemFile import BaseCoordsFile
from molop.io.base_models.Mixins import (
    DiskStorageMixin,
    FileMixin,
    MemoryStorageMixin,
)
from molop.io.logic.coords.frame_models.SDFFileFrame import SDFFileFrameDisk, SDFFileFrameMemory


if TYPE_CHECKING:
    from molop.io.codec_registry import Registry


class SDFFileMixin(FileMixin):
    file_frame_separator = "$$$$\n"


class SDFFileMemory(MemoryStorageMixin, SDFFileMixin, BaseCoordsFile[SDFFileFrameMemory]): ...


class SDFFileDisk(DiskStorageMixin, SDFFileMixin, BaseCoordsFile[SDFFileFrameDisk]): ...


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
        format_id="sdf",
        required_level=StructureLevel.GRAPH,
        domain="file",
        default_graph_policy="strict",
        priority=priority,
    )
    def _factory() -> WriterCodec:
        return cast(
            WriterCodec,
            FileRendererWriter(
                format_id="sdf",
                required_level=StructureLevel.GRAPH,
                file_cls=SDFFileDisk,
                frame_cls=SDFFileFrameDisk,
                priority=priority,
            ),
        )

    @registry.writer_factory(
        format_id="sdf",
        required_level=StructureLevel.GRAPH,
        domain="frame",
        default_graph_policy="strict",
        priority=priority,
    )
    def _frame_factory() -> WriterCodec:
        return cast(
            WriterCodec,
            FrameRendererWriter(
                format_id="sdf",
                required_level=StructureLevel.GRAPH,
                frame_cls=SDFFileFrameDisk,
                priority=priority,
            ),
        )
