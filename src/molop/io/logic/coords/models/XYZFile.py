"""
Author: TMJ
Date: 2025-07-29 16:54:42
LastEditors: TMJ
LastEditTime: 2026-02-04 11:42:13
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
from molop.io.logic.coords.frame_models.XYZFileFrame import XYZFileFrameDisk, XYZFileFrameMemory


if TYPE_CHECKING:
    from molop.io.codec_registry import Registry


class XYZFileMixin(FileMixin):
    file_frame_separator = "\n"


class XYZFileMemory(MemoryStorageMixin, XYZFileMixin, BaseCoordsFile[XYZFileFrameMemory]): ...


class XYZFileDisk(DiskStorageMixin, XYZFileMixin, BaseCoordsFile[XYZFileFrameDisk]): ...


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
        format_id="xyz",
        required_level=StructureLevel.COORDS,
        domain="file",
        default_graph_policy="coords",
        priority=priority,
    )
    def _factory() -> WriterCodec:
        return cast(
            WriterCodec,
            FileRendererWriter(
                format_id="xyz",
                required_level=StructureLevel.COORDS,
                file_cls=XYZFileDisk,
                frame_cls=XYZFileFrameDisk,
                priority=priority,
            ),
        )

    @registry.writer_factory(
        format_id="xyz",
        required_level=StructureLevel.COORDS,
        domain="frame",
        default_graph_policy="coords",
        priority=priority,
    )
    def _frame_factory() -> WriterCodec:
        return cast(
            WriterCodec,
            FrameRendererWriter(
                format_id="xyz",
                required_level=StructureLevel.COORDS,
                frame_cls=XYZFileFrameDisk,
                priority=priority,
            ),
        )
