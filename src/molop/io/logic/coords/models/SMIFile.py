"""
Author: TMJ
Date: 2025-12-14 23:26:19
LastEditors: TMJ
LastEditTime: 2026-02-05 19:53:16
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
from molop.io.logic.coords.frame_models.SMIFileFrame import SMIFileFrameDisk, SMIFileFrameMemory


if TYPE_CHECKING:
    from molop.io.codec_registry import Registry


class SMIFileMixin(FileMixin):
    file_frame_separator = "\n"


class SMIFileMemory(MemoryStorageMixin, SMIFileMixin, BaseCoordsFile[SMIFileFrameMemory]): ...


class SMIFileDisk(DiskStorageMixin, SMIFileMixin, BaseCoordsFile[SMIFileFrameDisk]): ...


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
        format_id="smi",
        required_level=StructureLevel.GRAPH,
        domain="file",
        default_graph_policy="strict",
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

    @registry.writer_factory(
        format_id="smi",
        required_level=StructureLevel.GRAPH,
        domain="frame",
        default_graph_policy="strict",
        priority=priority,
    )
    def _frame_factory() -> WriterCodec:
        return cast(
            WriterCodec,
            FrameRendererWriter(
                format_id="smi",
                required_level=StructureLevel.GRAPH,
                frame_cls=SMIFileFrameDisk,
                priority=priority,
            ),
        )
