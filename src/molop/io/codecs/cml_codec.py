"""CML writer codec.

This codec is not backed by a dedicated `coords_models` module; it uses the CML
file/frame adapters defined in `molop.io.codecs._shared.writer_helpers`.
"""

from __future__ import annotations

from typing import cast

from molop.io.codec_registry import Registry
from molop.io.codecs._shared.writer_helpers import (
    CMLFileDisk,
    CMLFileFrameDisk,
    FileRendererWriter,
    StructureLevel,
    WriterCodec,
)


def register(registry: Registry) -> None:
    priority = 100

    @registry.writer_factory(
        format_id="cml",
        required_level=StructureLevel.GRAPH,
        priority=priority,
    )
    def _factory() -> WriterCodec:
        return cast(
            WriterCodec,
            FileRendererWriter(
                format_id="cml",
                required_level=StructureLevel.GRAPH,
                file_cls=CMLFileDisk,
                frame_cls=CMLFileFrameDisk,
                priority=priority,
            ),
        )


__all__ = ["register"]
