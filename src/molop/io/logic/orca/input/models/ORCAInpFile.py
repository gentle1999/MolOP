from __future__ import annotations

from typing import TYPE_CHECKING

from molop.io.base_models.ChemFile import BaseQMInputFile
from molop.io.base_models.Mixins import (
    DiskStorageMixin,
    FileMixin,
    MemoryStorageMixin,
)
from molop.io.logic.orca.input.frame_models.ORCAInpFileFrame import (
    ORCAInpFileFrameDisk,
    ORCAInpFileFrameMemory,
)


if TYPE_CHECKING:
    from molop.io.codec_registry import Registry


class ORCAInpFileMixin(FileMixin):
    """Structured and renderable ORCA input file model."""

    file_frame_separator = "\n$new_job\n"


class ORCAInpFileMemory(
    MemoryStorageMixin, ORCAInpFileMixin, BaseQMInputFile[ORCAInpFileFrameMemory]
): ...


class ORCAInpFileDisk(
    DiskStorageMixin, ORCAInpFileMixin, BaseQMInputFile[ORCAInpFileFrameDisk]
): ...


def register(registry: Registry) -> None:
    """Register ORCA input file and frame writer codecs."""

    from typing import cast

    from molop.io.codecs._shared.writer_helpers import (
        FileRendererWriter,
        FrameRendererWriter,
        StructureLevel,
        WriterCodec,
    )

    priority = 100

    @registry.writer_factory(
        format_id="orcainp",
        required_level=StructureLevel.COORDS,
        domain="file",
        default_graph_policy="coords",
        priority=priority,
    )
    def _factory() -> WriterCodec:
        return cast(
            WriterCodec,
            FileRendererWriter(
                format_id="orcainp",
                output_extension="inp",
                required_level=StructureLevel.COORDS,
                file_cls=ORCAInpFileDisk,
                frame_cls=ORCAInpFileFrameDisk,
                priority=priority,
            ),
        )

    @registry.writer_factory(
        format_id="orcainp",
        required_level=StructureLevel.COORDS,
        domain="frame",
        default_graph_policy="coords",
        priority=priority,
    )
    def _frame_factory() -> WriterCodec:
        return cast(
            WriterCodec,
            FrameRendererWriter(
                format_id="orcainp",
                output_extension="inp",
                required_level=StructureLevel.COORDS,
                frame_cls=ORCAInpFileFrameDisk,
                priority=priority,
            ),
        )
