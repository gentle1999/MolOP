from __future__ import annotations

from typing import TYPE_CHECKING

from molop.io.base_models.ChemFile import BaseQMInputFile
from molop.io.base_models.Mixins import (
    DiskStorageMixin,
    MemoryStorageMixin,
)
from molop.io.logic.qminput_frame_models.ORCAInpFileFrame import (
    ORCAInpFileFrameDisk,
    ORCAInpFileFrameMemory,
)


if TYPE_CHECKING:
    from molop.io.codec_registry import Registry


class ORCAInpFileMixin:
    """Structured ORCA input file model.

    ORCA rendering is intentionally not registered until a structured renderer exists.
    """


class ORCAInpFileMemory(
    MemoryStorageMixin, ORCAInpFileMixin, BaseQMInputFile[ORCAInpFileFrameMemory]
): ...


class ORCAInpFileDisk(
    DiskStorageMixin, ORCAInpFileMixin, BaseQMInputFile[ORCAInpFileFrameDisk]
): ...


def register(registry: Registry) -> None:
    _ = registry
