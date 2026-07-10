from __future__ import annotations

from typing import TYPE_CHECKING, ClassVar

from pint._typing import UnitLike

from molop.io.base_models.ChemFile import BaseCalcFile
from molop.io.base_models.Mixins import DiskStorageMixin, MemoryStorageMixin
from molop.io.logic.orca.common import ORCAOutputQMFieldsMixin
from molop.io.logic.orca.log.frame_models.ORCALogFileFrame import (
    ORCALogFileFrameDisk,
    ORCALogFileFrameMemory,
)
from molop.unit import atom_ureg


if TYPE_CHECKING:
    from molop.io.codec_registry import Registry


class ORCALogFileMixin(ORCAOutputQMFieldsMixin):
    default_units: ClassVar[dict[str, UnitLike]] = {
        "running_time": atom_ureg.second,
    }


class ORCALogFileMemory(
    MemoryStorageMixin, ORCALogFileMixin, BaseCalcFile[ORCALogFileFrameMemory]
): ...


class ORCALogFileDisk(DiskStorageMixin, ORCALogFileMixin, BaseCalcFile[ORCALogFileFrameDisk]): ...


def register(registry: Registry) -> None:
    _ = registry
