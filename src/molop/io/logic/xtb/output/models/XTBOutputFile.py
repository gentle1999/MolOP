from __future__ import annotations

from typing import TYPE_CHECKING, ClassVar

from pint._typing import UnitLike

from molop.io.base_models.ChemFile import BaseCalcFile
from molop.io.base_models.Mixins import DiskStorageMixin, MemoryStorageMixin
from molop.io.logic.xtb.common import XTBOutputQMFieldsMixin
from molop.io.logic.xtb.output.frame_models.XTBOutputFileFrame import (
    XTBOutputFileFrameDisk,
    XTBOutputFileFrameMemory,
)
from molop.unit import atom_ureg


if TYPE_CHECKING:
    from molop.io.codec_registry import Registry


class XTBOutputFileMixin(XTBOutputQMFieldsMixin):
    default_units: ClassVar[dict[str, UnitLike]] = {
        **XTBOutputQMFieldsMixin.default_units,
        "running_time": atom_ureg.second,
        "temperature": atom_ureg.kelvin,
        "electron_temperature": atom_ureg.kelvin,
    }


class XTBOutputFileMemory(
    MemoryStorageMixin,
    XTBOutputFileMixin,
    BaseCalcFile[XTBOutputFileFrameMemory],
): ...


class XTBOutputFileDisk(
    DiskStorageMixin,
    XTBOutputFileMixin,
    BaseCalcFile[XTBOutputFileFrameDisk],
): ...


def register(registry: Registry) -> None:
    _ = registry


__all__ = ["XTBOutputFileDisk", "XTBOutputFileMemory"]
