from __future__ import annotations

from typing import ClassVar

from pint._typing import UnitLike

from molop.io.base_models.ChemFileFrame import BaseCalcFrame
from molop.io.base_models.Mixins import DiskStorageMixin, MemoryStorageMixin
from molop.io.logic.orca.common import ORCAOutputQMFieldsMixin
from molop.unit import atom_ureg


class ORCALogFileFrameMixin(ORCAOutputQMFieldsMixin):
    """Structured ORCA output frame.

    The model intentionally stores ORCA output data in the shared QM output
    containers from ``BaseCalcFrame``. ORCA-specific raw text remains available
    through ``frame_content``.
    """

    default_units: ClassVar[dict[str, UnitLike]] = {
        "coords": atom_ureg.angstrom,
        "forces": atom_ureg.Unit("hartree / bohr"),
        "running_time": atom_ureg.second,
    }


class ORCALogFileFrameMemory(
    MemoryStorageMixin, ORCALogFileFrameMixin, BaseCalcFrame["ORCALogFileFrameMemory"]
): ...


class ORCALogFileFrameDisk(
    DiskStorageMixin, ORCALogFileFrameMixin, BaseCalcFrame["ORCALogFileFrameDisk"]
): ...


ORCALogFileFrameMemory.model_rebuild()
ORCALogFileFrameDisk.model_rebuild()
