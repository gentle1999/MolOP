from __future__ import annotations

from typing import Any, ClassVar, cast

from pint._typing import UnitLike
from pydantic import model_validator
from typing_extensions import Self

from molop.io.base_models.ChemFileFrame import BaseCalcFrame
from molop.io.base_models.Mixins import DiskStorageMixin, MemoryStorageMixin
from molop.io.logic.xtb.common import XTBOutputQMFieldsMixin
from molop.unit import atom_ureg


class XTBOutputFileFrameMixin(XTBOutputQMFieldsMixin):
    default_units: ClassVar[dict[str, UnitLike]] = {
        **XTBOutputQMFieldsMixin.default_units,
        "coords": atom_ureg.angstrom,
        "running_time": atom_ureg.second,
        "rotation_constants": atom_ureg.gigahertz,
        "temperature": atom_ureg.kelvin,
        "electron_temperature": atom_ureg.kelvin,
    }

    @model_validator(mode="after")
    def physical_check(self) -> Self:
        """Allow property-only xTB logs whose command-line geometry is not embedded."""
        typed_self = cast(Any, self)
        num_atoms = len(typed_self.atoms)
        if (
            num_atoms
            and typed_self.vibrations
            and len(typed_self.vibrations)
            not in (
                num_atoms * 3 - 6,
                num_atoms * 3 - 5,
                num_atoms * 3 - 3,
            )
        ):
            raise ValueError(
                f"Invalid vibrational mode count: {len(typed_self.vibrations)} "
                f"for {num_atoms} atoms"
            )
        return self


class XTBOutputFileFrameMemory(
    MemoryStorageMixin,
    XTBOutputFileFrameMixin,
    BaseCalcFrame["XTBOutputFileFrameMemory"],
): ...


class XTBOutputFileFrameDisk(
    DiskStorageMixin,
    XTBOutputFileFrameMixin,
    BaseCalcFrame["XTBOutputFileFrameDisk"],
): ...


XTBOutputFileFrameMemory.model_rebuild()
XTBOutputFileFrameDisk.model_rebuild()


__all__ = ["XTBOutputFileFrameDisk", "XTBOutputFileFrameMemory"]
