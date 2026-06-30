from __future__ import annotations

from typing import ClassVar, cast

from pint._typing import UnitLike
from pydantic import Field, model_validator
from typing_extensions import Self

from molop.io.base_models.ChemFileFrame import BaseCalcFrame
from molop.io.base_models.Mixins import DiskStorageMixin, MemoryStorageMixin
from molop.unit import atom_ureg


class ORCALogFileFrameMixin:
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
    qm_software: str = Field(default="ORCA")
    input_file_name: str = Field(default="", description="ORCA input file name printed in output")
    auxiliary_basis_set: str = Field(
        default="",
        description="ORCA auxiliary basis keyword derived from the printed input",
    )
    dispersion_correction: str = Field(
        default="",
        description="ORCA dispersion correction keyword derived from the printed input",
    )

    @model_validator(mode="after")
    def _normalize_common_orca_fields(self) -> Self:
        typed_self = cast(BaseCalcFrame, self)
        typed_self.qm_software = "ORCA"
        typed_self.backfill_common_qm_containers_from_legacy()
        typed_self.project_common_qm_fields()
        return self


class ORCALogFileFrameMemory(
    MemoryStorageMixin, ORCALogFileFrameMixin, BaseCalcFrame["ORCALogFileFrameMemory"]
): ...


class ORCALogFileFrameDisk(
    DiskStorageMixin, ORCALogFileFrameMixin, BaseCalcFrame["ORCALogFileFrameDisk"]
): ...


ORCALogFileFrameMemory.model_rebuild()
ORCALogFileFrameDisk.model_rebuild()
