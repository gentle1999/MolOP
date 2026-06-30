from __future__ import annotations

from typing import TYPE_CHECKING, ClassVar, cast

from pint._typing import UnitLike
from pydantic import Field, model_validator
from typing_extensions import Self

from molop.io.base_models.ChemFile import BaseCalcFile
from molop.io.base_models.Mixins import DiskStorageMixin, MemoryStorageMixin
from molop.io.logic.QM_frame_models.ORCALogFileFrame import (
    ORCALogFileFrameDisk,
    ORCALogFileFrameMemory,
)
from molop.unit import atom_ureg


if TYPE_CHECKING:
    from molop.io.codec_registry import Registry


class ORCALogFileMixin:
    default_units: ClassVar[dict[str, UnitLike]] = {
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
        typed_self = cast(BaseCalcFile, self)
        typed_self.qm_software = "ORCA"
        typed_self.backfill_common_qm_containers_from_legacy()
        typed_self.project_common_qm_fields()
        return self


class ORCALogFileMemory(
    MemoryStorageMixin, ORCALogFileMixin, BaseCalcFile[ORCALogFileFrameMemory]
): ...


class ORCALogFileDisk(DiskStorageMixin, ORCALogFileMixin, BaseCalcFile[ORCALogFileFrameDisk]): ...


def register(registry: Registry) -> None:
    _ = registry
