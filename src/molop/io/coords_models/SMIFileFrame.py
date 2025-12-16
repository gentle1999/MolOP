from typing import TYPE_CHECKING, Protocol

from molop.io.base_models.Bases import BaseDataClassWithUnit
from molop.io.base_models.ChemFileFrame import BaseCoordsFrame
from molop.io.base_models.Mixins import DiskStorageMixin, MemoryStorageMixin
from molop.utils.types import OMol, RdMol


class SMIFileFrameProtocol(Protocol):
    rdmol: RdMol | None
    omol: OMol | None

    def to_canonical_SMILES(self) -> str: ...


if TYPE_CHECKING:

    class _SMIFileFrameProtocol(SMIFileFrameProtocol, BaseDataClassWithUnit): ...
else:

    class _SMIFileFrameProtocol(BaseDataClassWithUnit): ...


class SMIFileFrameMixin(_SMIFileFrameProtocol):
    def _render(self) -> str:
        return self.to_canonical_SMILES()


class SMIFileFrameMemory(
    MemoryStorageMixin, SMIFileFrameMixin, BaseCoordsFrame["SMIFileFrameMemory"]
): ...


class SMIFileFrameDisk(
    DiskStorageMixin, SMIFileFrameMixin, BaseCoordsFrame["SMIFileFrameDisk"]
): ...
