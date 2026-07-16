from __future__ import annotations

from molop.io.base_models.ChemFile import BaseCalcFile
from molop.io.base_models.Mixins import DiskStorageMixin, MemoryStorageMixin
from molop.io.logic.gaussian.fchk.output.frame_models.G16FchkFileFrame import (
    G16FchkFileFrameDisk,
    G16FchkFileFrameMemory,
)
from molop.io.logic.gaussian.input.GaussianRoute import GaussianRouteSemanticFieldsMixin


class G16FchkFileMixin(GaussianRouteSemanticFieldsMixin):
    """Gaussian formatted checkpoint file model."""


class G16FchkFileMemory(
    MemoryStorageMixin,
    G16FchkFileMixin,
    BaseCalcFile[G16FchkFileFrameMemory],
): ...


class G16FchkFileDisk(
    DiskStorageMixin,
    G16FchkFileMixin,
    BaseCalcFile[G16FchkFileFrameDisk],
): ...


G16FchkFileMemory.model_rebuild()
G16FchkFileDisk.model_rebuild()


__all__ = ["G16FchkFileDisk", "G16FchkFileMemory"]
