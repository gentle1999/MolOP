from molop.io.base_models.ChemFileFrame import BaseCalcFrame
from molop.io.base_models.Mixins import DiskStorageMixin, MemoryStorageMixin
from molop.io.logic.gaussian.input.GaussianRoute import GaussianRouteSemanticFieldsMixin


class G16FchkFileFrameMixin(GaussianRouteSemanticFieldsMixin):
    """Structured properties from one Gaussian formatted checkpoint."""


class G16FchkFileFrameMemory(
    MemoryStorageMixin,
    G16FchkFileFrameMixin,
    BaseCalcFrame["G16FchkFileFrameMemory"],
): ...


class G16FchkFileFrameDisk(
    DiskStorageMixin,
    G16FchkFileFrameMixin,
    BaseCalcFrame["G16FchkFileFrameDisk"],
): ...


G16FchkFileFrameMemory.model_rebuild()
G16FchkFileFrameDisk.model_rebuild()


__all__ = ["G16FchkFileFrameDisk", "G16FchkFileFrameMemory"]
