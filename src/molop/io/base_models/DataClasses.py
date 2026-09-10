"""Backward-compatible re-exports for structured calculation data classes.

Definitions live in domain-specific modules under :mod:`data_classes`; this
module remains as the stable import path for existing callers.
"""

from .data_classes.coordinates import (
    AtomInInternalCoords,  # noqa: F401
    CoordinateContainer,  # noqa: F401
    CoordinateParameter,  # noqa: F401
    CoordinateParameters,  # noqa: F401
    InternalCoords,  # noqa: F401
)
from .data_classes.electronic import (
    ElectronicConfiguration,  # noqa: F401
    ElectronicState,  # noqa: F401
    ElectronicStates,  # noqa: F401
    MultireferenceResult,  # noqa: F401
)
from .data_classes.orbitals import (  # noqa: F401
    MolecularOrbitals,
    MoleculeOrbital,
    NaturalAtomicOrbital,
    NaturalAtomicOrbitals,
    NaturalBondOrbital,
    NaturalBondOrbitals,
)
from .data_classes.populations import (  # noqa: F401
    AtomicPopulationSeries,
    ChargeSpinPopulations,
    TotalSpin,
)
from .data_classes.properties import (  # noqa: F401
    BondOrders,
    Dispersions,
    Polarizability,
    SinglePointProperties,
)
from .data_classes.qm_requests import (
    ActiveSpace,  # noqa: F401
    ExcitedStateRequest,  # noqa: F401
    ExplicitSolventRequest,  # noqa: F401
    MultireferenceRequest,  # noqa: F401
    MultireferenceStateBlock,  # noqa: F401
    QMBasisSet,  # noqa: F401
    QMModelChemistry,  # noqa: F401
    QMResourceRequest,  # noqa: F401
    QMTaskRequest,  # noqa: F401
)
from .data_classes.solvation import ImplicitSolvation  # noqa: F401
from .data_classes.spectroscopy import NMR, ShieldingTensor  # noqa: F401
from .data_classes.status import (  # noqa: F401
    GeometryOptimizationStatus,
    Status,
)
from .data_classes.thermochemistry import (  # noqa: F401
    Energies,
    EnergyObservation,
    ThermalInformations,
)
from .data_classes.vibrations import Vibration, Vibrations  # noqa: F401
