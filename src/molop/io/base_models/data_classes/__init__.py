"""Domain-specific structured calculation data classes."""

from .coordinates import (
    AtomInInternalCoords,
    CoordinateContainer,
    CoordinateParameter,
    CoordinateParameters,
    InternalCoords,
)
from .electronic import (
    ElectronicConfiguration,
    ElectronicState,
    ElectronicStates,
    MultireferenceResult,
)
from .orbitals import (
    MolecularOrbitals,
    MoleculeOrbital,
    NaturalAtomicOrbital,
    NaturalAtomicOrbitals,
    NaturalBondOrbital,
    NaturalBondOrbitals,
)
from .populations import AtomicPopulationSeries, ChargeSpinPopulations, TotalSpin
from .properties import BondOrders, Dispersions, Polarizability, SinglePointProperties
from .qm_requests import (
    ActiveSpace,
    ExcitedStateRequest,
    ExplicitSolventRequest,
    MultireferenceRequest,
    MultireferenceStateBlock,
    QMBasisSet,
    QMModelChemistry,
    QMResourceRequest,
    QMTaskRequest,
)
from .solvation import ImplicitSolvation
from .spectroscopy import NMR, ShieldingTensor
from .status import GeometryOptimizationStatus, Status
from .thermochemistry import Energies, EnergyObservation, ThermalInformations
from .vibrations import Vibration, Vibrations


__all__ = [
    "ActiveSpace",
    "AtomInInternalCoords",
    "AtomicPopulationSeries",
    "BondOrders",
    "ChargeSpinPopulations",
    "Dispersions",
    "CoordinateContainer",
    "CoordinateParameter",
    "CoordinateParameters",
    "ElectronicConfiguration",
    "ElectronicState",
    "ElectronicStates",
    "Energies",
    "EnergyObservation",
    "ExplicitSolventRequest",
    "ExcitedStateRequest",
    "GeometryOptimizationStatus",
    "ImplicitSolvation",
    "InternalCoords",
    "MultireferenceRequest",
    "MultireferenceResult",
    "MultireferenceStateBlock",
    "MolecularOrbitals",
    "MoleculeOrbital",
    "NMR",
    "NaturalAtomicOrbital",
    "NaturalAtomicOrbitals",
    "NaturalBondOrbital",
    "NaturalBondOrbitals",
    "Polarizability",
    "QMBasisSet",
    "QMModelChemistry",
    "QMResourceRequest",
    "QMTaskRequest",
    "SinglePointProperties",
    "Status",
    "ShieldingTensor",
    "ThermalInformations",
    "TotalSpin",
    "Vibration",
    "Vibrations",
]
