from __future__ import annotations

import os
import pickle
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from molop import BatchParseError, iter_parse_outcomes
from molop.io.base_models.ChemFileFrame import BaseQMInputFrame
from molop.io.base_models.DataClasses import QMModelChemistry
from molop.io.base_models.Molecule import Molecule
from molop.unit import atom_ureg


ROOT = Path(__file__).resolve().parents[1]
XYZ_FIXTURE = ROOT / "tests" / "test_files" / "xyz" / "dsgdb9nsd_125600-5" / "0.xyz"


def test_molecule_structure_edits_invalidate_cached_topology() -> None:
    molecule = Molecule.model_validate(
        {
            "atoms": [6, 1],
            "coords": np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]) * atom_ureg.angstrom,
            "bonds": [(0, 1, 1, 0)],
            "formal_charges": [0, 0],
            "formal_num_radicals": [0, 0],
        }
    )
    cached = molecule.rdmol
    assert cached is not None

    cached.GetAtomWithIdx(0).SetAtomicNum(8)
    assert molecule.rdmol is not None
    assert molecule.rdmol.GetAtomWithIdx(0).GetAtomicNum() == 6

    molecule.atoms[0] = 7
    assert molecule._rdmol is None
    assert molecule.rdmol is not None
    assert molecule.rdmol.GetAtomWithIdx(0).GetAtomicNum() == 7


def test_pickled_molecule_restores_tracking_state() -> None:
    molecule = Molecule.model_validate(
        {
            "atoms": [2],
            "coords": np.array([[0.0, 0.0, 0.0]]) * atom_ureg.angstrom,
        }
    )

    restored = pickle.loads(pickle.dumps(molecule))

    assert restored.rdmol is not None
    assert restored._structure_tracking_suspended == 0
    restored.atoms[0] = 1
    assert restored._rdmol is None


def test_vibration_domain_split_preserves_legacy_reexports() -> None:
    from molop.io.base_models.data_classes.vibrations import (
        Vibration as DomainVibration,
    )
    from molop.io.base_models.data_classes.vibrations import (
        Vibrations as DomainVibrations,
    )
    from molop.io.base_models.DataClasses import Vibration, Vibrations

    assert Vibration is DomainVibration
    assert Vibrations is DomainVibrations


def test_coordinate_domain_split_preserves_legacy_reexports() -> None:
    from molop.io.base_models.data_classes.coordinates import (
        AtomInInternalCoords as DomainAtomInInternalCoords,
    )
    from molop.io.base_models.data_classes.coordinates import (
        CoordinateContainer as DomainCoordinateContainer,
    )
    from molop.io.base_models.data_classes.coordinates import (
        CoordinateParameter as DomainCoordinateParameter,
    )
    from molop.io.base_models.data_classes.coordinates import (
        CoordinateParameters as DomainCoordinateParameters,
    )
    from molop.io.base_models.data_classes.coordinates import (
        InternalCoords as DomainInternalCoords,
    )
    from molop.io.base_models.DataClasses import (
        AtomInInternalCoords,
        CoordinateContainer,
        CoordinateParameter,
        CoordinateParameters,
        InternalCoords,
    )

    assert AtomInInternalCoords is DomainAtomInInternalCoords
    assert CoordinateContainer is DomainCoordinateContainer
    assert CoordinateParameter is DomainCoordinateParameter
    assert CoordinateParameters is DomainCoordinateParameters
    assert InternalCoords is DomainInternalCoords


def test_qm_request_domain_split_preserves_legacy_reexports() -> None:
    from molop.io.base_models.data_classes.qm_requests import (
        ActiveSpace as DomainActiveSpace,
    )
    from molop.io.base_models.data_classes.qm_requests import (
        ExcitedStateRequest as DomainExcitedStateRequest,
    )
    from molop.io.base_models.data_classes.qm_requests import (
        ExplicitSolventRequest as DomainExplicitSolventRequest,
    )
    from molop.io.base_models.data_classes.qm_requests import (
        MultireferenceRequest as DomainMultireferenceRequest,
    )
    from molop.io.base_models.data_classes.qm_requests import (
        MultireferenceStateBlock as DomainMultireferenceStateBlock,
    )
    from molop.io.base_models.data_classes.qm_requests import (
        QMBasisSet as DomainQMBasisSet,
    )
    from molop.io.base_models.data_classes.qm_requests import (
        QMModelChemistry as DomainQMModelChemistry,
    )
    from molop.io.base_models.data_classes.qm_requests import (
        QMResourceRequest as DomainQMResourceRequest,
    )
    from molop.io.base_models.data_classes.qm_requests import (
        QMTaskRequest as DomainQMTaskRequest,
    )
    from molop.io.base_models.DataClasses import (
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

    assert ActiveSpace is DomainActiveSpace
    assert ExplicitSolventRequest is DomainExplicitSolventRequest
    assert ExcitedStateRequest is DomainExcitedStateRequest
    assert MultireferenceRequest is DomainMultireferenceRequest
    assert MultireferenceStateBlock is DomainMultireferenceStateBlock
    assert QMBasisSet is DomainQMBasisSet
    assert QMModelChemistry is DomainQMModelChemistry
    assert QMResourceRequest is DomainQMResourceRequest
    assert QMTaskRequest is DomainQMTaskRequest


def test_electronic_domain_split_preserves_legacy_reexports() -> None:
    from molop.io.base_models.data_classes.electronic import (
        ElectronicConfiguration as DomainElectronicConfiguration,
    )
    from molop.io.base_models.data_classes.electronic import (
        ElectronicState as DomainElectronicState,
    )
    from molop.io.base_models.data_classes.electronic import (
        ElectronicStates as DomainElectronicStates,
    )
    from molop.io.base_models.data_classes.electronic import (
        MultireferenceResult as DomainMultireferenceResult,
    )
    from molop.io.base_models.DataClasses import (
        ElectronicConfiguration,
        ElectronicState,
        ElectronicStates,
        MultireferenceResult,
    )

    assert ElectronicConfiguration is DomainElectronicConfiguration
    assert ElectronicState is DomainElectronicState
    assert ElectronicStates is DomainElectronicStates
    assert MultireferenceResult is DomainMultireferenceResult


def test_thermochemistry_domain_split_preserves_legacy_reexports() -> None:
    from molop.io.base_models.data_classes.thermochemistry import (
        Energies as DomainEnergies,
    )
    from molop.io.base_models.data_classes.thermochemistry import (
        EnergyObservation as DomainEnergyObservation,
    )
    from molop.io.base_models.data_classes.thermochemistry import (
        ThermalInformations as DomainThermalInformations,
    )
    from molop.io.base_models.DataClasses import Energies, EnergyObservation, ThermalInformations

    assert Energies is DomainEnergies
    assert EnergyObservation is DomainEnergyObservation
    assert ThermalInformations is DomainThermalInformations


def test_orbital_domain_split_preserves_legacy_reexports() -> None:
    from molop.io.base_models.data_classes.orbitals import (
        MolecularOrbitals as DomainMolecularOrbitals,
    )
    from molop.io.base_models.data_classes.orbitals import (
        MoleculeOrbital as DomainMoleculeOrbital,
    )
    from molop.io.base_models.data_classes.orbitals import (
        NaturalAtomicOrbital as DomainNaturalAtomicOrbital,
    )
    from molop.io.base_models.data_classes.orbitals import (
        NaturalAtomicOrbitals as DomainNaturalAtomicOrbitals,
    )
    from molop.io.base_models.data_classes.orbitals import (
        NaturalBondOrbital as DomainNaturalBondOrbital,
    )
    from molop.io.base_models.data_classes.orbitals import (
        NaturalBondOrbitals as DomainNaturalBondOrbitals,
    )
    from molop.io.base_models.DataClasses import (
        MolecularOrbitals,
        MoleculeOrbital,
        NaturalAtomicOrbital,
        NaturalAtomicOrbitals,
        NaturalBondOrbital,
        NaturalBondOrbitals,
    )

    assert MolecularOrbitals is DomainMolecularOrbitals
    assert MoleculeOrbital is DomainMoleculeOrbital
    assert NaturalAtomicOrbital is DomainNaturalAtomicOrbital
    assert NaturalAtomicOrbitals is DomainNaturalAtomicOrbitals
    assert NaturalBondOrbital is DomainNaturalBondOrbital
    assert NaturalBondOrbitals is DomainNaturalBondOrbitals


def test_population_domain_split_preserves_legacy_reexports() -> None:
    from molop.io.base_models.data_classes.populations import (
        AtomicPopulationSeries as DomainAtomicPopulationSeries,
    )
    from molop.io.base_models.data_classes.populations import (
        ChargeSpinPopulations as DomainChargeSpinPopulations,
    )
    from molop.io.base_models.data_classes.populations import (
        TotalSpin as DomainTotalSpin,
    )
    from molop.io.base_models.DataClasses import (
        AtomicPopulationSeries,
        ChargeSpinPopulations,
        TotalSpin,
    )

    assert AtomicPopulationSeries is DomainAtomicPopulationSeries
    assert ChargeSpinPopulations is DomainChargeSpinPopulations
    assert TotalSpin is DomainTotalSpin


def test_property_domain_split_preserves_legacy_reexports() -> None:
    from molop.io.base_models.data_classes.properties import (
        BondOrders as DomainBondOrders,
    )
    from molop.io.base_models.data_classes.properties import (
        Dispersions as DomainDispersions,
    )
    from molop.io.base_models.data_classes.properties import (
        Polarizability as DomainPolarizability,
    )
    from molop.io.base_models.data_classes.properties import (
        SinglePointProperties as DomainSinglePointProperties,
    )
    from molop.io.base_models.DataClasses import (
        BondOrders,
        Dispersions,
        Polarizability,
        SinglePointProperties,
    )

    assert BondOrders is DomainBondOrders
    assert Dispersions is DomainDispersions
    assert Polarizability is DomainPolarizability
    assert SinglePointProperties is DomainSinglePointProperties


def test_status_domain_split_preserves_legacy_reexports() -> None:
    from molop.io.base_models.data_classes.status import (
        GeometryOptimizationStatus as DomainGeometryOptimizationStatus,
    )
    from molop.io.base_models.data_classes.status import (
        Status as DomainStatus,
    )
    from molop.io.base_models.DataClasses import GeometryOptimizationStatus, Status

    assert GeometryOptimizationStatus is DomainGeometryOptimizationStatus
    assert Status is DomainStatus


def test_spectroscopy_domain_split_preserves_legacy_reexports() -> None:
    from molop.io.base_models.data_classes.spectroscopy import (
        NMR as DomainNMR,
    )
    from molop.io.base_models.data_classes.spectroscopy import (
        ShieldingTensor as DomainShieldingTensor,
    )
    from molop.io.base_models.DataClasses import NMR, ShieldingTensor

    assert NMR is DomainNMR
    assert ShieldingTensor is DomainShieldingTensor


def test_solvation_domain_split_preserves_legacy_reexports() -> None:
    from molop.io.base_models.data_classes.solvation import (
        ImplicitSolvation as DomainImplicitSolvation,
    )
    from molop.io.base_models.DataClasses import ImplicitSolvation

    assert ImplicitSolvation is DomainImplicitSolvation


def test_qm_refresh_projects_structured_values_after_legacy_backfill() -> None:
    frame = BaseQMInputFrame(
        model_chemistry=QMModelChemistry(method_family="DFT", functional="B3LYP"),
        functional="PBE0",
    )

    frame.refresh_common_qm_containers()

    assert frame.model_chemistry.functional == "B3LYP"
    assert frame.functional == "B3LYP"


def test_molop_import_has_no_file_log_or_host_recursion_side_effect(tmp_path: Path) -> None:
    source_root = ROOT / "src"
    environment = os.environ.copy()
    environment["PYTHONPATH"] = os.pathsep.join(
        filter(None, (str(source_root), environment.get("PYTHONPATH")))
    )
    script = """
import os
import sys

before = sys.getrecursionlimit()
import molop

assert sys.getrecursionlimit() == before
assert not os.path.exists("molop.log")
assert "rdkit_dof" not in sys.modules
assert molop.molopconfig.log_to_file is False
"""

    result = subprocess.run(
        [sys.executable, "-c", script],
        cwd=tmp_path,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr


def test_iter_parse_outcomes_preserves_repeated_inputs_and_reports_failures(tmp_path: Path) -> None:
    missing = tmp_path / "missing.xyz"
    outcomes = list(
        iter_parse_outcomes(
            [XYZ_FIXTURE, XYZ_FIXTURE, missing],
            n_jobs=1,
            parser_detection="xyz",
        )
    )

    assert [outcome.input_index for outcome in outcomes] == [0, 1, 2]
    assert [outcome.status for outcome in outcomes] == ["ok", "ok", "missing"]
    assert outcomes[0].file_path == outcomes[1].file_path

    with pytest.raises(BatchParseError) as error:
        list(
            iter_parse_outcomes(
                [missing, XYZ_FIXTURE],
                n_jobs=1,
                parser_detection="xyz",
                fail_fast=True,
            )
        )
    assert error.value.outcome.status == "missing"
    assert error.value.outcome.input_index == 0
