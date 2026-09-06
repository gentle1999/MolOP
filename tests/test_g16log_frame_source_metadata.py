from pathlib import Path

import numpy as np
import pytest

from molop.io.base_models.DataClasses import EnergyObservation, Vibrations
from molop.io.logic.gaussian.log.frame_models.G16LogFileFrame import G16LogFileFrameMemory
from molop.io.logic.gaussian.log.frame_parsers._g16_extractors import (
    merge_g16_energy_payloads,
)
from molop.io.logic.gaussian.log.frame_parsers.G16LogFileFrameParser import (
    G16LogFileFrameParserMemory,
)
from molop.io.logic.gaussian.log.locators import (
    locate_g16_section_frames,
    locate_g16_sections,
)
from molop.unit import atom_ureg


FIXTURE_ROOT = Path(__file__).resolve().parent / "test_files" / "g16log"


def _frame_blocks(name: str) -> list[str]:
    text = (FIXTURE_ROOT / name).read_text()
    return [
        frame.text(text)
        for section in locate_g16_sections(text)
        for frame in locate_g16_section_frames(text, section)
    ]


def test_gaussian_frame_source_metadata_is_opt_in() -> None:
    block = next(frame for frame in _frame_blocks("3-m-Py_anion_Opt.log") if "Converged?" in frame)

    plain = G16LogFileFrameParserMemory().parse(block)
    captured = G16LogFileFrameParserMemory(capture_source_evidence=True).parse(block)

    assert plain.energies is not None
    assert plain.energies.observations == []
    assert plain.coordinate_source is None
    assert plain.geometry_optimization_status is not None
    assert plain.geometry_optimization_status.source_converged is None
    assert plain.geometry_optimization_status.source_labels is None
    assert "coordinate_source" not in plain.model_dump()
    assert "observations" not in plain.energies.model_dump()
    assert "source_converged" not in plain.geometry_optimization_status.model_dump()
    assert "source_labels" not in plain.geometry_optimization_status.model_dump()

    assert captured.energies is not None
    assert captured.energies.observations
    assert captured.coordinate_source == "observed"
    assert captured.geometry_optimization_status is not None
    assert captured.geometry_optimization_status.source_converged
    assert captured.geometry_optimization_status.source_labels
    assert captured.geometry_optimization_status.max_force is not None
    assert str(captured.geometry_optimization_status.max_force.units) == "hartree / bohr"
    assert captured.geometry_optimization_status.max_displacement is not None
    assert str(captured.geometry_optimization_status.max_displacement.units) == "bohr"
    assert captured.model_dump()["coordinate_source"] == "observed"
    assert captured.energies.model_dump()["observations"]


def test_gaussian_vibration_orientation_transform_does_not_translate_modes() -> None:
    coordinates = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 2.0, 0.0],
            [0.0, 0.0, 3.0],
        ]
    )
    rotation = np.array(
        [
            [0.0, -1.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    translation = np.array([10.0, 20.0, 30.0])
    standard_coordinates = (rotation @ coordinates.T).T + translation
    standard_mode = np.array(
        [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [-1.0, -1.0, -1.0],
        ]
    )

    frame = G16LogFileFrameMemory(
        atoms=[1, 6, 8, 1],
        coords=coordinates * atom_ureg.angstrom,
        standard_coords=standard_coordinates * atom_ureg.angstrom,
        vibrations=Vibrations(
            frequencies=np.full(6, 100.0) * atom_ureg.cm_1,
            vibration_modes=[standard_mode * atom_ureg.angstrom for _ in range(6)],
        ),
    )

    assert frame.vibrations is not None
    assert frame.standard_orientation_transformation_matrix is not None
    np.testing.assert_allclose(
        frame.standard_orientation_transformation_matrix[:3, :3],
        rotation,
        atol=1.0e-12,
    )
    np.testing.assert_allclose(
        frame.standard_orientation_transformation_matrix[:3, 3],
        translation,
        atol=1.0e-12,
    )
    assert frame.standard_coords is not None
    np.testing.assert_allclose(frame.standard_coords.magnitude, standard_coordinates, atol=1.0e-12)
    np.testing.assert_allclose(
        frame.vibrations.vibration_modes[0].magnitude,
        standard_mode @ rotation,
        atol=1.0e-12,
    )


def test_gaussian_archive_observations_reuse_parsed_energies() -> None:
    block = _frame_blocks("1.log")[-1]

    plain = G16LogFileFrameParserMemory().parse(block)
    captured = G16LogFileFrameParserMemory(capture_source_evidence=True).parse(block)

    assert plain.energies is not None
    assert plain.energies.observations == []
    assert captured.energies is not None
    assert {observation.method for observation in captured.energies.observations} >= {
        "reference",
        "CCSD",
    }
    assert all(
        observation.source_label.startswith("Gaussian archive")
        for observation in captured.energies.observations
    )


def test_gaussian_terminal_archive_energy_has_matching_typed_observation() -> None:
    block = _frame_blocks("000000000000_000016928457_00_conf_01_ts.107c60f3cfcb.log")[-1]

    frame = G16LogFileFrameParserMemory(capture_source_evidence=True).parse(block)

    assert frame.energies is not None
    assert frame.energies.reference_energy is not None
    reference_observations = [
        observation
        for observation in frame.energies.observations
        if observation.method == "reference" and observation.quantity_semantics == "total_energy"
    ]
    assert len(reference_observations) == 1
    assert reference_observations[0].source_label == "Gaussian archive reference"
    assert reference_observations[0].value.to("hartree").m == pytest.approx(
        frame.energies.reference_energy.to("hartree").m
    )


def test_gaussian_live_and_archive_energy_observations_merge_deterministically() -> None:
    frame = G16LogFileFrameParserMemory(capture_source_evidence=True).parse(
        _frame_blocks("test_ccsd_t.log")[-1]
    )

    assert frame.energies is not None
    assert frame.energies.reference_energy is not None
    assert frame.energies.reference_energy.to("hartree").m == pytest.approx(-0.499821176024)
    assert frame.energies.ccsd_energy is not None
    assert frame.energies.ccsd_t_energy is not None
    assert [
        (observation.method, observation.source_label)
        for observation in frame.energies.observations
    ] == [
        ("reference", "SCF Done"),
        ("MP2", "EUMP2"),
        ("MP3", "EUMP3"),
        ("MP4", "EUMP4"),
        ("CCSD", "Wavefunction amplitudes converged. E(Corr)"),
        ("CCSD(T)", "CCSD(T)"),
        ("reference", "Gaussian archive reference"),
        ("MP2", "Gaussian archive MP2"),
        ("MP3", "Gaussian archive MP3"),
        ("MP4", "Gaussian archive MP4"),
        ("CCSD", "Gaussian archive CCSD"),
        ("CCSD(T)", "Gaussian archive CCSD(T)"),
    ]

    method_by_field = {
        "reference_energy": "reference",
        "mp2_energy": "MP2",
        "mp3_energy": "MP3",
        "mp4_energy": "MP4",
        "ccsd_energy": "CCSD",
        "ccsd_t_energy": "CCSD(T)",
    }
    for field_name, method in method_by_field.items():
        scalar = getattr(frame.energies, field_name)
        assert scalar is not None
        assert any(
            observation.method == method
            and observation.quantity_semantics == "total_energy"
            and observation.value.to("hartree").m == pytest.approx(scalar.to("hartree").m)
            for observation in frame.energies.observations
        )


def test_g16_energy_payload_merge_deduplicates_full_observation_identity() -> None:
    live_observation = EnergyObservation(
        method="reference",
        quantity_semantics="total_energy",
        value=-1.0 * atom_ureg.hartree,
        source_label="SCF Done",
    )
    archive_observation = EnergyObservation(
        method="reference",
        quantity_semantics="total_energy",
        value=-1.1 * atom_ureg.hartree,
        source_label="Gaussian archive reference",
    )
    fallback = {
        "reference_energy": -1.1 * atom_ureg.hartree,
        "mp2_energy": -1.2 * atom_ureg.hartree,
        "observations": [live_observation.model_copy(), archive_observation],
    }

    merged = merge_g16_energy_payloads(
        {
            "reference_energy": -1.0 * atom_ureg.hartree,
            "observations": [live_observation],
        },
        fallback,
    )
    merged_twice = merge_g16_energy_payloads(merged, fallback)

    assert merged["reference_energy"].to("hartree").m == pytest.approx(-1.0)
    assert merged["mp2_energy"].to("hartree").m == pytest.approx(-1.2)
    assert [observation.source_label for observation in merged["observations"]] == [
        "SCF Done",
        "Gaussian archive reference",
    ]
    assert merged_twice["observations"] == merged["observations"]
