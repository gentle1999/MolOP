from __future__ import annotations

from pathlib import Path
from typing import cast

import numpy as np
import pytest

from molop.io import AutoParser  # type: ignore[reportMissingImports]
from molop.io.codec_exceptions import UnsupportedFormatError
from molop.io.logic.xtb.output.models.XTBOutputFile import XTBOutputFileDisk
from molop.io.logic.xtb.output.parsers.XTBOutputFileParser import (
    XTBOutputFileParserDisk,
    XTBOutputFileParserMemory,
)


FIXTURE_DIR = Path(__file__).resolve().parent / "test_files" / "xtbout"


def test_xtbout_all_maintained_fixtures_parse_with_supported_major_versions() -> None:
    fixture_paths = sorted(FIXTURE_DIR.glob("*.out"))

    assert len(fixture_paths) >= 14
    for fixture_path in fixture_paths:
        parsed = XTBOutputFileParserDisk().parse(str(fixture_path), release_file_content=True)

        assert len(parsed.frames) == 1, fixture_path.name
        assert parsed.qm_software == "xTB", fixture_path.name
        assert parsed.qm_software_version.split(".", 1)[0] in {"5", "6"}, fixture_path.name
        assert parsed.frames[0].energies is not None, fixture_path.name
        assert parsed.frames[0].energies.total_energy is not None, fixture_path.name


def test_xtbout_legacy_v5_contract_parses_metadata_geometry_and_results() -> None:
    fixture_path = FIXTURE_DIR / "xtb_5_8_1_legacy_contract.out"
    parsed = XTBOutputFileParserDisk().parse(str(fixture_path), release_file_content=False)
    frame = parsed.frames[0]

    assert parsed.qm_software_version == "5.8.1 (legacy-contract)"
    assert parsed.input_file_name == "methane.xyz"
    assert parsed.keywords == "xtb methane.xyz --opt"
    assert parsed.method == "GFN-xTB"
    assert parsed.model_chemistry.method_family == "SEMIEMPIRICAL"
    assert parsed.request_num_cpu == 2
    assert parsed.charge == 0
    assert parsed.multiplicity == 1
    assert [task.task_type for task in parsed.task_requests] == ["opt"]
    assert frame.atoms == [6, 1, 1, 1, 1]
    assert frame.coords.shape == (5, 3)
    assert frame.coords[1, 0].m_as("angstrom") == pytest.approx(1.0821673963, abs=1.0e-9)
    assert frame.energies is not None
    assert frame.energies.total_energy.m_as("hartree") == pytest.approx(-4.1752185)
    assert frame.molecular_orbitals is not None
    assert len(frame.molecular_orbitals.alpha_energies) == 8
    assert frame.geometry_optimization_status is not None
    assert frame.geometry_optimization_status.geometry_optimized is True
    assert frame.status is not None
    assert frame.status.scf_converged is True
    assert frame.status.normal_terminated is True
    assert frame.running_time.m_as("second") == pytest.approx(0.125)


def test_xtbout_major_seven_uses_the_modern_six_contract() -> None:
    source = (FIXTURE_DIR / "dsgdb9nsd_130336-5_opt.out").read_text(encoding="utf-8")
    future_source = source.replace("* xtb version 6.6.1", "* xtb version 7.0.0", 1)

    parsed_v6 = XTBOutputFileParserMemory().parse(source)
    parsed_v7 = XTBOutputFileParserMemory().parse(future_source)
    frame_v6 = parsed_v6.frames[0]
    frame_v7 = parsed_v7.frames[0]

    assert parsed_v7.qm_software_version.startswith("7.0.0")
    assert parsed_v7.method == parsed_v6.method
    assert parsed_v7.charge == parsed_v6.charge
    assert parsed_v7.multiplicity == parsed_v6.multiplicity
    assert parsed_v7.task_requests == parsed_v6.task_requests
    assert frame_v7.atoms == frame_v6.atoms
    assert np.array_equal(frame_v7.coords.magnitude, frame_v6.coords.magnitude)
    assert frame_v7.energies == frame_v6.energies
    assert frame_v7.molecular_orbitals is not None
    assert frame_v6.molecular_orbitals is not None
    np.testing.assert_allclose(
        frame_v7.molecular_orbitals.alpha_energies.magnitude,
        frame_v6.molecular_orbitals.alpha_energies.magnitude,
        equal_nan=True,
    )
    np.testing.assert_allclose(
        frame_v7.molecular_orbitals.alpha_occupancies,
        frame_v6.molecular_orbitals.alpha_occupancies,
        equal_nan=True,
    )
    assert frame_v7.charge_spin_populations is not None
    assert frame_v6.charge_spin_populations is not None
    np.testing.assert_allclose(
        frame_v7.charge_spin_populations["mulliken_charges"].values,
        frame_v6.charge_spin_populations["mulliken_charges"].values,
    )
    assert frame_v7.geometry_optimization_status == frame_v6.geometry_optimization_status
    assert frame_v7.gradient_norm == frame_v6.gradient_norm
    assert frame_v7.status == frame_v6.status


def test_xtbout_modern_single_point_keeps_results_without_embedded_geometry() -> None:
    fixture_path = FIXTURE_DIR / "dsgdb9nsd_130336-5.out"
    parsed = XTBOutputFileParserDisk().parse(str(fixture_path))
    frame = parsed.frames[0]

    assert parsed.qm_software_version == "6.6.1"
    assert parsed.method == "GFN1-xTB"
    assert parsed.multiplicity == 2
    assert [task.task_type for task in parsed.task_requests] == ["sp"]
    assert frame.atoms == []
    assert frame.energies is not None
    assert frame.energies.total_energy.m_as("hartree") == pytest.approx(-26.296158842419)
    assert frame.molecular_orbitals is not None
    assert len(frame.molecular_orbitals.alpha_energies) == 44
    assert frame.charge_spin_populations is not None
    assert len(frame.charge_spin_populations["mulliken_charges"].values) == 13
    assert frame.charge_spin_populations["mulliken_charges"].source_label == "Mulliken/CM5 charges"
    assert frame.charge_spin_populations["cm5_charges"].scheme == "cm5"
    assert frame.single_point_properties is not None
    assert frame.single_point_properties.vip.m_as("eV/particle") == pytest.approx(7.9982)


@pytest.mark.parametrize(
    ("fixture_name", "expected_energy", "expected_gradient_norm"),
    [
        ("dsgdb9nsd_130336-5_opt.out", -25.324456612331, 0.0005850),
        ("dsgdb9nsd_130336-5_sdf_opt.out", -25.324456688174, 0.0005531),
    ],
)
def test_xtbout_modern_optimization_parses_xyz_and_sdf_final_structures(
    fixture_name: str,
    expected_energy: float,
    expected_gradient_norm: float,
) -> None:
    parsed = XTBOutputFileParserDisk().parse(str(FIXTURE_DIR / fixture_name))
    frame = parsed.frames[0]

    assert frame.atoms == [7, 1, 1, 1, 6, 6, 6, 6, 6, 7, 7, 8, 1]
    assert frame.coords.shape == (13, 3)
    assert np.isfinite(frame.coords.magnitude).all()
    assert frame.energies is not None
    assert frame.energies.total_energy.m_as("hartree") == pytest.approx(expected_energy)
    assert frame.geometry_optimization_status is not None
    assert frame.geometry_optimization_status.geometry_optimized is True
    assert frame.gradient_norm is not None
    assert frame.gradient_norm.m_as("hartree/bohr") == pytest.approx(expected_gradient_norm)
    assert frame.polarizability is not None
    assert frame.polarizability.dipole.shape == (3,)
    assert frame.rotation_constants.shape == (3,)
    assert frame.thermal_informations is None


def test_xtbout_frequency_and_thermochemistry_are_structured() -> None:
    parsed = XTBOutputFileParserDisk().parse(str(FIXTURE_DIR / "6-6-1-hess.out"))
    frame = parsed.frames[0]

    assert [task.task_type for task in parsed.task_requests] == ["opt", "freq"]
    assert frame.vibrations is not None
    assert len(frame.vibrations) == 36
    assert frame.vibrations.num_imaginary == 1
    assert frame.vibrations.frequencies[0].m_as("cm^-1") == pytest.approx(-193.94)
    assert len(frame.vibrations.reduced_masses) == 36
    assert len(frame.vibrations.IR_intensities) == 36
    assert frame.thermal_informations is not None
    assert frame.thermal_informations.ZPVE.m_as("hartree/particle") == pytest.approx(0.094224490746)
    assert frame.thermal_informations.H_T.m_as("hartree/particle") == pytest.approx(
        -26.415301029112
    )
    assert frame.thermal_informations.G_T.m_as("hartree/particle") == pytest.approx(
        -26.455552821913
    )
    assert frame.thermal_informations.C_V.m_as("cal/mol/K") == pytest.approx(29.8614)
    assert frame.thermal_informations.S.m_as("cal/mol/K") == pytest.approx(84.7170)


def test_xtbout_fukui_indices_are_structured() -> None:
    parsed = XTBOutputFileParserDisk().parse(str(FIXTURE_DIR / "dsgdb9nsd_130336-5_fukui.out"))
    properties = parsed.frames[0].single_point_properties

    assert properties is not None
    assert len(properties.fukui_positive) == 13
    assert len(properties.fukui_negative) == 13
    assert len(properties.fukui_zero) == 13
    assert properties.fukui_positive[0] == pytest.approx(0.004)


def test_xtbout_vipea_and_gei_properties_are_structured() -> None:
    parsed = XTBOutputFileParserDisk().parse(str(FIXTURE_DIR / "dsgdb9nsd_123127-1-vipea.out"))
    properties = parsed.frames[0].single_point_properties

    assert properties is not None
    assert properties.vip.m_as("eV/particle") == pytest.approx(4.3218)
    assert properties.vea.m_as("eV/particle") == pytest.approx(-2.5636)

    gei_text = """
     |                           x T B                           |
     |               Version 6.6.1 (property-contract)          |
program call               : xtb molecule.xyz --vipea
coordinate file            : molecule.xyz
   *** convergence criteria satisfied after 3 iterations ***
total E       :      -4.0000000
Global electrophilicity index (eV):    1.2345
 * finished run on 2024/01/01 at 00:00:00.000
"""
    gei_frame = XTBOutputFileParserMemory().parse(gei_text).frames[0]

    assert gei_frame.single_point_properties is not None
    assert gei_frame.single_point_properties.gei.m_as("eV/particle") == pytest.approx(1.2345)


def test_autoparser_detects_xtbout_and_preserves_source_evidence() -> None:
    fixture_path = FIXTURE_DIR / "xtb_6_3_2_opt.out"
    source_text = fixture_path.read_text(encoding="utf-8")
    batch = AutoParser(str(fixture_path), capture_source_evidence=True)
    parsed = cast(XTBOutputFileDisk, batch[0])
    frame = parsed.frames[0]

    assert parsed.detected_format_id == "xtbout"
    assert parsed.source_format == "xtbout"
    assert parsed.source_segments[0].source_span.start_char == 0
    assert parsed.source_segments[0].source_span.end_char == len(source_text)
    assert frame.source_span.start_char == 0
    assert frame.source_span.end_char == len(source_text)
    assert frame.coordinate_source == "observed"
    assert frame.energies is not None
    assert len(frame.energies.observations) == 1
    assert frame.energies.observations[0].source_label.startswith("| TOTAL ENERGY")


@pytest.mark.parametrize("major", [0, 4])
def test_xtbout_explicitly_rejects_versions_before_five(major: int) -> None:
    text = f"""
     |                           x T B                           |
     |               Version {major}.0.0 (unsupported)               |
"""

    with pytest.raises(UnsupportedFormatError, match="below the supported minimum 5"):
        XTBOutputFileParserMemory().parse(text)


def test_xtbout_probe_does_not_claim_gaussian_or_orca_outputs() -> None:
    g16_path = Path(__file__).resolve().parent / "test_files" / "g16log" / "H2O.log"
    orca_path = Path(__file__).resolve().parent / "test_files" / "orca" / "H2_sp_orca.out"

    assert XTBOutputFileParserDisk.probe_file_format(g16_path) is False
    assert XTBOutputFileParserDisk.probe_file_format(orca_path) is False
