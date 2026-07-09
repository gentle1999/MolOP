from pathlib import Path

import pytest

from molop import AutoParser
from molop.io.logic.QM_frame_models.G16V3Components import (
    _extract_labeled_float_tokens as _extract_v3_labeled_float_tokens,
)
from molop.io.logic.QM_frame_parsers._g16_v2_extractors import (
    ParseState,
    extract_archive_tail_payload_from_state,
)
from molop.io.logic.QM_frame_parsers._g16_v2_shared import _extract_labeled_float_tokens
from molop.io.logic.QM_parsers.G16LogFileParser import G16LogFileParserMemory


MERGED_FREQUENCY_FIXTURE = (
    Path(__file__).resolve().parent / "test_files" / "g16log" / "3-m-Py_anion_Opt.log"
)
TERMINAL_ARCHIVE_ENERGY_FIXTURE = (
    Path(__file__).resolve().parent
    / "test_files"
    / "g16log"
    / "000000000000_000016928457_00_conf_01_ts.107c60f3cfcb.log"
)


def test_default_g16log_parser_handles_concatenated_frequency_values() -> None:
    batch = AutoParser(str(MERGED_FREQUENCY_FIXTURE), n_jobs=1)
    frame = batch[0][-1]
    vibrations = frame.vibrations

    assert vibrations is not None
    assert len(vibrations.frequencies) == 51
    assert len(vibrations.reduced_masses) == 51
    assert len(vibrations.force_constants) == 51
    assert len(vibrations.IR_intensities) == 51


def test_extract_labeled_float_tokens_handles_concatenated_polarizability_values() -> None:
    approx_line = "Approx polarizability: 907.677  15.5261217.377 -79.003-102.3351272.481"
    before_force_line = "Polarizability= 7.08599298D+00-1.00725577D+00 6.37248853D+00"

    assert _extract_labeled_float_tokens(
        approx_line, "Approx polarizability:", expected_count=6, decimal_places=3
    ) == [
        907.677,
        15.526,
        1217.377,
        -79.003,
        -102.335,
        1272.481,
    ]
    assert _extract_labeled_float_tokens(
        before_force_line, "Polarizability=", expected_count=3, decimal_places=8
    ) == [
        7.08599298,
        -1.00725577,
        6.37248853,
    ]
    assert _extract_v3_labeled_float_tokens(
        before_force_line, "Polarizability=", expected_count=3, decimal_places=8
    ) == [
        7.08599298,
        -1.00725577,
        6.37248853,
    ]


def test_g16log_file_running_time_is_accumulated_into_model_field() -> None:
    fixture = Path(__file__).resolve().parent / "test_files" / "g16log" / "H2O.log"
    parser = G16LogFileParserMemory()
    parsed = parser.parse(fixture.read_text())

    assert parsed.running_time is not None
    assert parsed.running_time.to("second").m == 102.5


def test_default_g16log_parser_uses_archive_energy_for_terminal_frequency_frame() -> None:
    batch = AutoParser(str(TERMINAL_ARCHIVE_ENERGY_FIXTURE), parser_detection="g16log", n_jobs=1)
    file_model = batch[0]
    last_frame = file_model[-1]

    assert file_model.is_error is False
    assert last_frame.is_error is False
    assert last_frame.energies is not None
    assert last_frame.energies.reference_energy is not None
    assert last_frame.energies.reference_energy.to("hartree").m == pytest.approx(-502.1233491)


def test_g16log_archive_tail_payload_extracts_thermal_and_polar_fields() -> None:
    frames = G16LogFileParserMemory()._split_file(TERMINAL_ARCHIVE_ENERGY_FIXTURE.read_text())

    optimization_tail = extract_archive_tail_payload_from_state(ParseState(frames[25]))
    optimization_polar = optimization_tail["polarizability"]
    assert optimization_polar["dipole"].to("debye").magnitude.tolist() == pytest.approx(
        [-0.9067916, 0.1321581, -0.9767749]
    )
    assert optimization_polar["quadrupole"].to("debye * angstrom").magnitude.tolist() == (
        pytest.approx([-1.3650278, 3.8055608, -2.4405331, -0.0883569, -4.3998293, -0.0069285])
    )

    frequency_tail = extract_archive_tail_payload_from_state(ParseState(frames[-1]))
    frequency_thermal = frequency_tail["thermal_informations"]
    assert frequency_thermal["ZPVE"].to("hartree / particle").magnitude == pytest.approx(0.2847577)
    assert frequency_thermal["TCE"].to("hartree / particle").magnitude == pytest.approx(0.2990909)

    frequency_polar = frequency_tail["polarizability"]
    assert frequency_polar["dipole"].to("debye").magnitude.tolist() == pytest.approx(
        [-0.6597763, -0.0350057, -1.2277048]
    )
    assert frequency_polar["polarizability_tensor"].to("bohr ** 3").magnitude.tolist() == (
        pytest.approx([127.3525155, -2.1083289, 150.7303896, -6.5742126, 3.3900566, 119.0969191])
    )
