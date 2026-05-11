from pathlib import Path

from molop import AutoParser
from molop.io.logic.QM_frame_models.G16V3Components import (
    _extract_labeled_float_tokens as _extract_v3_labeled_float_tokens,
)
from molop.io.logic.QM_frame_parsers._g16_v2_shared import _extract_labeled_float_tokens
from molop.io.logic.QM_parsers.G16LogFileParser import G16LogFileParserMemory


MERGED_FREQUENCY_FIXTURE = Path(__file__).resolve().parent / "test_files" / "g16log" / "3-m-Py_anion_Opt.log"


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
