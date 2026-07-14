from pathlib import Path

from molop.io.base_models.FileParser import BaseFileParser
from molop.io.logic.orca.log.parsers.ORCALogFileParser import (
    ORCALogFileParserMemory,
    ORCALogFileParserMixin,
)


ORCA_SINGLE_POINT = (
    Path(__file__).resolve().parent
    / "test_files"
    / "orca"
    / "output_files"
    / "local"
    / "H2_sp_orca.out"
)
ORCA_OPTIMIZATION = (
    Path(__file__).resolve().parent
    / "test_files"
    / "orca"
    / "output_files"
    / "local"
    / "opt_orca.out"
)


def _multi_job_output() -> str:
    first_job = ORCA_SINGLE_POINT.read_text(encoding="utf-8")
    second_job = first_job.replace("H2_sp.inp", "second.inp").replace(
        "! SP PBE def2-SVP",
        "! SP HF def2-TZVP",
    )
    second_job = "".join(
        line for line in second_job.splitlines(keepends=True) if "Program Version" not in line
    )
    return f"{first_job}\n********** JOB NUMBER 2 **********\n{second_job}"


def _multi_job_output_with_final_job_without_coordinates() -> str:
    first_job = ORCA_SINGLE_POINT.read_text(encoding="utf-8")
    second_job = first_job.replace("H2_sp.inp", "second.inp").replace(
        "! SP PBE def2-SVP",
        "! SP HF def2-TZVP",
    )
    second_job = "".join(
        line
        for line in second_job.splitlines(keepends=True)
        if "Program Version" not in line and not line.startswith("  H      0.000000    0.000000")
    )
    return f"{first_job}\n********** JOB NUMBER 2 **********\n{second_job}"


def test_orca_log_file_parser_uses_base_parse_lifecycle() -> None:
    assert "_parse" not in ORCALogFileParserMixin.__dict__
    assert "_parse_metadata" not in ORCALogFileParserMixin.__dict__
    assert "_split_file" not in ORCALogFileParserMixin.__dict__
    assert ORCALogFileParserMemory._parse is BaseFileParser._parse


def test_orca_multi_job_frames_keep_job_local_metadata_and_artifact_version() -> None:
    parsed = ORCALogFileParserMemory(capture_source_evidence=True).parse(_multi_job_output())

    assert len(parsed) == 2
    assert parsed.input_file_name == "H2_sp.inp"
    assert parsed.qm_software_version == "4.1.1"

    first, second = parsed.frames
    assert (first.input_file_name, first.method, first.basis_set) == (
        "H2_sp.inp",
        "DFT",
        "def2-SVP",
    )
    assert (second.input_file_name, second.method, second.basis_set) == (
        "second.inp",
        "HF",
        "def2-TZVP",
    )
    assert first.qm_software_version == second.qm_software_version == "4.1.1"
    assert (first.segment_index, first.segment_frame_index, first.frame_role) == (
        0,
        0,
        "single_point",
    )
    assert first.file_frame_index == 0
    assert (second.segment_index, second.segment_frame_index, second.frame_role) == (
        1,
        0,
        "single_point",
    )
    assert second.file_frame_index == 1
    assert parsed.source_segments[0].protocol is not None
    assert parsed.source_segments[0].protocol["functional"] == "PBE"
    assert parsed.source_segments[0].task_requests[0]["task_type"] == "sp"
    assert parsed.source_segments[1].protocol is not None
    assert parsed.source_segments[1].protocol["method_family"] == "HF"
    assert parsed.source_segments[1].protocol["basis_set"] == "def2-TZVP"
    assert parsed.source_segments[1].task_requests[0]["task_type"] == "sp"


def test_orca_only_last_uses_last_job_metadata_and_true_segment_index() -> None:
    parsed = ORCALogFileParserMemory(
        capture_source_evidence=True,
        only_last_frame=True,
    ).parse(_multi_job_output())

    assert len(parsed) == 1
    assert parsed.input_file_name == "second.inp"
    assert parsed.qm_software_version == "4.1.1"
    assert parsed[0].input_file_name == "second.inp"
    assert parsed[0].segment_index == 1
    assert parsed[0].segment_frame_index == 0
    assert parsed[0].file_frame_index == 1
    assert parsed.source_segments[0]["segment_index"] == 1


def test_orca_multi_job_retains_final_segment_without_parseable_coordinates() -> None:
    parsed = ORCALogFileParserMemory(capture_source_evidence=True).parse(
        _multi_job_output_with_final_job_without_coordinates()
    )

    assert len(parsed) == 1
    assert len(parsed.source_segments) == 2
    assert parsed.qm_software_version == "4.1.1"

    frame = parsed[0]
    assert (frame.input_file_name, frame.method, frame.basis_set) == (
        "H2_sp.inp",
        "DFT",
        "def2-SVP",
    )
    assert frame.qm_software_version == "4.1.1"

    failed_segment = parsed.source_segments[1]
    assert failed_segment.frame_count == 0
    assert failed_segment.captured_frame_indices == []
    assert failed_segment.task_types == ["sp"]
    assert failed_segment.protocol is not None
    assert failed_segment.protocol["method_family"] == "HF"
    assert failed_segment.protocol["basis_set"] == "def2-TZVP"
    assert failed_segment.task_requests[0]["task_type"] == "sp"
    assert failed_segment.termination_status is True
    assert failed_segment.scf_status is True
    assert failed_segment.qm_software_version == "4.1.1"
    assert failed_segment.parse_presence["geometry"] == "absent_in_source"
    assert any(
        diagnostic.code == "MOL.PARSE.SEGMENT_FRAMES_ABSENT"
        for diagnostic in failed_segment.diagnostics
    )


def test_orca_single_point_request_rejects_optimization_evidence() -> None:
    source = ORCA_OPTIMIZATION.read_text(encoding="utf-8").replace(
        "! Opt PBE0 RIJCOSX D3BJ def2-SVP def2/J CPCM(Water)",
        "! SP PBE0 RIJCOSX D3BJ def2-SVP def2/J CPCM(Water)",
    )

    parsed = ORCALogFileParserMemory(capture_source_evidence=True).parse(source)

    assert parsed.source_segments[0].task_types == ["sp"]
    optimized_frames = [frame for frame in parsed if frame.geometry_optimization_status is not None]
    assert optimized_frames
    for frame in optimized_frames:
        diagnostic = next(
            diagnostic
            for diagnostic in frame.parse_diagnostics
            if diagnostic.code == "MOL.CALC.UNEXPECTED_OPTIMIZATION"
        )
        assert diagnostic.severity == "error"
        assert diagnostic.scope == "frame"
