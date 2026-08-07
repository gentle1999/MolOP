from __future__ import annotations

from hashlib import sha256
from pathlib import Path
from typing import Any

import pytest

from molop.io.base_models.source import SourceSpan
from molop.io.logic.gaussian.log.parsers.G16LogFileParser import G16LogFileParserMemory
from molop.io.logic.orca.log.parsers.ORCALogFileParser import ORCALogFileParserMemory


FIXTURE_ROOT = Path(__file__).resolve().parent / "test_files"
FORMAT_CASES = (
    (G16LogFileParserMemory, FIXTURE_ROOT / "g16log" / "H2O.log"),
    (ORCALogFileParserMemory, FIXTURE_ROOT / "orca" / "opt_orca.out"),
)

_FRAME_EVIDENCE_FIELDS = {
    "source_span",
    "source_block_sha256",
    "segment_index",
    "segment_frame_index",
    "file_frame_index",
    "parse_presence",
    "parse_diagnostics",
    "parse_completeness",
    "coordinate_source",
    "coordinate_provenance",
    "coordinate_decimal_places",
    "frame_role",
    "force_source_field",
    "force_transformation",
    "source_to_topology_atom_permutation",
    "topology_reconstruction_backend",
    "topology_make_dative_bonds",
    "topology_reconstruction_config_sha256",
    "topology_reconstruction_status",
}


def _scientific_frame_dump(frame: Any) -> dict[str, Any]:
    include = set(type(frame).model_fields) - _FRAME_EVIDENCE_FIELDS
    return frame.to_unitless_dump_with_unit_keys(
        include=include,
        exclude={
            "energies": {"observations"},
            "geometry_optimization_status": {"source_converged", "source_labels"},
        },
        exclude_none=True,
    )


def _span_dump(span: SourceSpan | dict[str, int]) -> dict[str, int]:
    return span.model_dump() if isinstance(span, SourceSpan) else span


def _assert_hashed_span(
    raw_bytes: bytes,
    span: SourceSpan | dict[str, int],
    digest: str,
) -> None:
    dumped = _span_dump(span)
    assert dumped["start_byte"] < dumped["end_byte"] <= len(raw_bytes)
    assert sha256(raw_bytes[dumped["start_byte"] : dumped["end_byte"]]).hexdigest() == digest


@pytest.mark.parametrize(("parser_cls", "fixture"), FORMAT_CASES)
def test_calculation_source_capture_is_disabled_by_default(
    parser_cls: type[Any],
    fixture: Path,
) -> None:
    parsed = parser_cls().parse(fixture.read_bytes().decode("utf-8"))

    assert parsed.artifact_sha256 is None
    assert parsed.artifact_size_bytes is None
    assert parsed.source_encoding is None
    assert parsed.source_segments == []
    assert parsed.source_complete is None
    assert all(frame.source_span is None for frame in parsed.frames)
    assert all(frame.source_block_sha256 is None for frame in parsed.frames)
    assert all(frame.segment_index is None for frame in parsed.frames)
    assert all(frame.segment_frame_index is None for frame in parsed.frames)
    assert all(frame.file_frame_index is None for frame in parsed.frames)
    assert all(frame.frame_role is None for frame in parsed.frames)


@pytest.mark.parametrize(("parser_cls", "fixture"), FORMAT_CASES)
def test_source_capture_preserves_frame_boundaries_and_scientific_payload(
    parser_cls: type[Any],
    fixture: Path,
) -> None:
    text = fixture.read_text(encoding="utf-8")
    plain = parser_cls().parse(text)
    captured = parser_cls(capture_source_evidence=True).parse(text)

    assert len(plain.frames) == len(captured.frames)
    assert [frame.frame_content for frame in plain.frames] == [
        frame.frame_content for frame in captured.frames
    ]
    assert [_scientific_frame_dump(frame) for frame in plain.frames] == [
        _scientific_frame_dump(frame) for frame in captured.frames
    ]


@pytest.mark.parametrize(("parser_cls", "fixture"), FORMAT_CASES)
def test_calculation_source_capture_uses_exact_parser_blocks(
    parser_cls: type[Any],
    fixture: Path,
) -> None:
    raw_bytes = fixture.read_bytes()
    parsed = parser_cls(capture_source_evidence=True).parse(raw_bytes.decode("utf-8"))

    assert parsed.artifact_sha256 == sha256(raw_bytes).hexdigest()
    assert parsed.artifact_size_bytes == len(raw_bytes)
    assert parsed.source_encoding == "utf-8"
    assert parsed.source_complete is True
    assert parsed.source_segments

    for expected_segment_index, segment in enumerate(parsed.source_segments):
        assert segment["segment_index"] == expected_segment_index
        _assert_hashed_span(
            raw_bytes,
            segment["source_span"],
            segment["source_block_sha256"],
        )

    for expected_file_frame_index, frame in enumerate(parsed.frames):
        assert frame.source_span is not None
        assert frame.source_block_sha256 is not None
        assert frame.segment_index is not None
        assert frame.segment_frame_index is not None
        assert frame.file_frame_index == expected_file_frame_index
        assert frame.frame_role in {"initial", "intermediate", "terminal", "single_point"}
        _assert_hashed_span(raw_bytes, frame.source_span, frame.source_block_sha256)
        frame_span = _span_dump(frame.source_span)
        start_char = frame_span["start_char"]
        end_char = frame_span["end_char"]
        assert frame.frame_content == parsed.file_content[start_char:end_char]


def test_gaussian_crlf_parsing_preserves_exact_source_evidence() -> None:
    lf_text = (FIXTURE_ROOT / "g16log" / "H2O.log").read_bytes().decode("utf-8")
    lf_text = lf_text.replace("\r\n", "\n").replace("\r", "\n")
    crlf_text = lf_text.replace("\n", "\r\n")
    raw_bytes = crlf_text.encode("utf-8")

    parsed = G16LogFileParserMemory(capture_source_evidence=True).parse(crlf_text)

    assert parsed.artifact_sha256 == sha256(raw_bytes).hexdigest()
    assert parsed.artifact_size_bytes == len(raw_bytes)
    assert parsed.frames
    for frame in parsed.frames:
        assert frame.source_span is not None
        assert frame.source_block_sha256 is not None
        _assert_hashed_span(raw_bytes, frame.source_span, frame.source_block_sha256)
        assert (
            frame.frame_content
            == crlf_text[frame.source_span.start_char : frame.source_span.end_char]
        )
        if frame.forces is not None:
            assert frame.forces.shape == (len(frame.atoms), 3)


@pytest.mark.parametrize(("parser_cls", "fixture"), FORMAT_CASES)
def test_only_last_frame_retains_true_source_indices(
    parser_cls: type[Any],
    fixture: Path,
) -> None:
    text = fixture.read_bytes().decode("utf-8")
    complete = parser_cls(capture_source_evidence=True).parse(text)
    partial = parser_cls(
        capture_source_evidence=True,
        only_last_frame=True,
    ).parse(text, release_file_content=True)

    expected = complete.frames[-1]
    actual = partial.frames[0]
    assert partial.source_complete is False
    assert len(partial.frames) == 1
    assert actual.source_span == expected.source_span
    assert actual.source_block_sha256 == expected.source_block_sha256
    assert actual.segment_index == expected.segment_index
    assert actual.segment_frame_index == expected.segment_frame_index
    assert actual.file_frame_index == expected.file_frame_index
    assert actual.frame_role == expected.frame_role
    partial_segment = partial.source_segments[0]
    complete_segment = next(
        segment
        for segment in complete.source_segments
        if segment.segment_index == actual.segment_index
    )
    assert partial_segment.captured_frame_indices == [actual.segment_frame_index]
    assert partial_segment.model_dump(exclude={"captured_frame_indices"}) == (
        complete_segment.model_dump(exclude={"captured_frame_indices"})
    )
    assert partial.file_content == ""
    assert actual.frame_content == ""
