from __future__ import annotations

from hashlib import sha256
from pathlib import Path
from typing import Any

import pytest
from rdkit import Chem

from molop.io.base_models.source import LocatedTextBlock
from molop.io.codec_exceptions import ParseError
from molop.io.logic.coords.parsers.SDFFileParser import SDFFileParserMemory
from molop.io.logic.coords.parsers.SMIFileParser import SMIFileParserMemory
from molop.io.logic.coords.parsers.XYZFileParser import XYZFileParserMemory
from molop.io.logic.gaussian.input.parsers._gjf_file_extractors import (
    GJF_INCLUDE_PROVENANCE_DIAGNOSTIC,
    locate_gjf_link1_frames,
)
from molop.io.logic.gaussian.input.parsers.GJFFileParser import (
    GJFFileParserDisk,
    GJFFileParserMemory,
)
from molop.io.logic.orca.input.parsers._orca_inp_file_extractors import (
    locate_orca_input_frames,
)
from molop.io.logic.orca.input.parsers.ORCAInpFileParser import ORCAInpFileParserMemory


def _assert_exact_frame_blocks(
    parsed: Any,
    source: str,
    blocks: tuple[LocatedTextBlock, ...],
) -> None:
    raw_bytes = source.encode("utf-8")
    assert parsed.artifact_sha256 == sha256(raw_bytes).hexdigest()
    assert parsed.artifact_size_bytes == len(raw_bytes)
    for source_segment in parsed.source_segments:
        segment_span = source_segment.source_span
        segment_bytes = raw_bytes[segment_span.start_byte : segment_span.end_byte]
        assert source_segment.source_block_sha256 == sha256(segment_bytes).hexdigest()

    frames = parsed.frames
    assert len(frames) == len(blocks)
    for frame, block in zip(frames, blocks, strict=True):
        expected_text = block.text(source)
        span = frame.source_span
        assert span is not None
        assert frame.frame_content == expected_text
        assert source[span.start_char : span.end_char] == expected_text
        expected_bytes = expected_text.encode("utf-8")
        assert raw_bytes[span.start_byte : span.end_byte] == expected_bytes
        assert frame.source_block_sha256 == sha256(expected_bytes).hexdigest()


def _assert_basic_frame_parity(plain: Any, captured: Any) -> None:
    assert plain.frame_content == captured.frame_content
    assert plain.atoms == captured.atoms
    assert plain.charge == captured.charge
    assert plain.multiplicity == captured.multiplicity
    assert plain.coords.units == captured.coords.units
    assert plain.coords.m.tolist() == captured.coords.m.tolist()
    for field in (
        "keywords",
        "method",
        "basis_set",
        "functional",
        "model_chemistry",
        "task_requests",
        "resource_request",
    ):
        if hasattr(plain, field):
            assert getattr(plain, field) == getattr(captured, field)


def _sdf_record(smiles: str, label: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None
    mol_block = Chem.MolToMolBlock(mol).replace("\n", "\r\n")
    return f"{mol_block}>  <label>\r\n{label}\r\n\r\n$$$$\r\n"


@pytest.mark.parametrize(
    ("parser_cls", "source"),
    [
        (
            XYZFileParserMemory,
            "2\r\n\u7b2c\u4e00\u5e27\r\nH 0.0 0.0 0.0\r\nH 0.0 0.0 0.7\r\n"
            "1\r\nlast\r\nHe 1.0 2.0 3.0\r\n",
        ),
        (SMIFileParserMemory, "C methane\r\n\r\nO water\n"),
        (SDFFileParserMemory, _sdf_record("CO", "\u7532") + _sdf_record("N", "\u4e59")),
    ],
)
def test_coordinate_source_blocks_are_exact_and_only_last_keeps_true_index(
    parser_cls: type,
    source: str,
) -> None:
    parser = parser_cls(capture_source_evidence=True)
    segments = parser._locate_segments(source)
    assert segments is not None
    assert len(segments) == 1

    complete = parser.parse(source)
    plain = parser_cls().parse(source)
    assert len(plain.frames) == len(complete.frames)
    for plain_frame, captured_frame in zip(plain.frames, complete.frames, strict=True):
        _assert_basic_frame_parity(plain_frame, captured_frame)
    _assert_exact_frame_blocks(complete, source, segments[0].frames)
    assert [frame.segment_index for frame in complete.frames] == [0] * len(complete)
    assert [frame.segment_frame_index for frame in complete.frames] == list(range(len(complete)))
    assert [frame.file_frame_index for frame in complete.frames] == list(range(len(complete)))

    partial = parser_cls(
        capture_source_evidence=True,
        only_last_frame=True,
    ).parse(source)
    assert len(partial) == 1
    assert partial.source_complete is False
    assert partial[0].source_span == complete[-1].source_span
    assert partial[0].source_block_sha256 == complete[-1].source_block_sha256
    assert partial[0].segment_index == 0
    assert partial[0].segment_frame_index == len(complete) - 1
    assert partial[0].file_frame_index == len(complete) - 1


def test_sdf_located_records_with_dollar_delimiters_reach_frame_parser() -> None:
    source = _sdf_record("CO", "first") + _sdf_record("N", "second")
    parsed = SDFFileParserMemory(capture_source_evidence=True).parse(source)

    assert len(parsed) == 2
    assert [frame.atoms for frame in parsed.frames] == [[6, 8], [7]]
    assert all(frame.frame_content.endswith("$$$$\r\n") for frame in parsed.frames)


def _gjf_job(title: str, distance: str) -> str:
    return f"#p hf/3-21g\n\n{title}\n\n0 1\nH 0.0 0.0 0.0\nH 0.0 0.0 {distance}\n\n"


def test_gjf_link1_blocks_are_exact_and_only_last_keeps_segment_index() -> None:
    first_job = _gjf_job("\u7b2c\u4e00\u6b65", "0.7")
    source = f"{first_job}--Link1--\n{_gjf_job('second', '0.8')}"
    blocks = locate_gjf_link1_frames(source)

    complete = GJFFileParserMemory(capture_source_evidence=True).parse(source)
    plain = GJFFileParserMemory().parse(source)
    assert len(plain.frames) == len(complete.frames)
    for plain_frame, captured_frame in zip(plain.frames, complete.frames, strict=True):
        _assert_basic_frame_parity(plain_frame, captured_frame)
    _assert_exact_frame_blocks(complete, source, blocks)
    assert len(complete.source_segments) == 2
    assert [frame.segment_index for frame in complete.frames] == [0, 1]
    assert [frame.segment_frame_index for frame in complete.frames] == [0, 0]
    assert [frame.file_frame_index for frame in complete.frames] == [0, 1]

    partial = GJFFileParserMemory(
        capture_source_evidence=True,
        only_last_frame=True,
    ).parse(source)
    assert len(partial) == 1
    assert partial[0].source_span == complete[-1].source_span
    assert partial[0].segment_index == 1
    assert partial[0].segment_frame_index == 0
    assert partial[0].file_frame_index == 1


@pytest.mark.parametrize(
    "source",
    [
        f"{_gjf_job('first', '0.7')}--Link1--\n",
        (f"{_gjf_job('first', '0.7')}--Link1--\n\n--Link1--\n{_gjf_job('third', '0.9')}"),
        f"{_gjf_job('first', '0.7')}--Link1--\n".replace("\n", "\r\n"),
    ],
)
def test_gjf_rejects_link1_markers_with_empty_jobs(source: str) -> None:
    with pytest.raises(ValueError, match="non-empty Gaussian input jobs"):
        GJFFileParserMemory().parse(source)


def test_gjf_crlf_link1_blocks_keep_exact_frame_spans() -> None:
    source = (f"{_gjf_job('first', '0.7')}--Link1--\n{_gjf_job('second', '0.8')}").replace(
        "\n", "\r\n"
    )
    parser = GJFFileParserMemory(capture_source_evidence=True)

    segments = parser._locate_segments(source)
    parsed = parser.parse(source)
    blocks = tuple(segment.frames[0] for segment in segments)

    assert len(parsed) == 2
    _assert_exact_frame_blocks(parsed, source, blocks)
    assert [frame.segment_index for frame in parsed.frames] == [0, 1]


def test_gjf_include_is_rejected_until_multi_artifact_spans_are_supported(
    tmp_path: Path,
) -> None:
    included = tmp_path / "included.gjf"
    included.write_text(_gjf_job("included", "0.7"), encoding="utf-8")
    parent = tmp_path / "parent.gjf"
    parent.write_text("@included.gjf\n", encoding="utf-8")

    with pytest.raises(ParseError, match=GJF_INCLUDE_PROVENANCE_DIAGNOSTIC):
        GJFFileParserDisk(capture_source_evidence=True).parse(str(parent))


def test_orca_new_job_inside_geometry_is_not_a_delimiter() -> None:
    source = (
        "! SP\n* xyz 0 1\nH 0.0 0.0 0.0\n$new_job\n*\n$new_job\n! SP\n* xyz 0 1\nH 0.0 0.0 0.0\n*\n"
    )
    blocks = locate_orca_input_frames(source)

    assert len(blocks) == 2
    assert "$new_job" in blocks[0].text(source)
    assert "$new_job" not in blocks[1].text(source)
    separator = source.rindex("$new_job\n")
    assert [block.text(source) for block in blocks] == [
        source[:separator],
        source[separator + len("$new_job\n") :],
    ]


def test_orca_new_job_inside_percent_coords_is_not_a_delimiter() -> None:
    source = (
        "! SP\n"
        "%coords\n"
        "  CTyp xyz\n"
        "  Charge 0\n"
        "  Mult 1\n"
        "  coords\n"
        "    H 0.0 0.0 0.0\n"
        "    $new_job\n"
        "  end\n"
        "end\n"
        "$new_job\n"
        "! SP\n"
        "* xyz 0 1\n"
        "H 0.0 0.0 0.0\n"
        "*\n"
    )
    blocks = locate_orca_input_frames(source)

    assert len(blocks) == 2
    assert "$new_job" in blocks[0].text(source)
    separator = source.rindex("$new_job\n")
    assert [block.text(source) for block in blocks] == [
        source[:separator],
        source[separator + len("$new_job\n") :],
    ]


def test_orca_external_geometry_does_not_hide_following_new_job() -> None:
    source = "! SP\n* xyzfile 0 1 geometry.xyz\n$new_job\n! SP\n* xyz 0 1\nH 0.0 0.0 0.0\n*\n"
    blocks = locate_orca_input_frames(source)
    assert len(blocks) == 2

    complete = ORCAInpFileParserMemory(capture_source_evidence=True).parse(source)
    plain = ORCAInpFileParserMemory().parse(source)
    assert len(plain.frames) == len(complete.frames)
    for plain_frame, captured_frame in zip(plain.frames, complete.frames, strict=True):
        _assert_basic_frame_parity(plain_frame, captured_frame)
    _assert_exact_frame_blocks(complete, source, blocks)
    assert complete[0].geometry is not None
    assert complete[0].geometry.ctype == "xyzfile"
    assert complete[0].geometry.external_path == "geometry.xyz"

    partial = ORCAInpFileParserMemory(
        capture_source_evidence=True,
        only_last_frame=True,
    ).parse(source)
    assert len(partial) == 1
    assert partial[0].source_span == complete[-1].source_span
    assert partial[0].segment_index == 1
    assert partial[0].segment_frame_index == 0
    assert partial[0].file_frame_index == 1
