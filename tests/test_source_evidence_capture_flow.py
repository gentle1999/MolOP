from __future__ import annotations

import importlib
from collections.abc import Mapping
from hashlib import sha256
from pathlib import Path
from typing import Any, ClassVar, cast

import pytest

import molop.io as io_module
from molop.io.base_models.FrameParser import FrameParseContext
from molop.io.base_models.source import DecodedSource, LocatedSourceSegment, LocatedTextBlock
from molop.io.codec_types import ParseOptions, ParseResult, StructureLevel
from molop.io.codecs._shared.reader_helpers import ParserDiskReader
from molop.io.FileBatchModelDisk import FileBatchModelDisk
from molop.io.FileBatchParserDisk import FileBatchParserDisk, single_file_parser
from molop.io.logic.coords.frame_parsers.XYZFileFrameParser import (
    XYZFileFrameParserMemory,
)
from molop.io.logic.coords.parsers.XYZFileParser import XYZFileParserDisk, XYZFileParserMemory


file_batch_parser_module = importlib.import_module("molop.io.FileBatchParserDisk")

_XYZ_TEXT = "1\ncharge 1 multiplicity 2\nH 0.0 0.0 0.0\n"
_SECOND_XYZ_TEXT = "1\ncharge 0 multiplicity 1\nH 1.0 0.0 0.0\n"
_TWO_XYZ_TEXT = _XYZ_TEXT + _SECOND_XYZ_TEXT


class _FakeDiskFile:
    def __init__(self, file_path: str) -> None:
        self.file_path = file_path
        self.filename = Path(file_path).name
        self.file_format = Path(file_path).suffix
        self.detected_format_id: str | None = None

    def __len__(self) -> int:
        return 1

    def __getitem__(self, index: int) -> object:
        _ = index
        return object()

    def format_transform(self, *_args: Any, **_kwargs: Any) -> str:
        return ""

    def to_summary_series(self, **_kwargs: Any) -> object:
        return {}

    def release_file_content(self) -> None:
        return None


class _TrackingXYZFrameParser(XYZFileFrameParserMemory):
    seen_capture_flags: ClassVar[list[bool]] = []
    seen_structure_flags: ClassVar[list[bool]] = []
    seen_additional_data: ClassVar[list[Mapping[str, Any]]] = []

    def parse(self, block: str, *, additional_data: dict[str, Any] | None = None) -> Any:
        self.seen_capture_flags.append(self.capture_source_evidence)
        self.seen_structure_flags.append(self.only_extract_structure)
        return super().parse(block, additional_data=additional_data)

    def _parse_frame(self, block: str, *, context: FrameParseContext) -> Mapping[str, Any]:
        self.seen_additional_data.append(context.additional_data)
        return super()._parse_frame(block, context=context)


class _TrackingXYZFileParser(XYZFileParserMemory):
    pre_append_content: ClassVar[list[tuple[bool, bool]]] = []

    def _postprocess_parsed_frame(
        self,
        frame: Any,
        chem_file: Any,
        *,
        segment_index: int | None,
        segment_frame_index: int,
        segment_frame_count: int,
    ) -> None:
        _ = segment_index, segment_frame_index, segment_frame_count
        self.pre_append_content.append((bool(chem_file.file_content), bool(frame.frame_content)))

    def _source_frame_fields(self, frame: Any, **kwargs: Any) -> Mapping[str, Any]:
        return super()._source_frame_fields(frame, **kwargs) | {
            "coordinate_source": "xyz",
        }


class _InvalidXYZFileParser(XYZFileParserMemory):
    def _locate_segments(self, file_content: str) -> None:
        _ = file_content
        return None


class _UncoveredXYZFileParser(XYZFileParserMemory):
    def _locate_segments(self, file_content: str) -> tuple[LocatedSourceSegment, ...]:
        prefix = "unclaimed\n"
        return (
            LocatedSourceSegment(
                segment=LocatedTextBlock(0, len(file_content)),
                frames=(LocatedTextBlock(len(prefix), len(file_content)),),
            ),
        )


class _SegmentedXYZFileParser(XYZFileParserMemory):
    parsed_segment_contents: ClassVar[list[str]] = []

    def _locate_segments(self, file_content: str) -> tuple[LocatedSourceSegment, ...]:
        located = super()._locate_segments(file_content)
        return tuple(
            LocatedSourceSegment(segment=frame, frames=(frame,)) for frame in located[0].frames
        )

    def _parse_segment_metadata(
        self,
        segment_content: str,
        *,
        artifact_metadata: Mapping[str, Any],
    ) -> dict[str, Any] | None:
        self.parsed_segment_contents.append(segment_content)
        return super()._parse_segment_metadata(
            segment_content,
            artifact_metadata=artifact_metadata,
        )


def test_base_parser_capture_is_opt_in() -> None:
    parsed = XYZFileParserMemory().parse(_XYZ_TEXT)

    assert len(parsed.frames) == 1
    assert parsed.artifact_sha256 is None
    assert parsed.frames[0].source_span is None
    assert parsed.frames[0].file_frame_index is None
    assert "file_frame_index" not in parsed.frames[0].model_dump()
    assert "file_frame_index" not in parsed.frames[0].to_unitless_dump_with_unit_keys(
        exclude_none=True
    )


def test_memory_parser_accepts_preloaded_bytes_with_exact_source_identity() -> None:
    raw_bytes = _XYZ_TEXT.replace("\n", "\r\n").encode("utf-8")

    parsed = XYZFileParserMemory(capture_source_evidence=True).parse_bytes(raw_bytes)

    assert parsed.artifact_sha256 == sha256(raw_bytes).hexdigest()
    assert parsed.file_content == raw_bytes.decode("utf-8")
    assert parsed[0].source_span is not None
    assert parsed[0].source_span.end_byte == len(raw_bytes)


def test_disk_parser_accepts_decoded_source_without_reopening_path(tmp_path: Path) -> None:
    virtual_path = tmp_path / "preloaded.xyz"
    source = DecodedSource.from_bytes(_XYZ_TEXT.encode("utf-8"))

    parsed = XYZFileParserDisk(capture_source_evidence=True).parse_decoded_source(
        source,
        file_path=str(virtual_path),
    )

    assert parsed.file_path == str(virtual_path.resolve())
    assert parsed.artifact_sha256 == sha256(source.raw_bytes).hexdigest()


def test_locator_contract_has_no_splitter_fallback() -> None:
    with pytest.raises(TypeError, match=r"_locate_segments\(\).*non-empty sequence"):
        _InvalidXYZFileParser(capture_source_evidence=True).parse(_XYZ_TEXT)


@pytest.mark.parametrize("capture_source_evidence", [False, True])
def test_locator_contract_rejects_uncovered_non_whitespace_source(
    capture_source_evidence: bool,
) -> None:
    with pytest.raises(ValueError, match="cover every non-whitespace character"):
        _UncoveredXYZFileParser(capture_source_evidence=capture_source_evidence).parse(
            "unclaimed\n" + _XYZ_TEXT
        )


def test_base_parser_captures_before_release_and_propagates_to_frame_parser() -> None:
    _TrackingXYZFrameParser.seen_capture_flags.clear()
    _TrackingXYZFrameParser.seen_structure_flags.clear()
    _TrackingXYZFrameParser.seen_additional_data.clear()
    _TrackingXYZFileParser.pre_append_content.clear()
    parser = _TrackingXYZFileParser(
        capture_source_evidence=True,
        forced_charge=2,
        forced_multiplicity=3,
        only_extract_structure=True,
    )
    parser._frame_parser = _TrackingXYZFrameParser

    parsed = parser.parse(_XYZ_TEXT, release_file_content=True)

    expected_hash = sha256(_XYZ_TEXT.encode("utf-8")).hexdigest()
    assert _TrackingXYZFrameParser.seen_capture_flags == [True]
    assert _TrackingXYZFrameParser.seen_structure_flags == [True]
    assert _TrackingXYZFrameParser.seen_additional_data[0]["charge"] == 2
    assert _TrackingXYZFrameParser.seen_additional_data[0]["multiplicity"] == 3
    assert parsed.artifact_sha256 == expected_hash
    assert parsed.artifact_size_bytes == len(_XYZ_TEXT.encode("utf-8"))
    assert parsed.source_encoding == "utf-8"
    assert parsed.source_complete is True
    assert parsed.parse_completeness == "not_assessed"
    assert parsed.source_segments[0].parse_completeness == "not_assessed"
    assert parsed.source_segments[0]["source_block_sha256"] == expected_hash
    assert parsed.frames[0].source_block_sha256 == expected_hash
    assert parsed.frames[0].source_span is not None
    assert parsed.frames[0].file_frame_index == 0
    assert (
        parsed.frames[0].to_unitless_dump_with_unit_keys(exclude_none=True)["file_frame_index"] == 0
    )
    assert parsed.frames[0].source_span.end_byte == len(_XYZ_TEXT.encode("utf-8"))
    assert parsed.frames[0].coordinate_source == "xyz"
    assert parsed.frames[0].charge == 2
    assert parsed.frames[0].multiplicity == 3
    assert _TrackingXYZFileParser.pre_append_content == [(True, True)]
    assert parsed.file_content == ""
    assert parsed.frames[0].frame_content == ""


def test_only_last_filter_preserves_the_original_frame_index() -> None:
    parsed = XYZFileParserMemory(
        capture_source_evidence=True,
        only_last_frame=True,
    ).parse(_TWO_XYZ_TEXT)

    frame = parsed.frames[0]
    assert frame.segment_index == 0
    assert frame.segment_frame_index == 1
    assert frame.file_frame_index == 1
    assert frame.source_span is not None
    assert (
        frame.frame_content
        == _TWO_XYZ_TEXT[frame.source_span.start_char : frame.source_span.end_char]
    )


def test_only_last_parses_metadata_only_for_the_selected_segment() -> None:
    _SegmentedXYZFileParser.parsed_segment_contents.clear()

    parsed = _SegmentedXYZFileParser(
        capture_source_evidence=True,
        only_last_frame=True,
    ).parse(_TWO_XYZ_TEXT)

    assert len(parsed) == 1
    assert parsed[0].segment_index == 1
    assert parsed[0].file_frame_index == 1
    assert _SegmentedXYZFileParser.parsed_segment_contents == [_SECOND_XYZ_TEXT]


def test_single_file_parser_passes_capture_flag_to_reader() -> None:
    captured: dict[str, Any] = {}

    class RecordingReader:
        format_id = "xyz"

        def read(self, path: str | Path, **kwargs: Any) -> ParseResult[object]:
            captured.update(kwargs)
            return ParseResult(
                value=_FakeDiskFile(str(path)),
                level=StructureLevel.COORDS,
                detected_format=self.format_id,
            )

    parsed = single_file_parser(
        "/tmp/source.xyz",
        possible_readers=(RecordingReader(),),
        capture_source_evidence=True,
    )

    assert parsed is not None
    assert captured["capture_source_evidence"] is True


def test_batch_parser_includes_capture_flag_in_each_task(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    path = tmp_path / "source.xyz"
    path.write_text(_XYZ_TEXT, encoding="utf-8")
    captured: dict[str, Any] = {}

    class Reader:
        format_id = "xyz"

    def fake_single_file_parser(**task: Any) -> _FakeDiskFile:
        captured.update(task)
        return _FakeDiskFile(task["file_path"])

    monkeypatch.setattr(
        file_batch_parser_module.codec_registry,
        "select_reader",
        lambda _path, hint_format=None: (Reader(),),
    )
    monkeypatch.setattr(file_batch_parser_module, "single_file_parser", fake_single_file_parser)

    batch = FileBatchParserDisk(n_jobs=1).parse(
        [path],
        capture_source_evidence=True,
    )

    assert len(batch) == 1
    assert captured["capture_source_evidence"] is True
    assert captured["parse_options"] == ParseOptions(capture_source_evidence=True).resolved()


def test_parser_disk_reader_configures_parser_capture_flag() -> None:
    captured: dict[str, Any] = {}

    class Parser:
        def __init__(self, **kwargs: Any) -> None:
            captured["init"] = kwargs

        def parse(self, path: str, **kwargs: Any) -> object:
            captured["parse"] = {"path": path, **kwargs}
            return object()

    reader = ParserDiskReader(
        format_id="xyz",
        extensions=frozenset({".xyz"}),
        level=StructureLevel.COORDS,
        parser_cls=cast(Any, Parser),
        priority=1,
    )

    reader.read("/tmp/source.xyz", capture_source_evidence=True)

    assert captured["init"]["capture_source_evidence"] is True


def test_parser_disk_reader_accepts_immutable_parse_options() -> None:
    captured: dict[str, Any] = {}

    class Parser:
        def __init__(self, **kwargs: Any) -> None:
            captured["init"] = kwargs

        def parse(self, path: str, **kwargs: Any) -> object:
            captured["parse"] = {"path": path, **kwargs}
            return object()

    reader = ParserDiskReader(
        format_id="xyz",
        extensions=frozenset({".xyz"}),
        level=StructureLevel.COORDS,
        parser_cls=cast(Any, Parser),
        priority=1,
    )
    options = ParseOptions(
        total_charge=2,
        total_multiplicity=3,
        capture_source_evidence=True,
        release_file_content=False,
    )

    reader.read("/tmp/source.xyz", parse_options=options)

    assert captured["init"]["forced_charge"] == 2
    assert captured["init"]["forced_multiplicity"] == 3
    assert captured["init"]["capture_source_evidence"] is True
    assert captured["parse"]["release_file_content"] is False


def test_autoparser_passes_capture_flag_to_batch_parser(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    path = tmp_path / "source.xyz"
    path.write_text(_XYZ_TEXT, encoding="utf-8")
    captured: dict[str, Any] = {}

    class BatchParser:
        def __init__(self, n_jobs: int) -> None:
            captured["n_jobs"] = n_jobs

        def parse(self, file_paths: Any, **kwargs: Any) -> FileBatchModelDisk[Any]:
            captured["file_paths"] = list(file_paths)
            captured["kwargs"] = kwargs
            return FileBatchModelDisk()

    monkeypatch.setattr(io_module, "FileBatchParserDisk", BatchParser)

    io_module.AutoParser(path, n_jobs=1, capture_source_evidence=True)

    assert captured["kwargs"]["capture_source_evidence"] is True
