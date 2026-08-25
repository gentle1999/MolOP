from __future__ import annotations

from hashlib import sha256
from io import BytesIO, StringIO
from pathlib import Path

import pytest

from molop import AutoBytesParser, AutoMemoryParser, AutoTextParser
from molop.io.base_models.ChemFile import BaseChemFile
from molop.io.codec_exceptions import FormatMismatchError


XYZ_TEXT = "1\nwater\nH 0.0 0.0 0.0\n"


def test_auto_text_parser_returns_a_memory_file_model() -> None:
    parsed = AutoTextParser(XYZ_TEXT, parser_detection="xyz", release_file_content=False)

    assert isinstance(parsed, BaseChemFile)
    assert type(parsed).__name__ == "XYZFileMemory"
    assert parsed.source_format == "xyz"
    assert parsed.file_content == XYZ_TEXT


def test_auto_bytes_parser_preserves_loaded_source_identity() -> None:
    raw_bytes = XYZ_TEXT.replace("\n", "\r\n").encode("utf-8")

    parsed = AutoBytesParser(
        raw_bytes,
        parser_detection="xyz",
        capture_source_evidence=True,
        release_file_content=False,
    )

    assert parsed.source_format == "xyz"
    assert parsed.artifact_sha256 == sha256(raw_bytes).hexdigest()
    assert parsed.artifact_size_bytes == len(raw_bytes)
    assert parsed.file_content == raw_bytes.decode("utf-8")
    assert parsed[0].source_block_sha256 == sha256(raw_bytes).hexdigest()
    assert parsed[0].source_span.end_byte == len(raw_bytes)


def test_auto_memory_parser_accepts_text_and_binary_streams() -> None:
    text_parsed = AutoMemoryParser(StringIO(XYZ_TEXT), parser_detection="xyz")
    bytes_parsed = AutoMemoryParser(BytesIO(XYZ_TEXT.encode()), parser_detection="xyz")

    assert text_parsed.source_format == "xyz"
    assert bytes_parsed.source_format == "xyz"


def test_auto_memory_parser_can_detect_simple_coordinate_formats_without_a_path() -> None:
    parsed = AutoTextParser(XYZ_TEXT, release_file_content=False)

    assert parsed.source_format == "xyz"
    assert len(parsed) == 1


def test_auto_memory_parser_can_detect_output_formats_without_a_path() -> None:
    source = Path("tests/test_files/g16log/H2O.log").read_text(encoding="utf-8")

    parsed = AutoTextParser(source, only_extract_structure=True)

    assert parsed.source_format == "g16log"
    assert len(parsed) > 0


def test_auto_memory_parser_reports_format_mismatch() -> None:
    with pytest.raises(FormatMismatchError):
        AutoTextParser("not a supported chemical source", parser_detection="xyz")
