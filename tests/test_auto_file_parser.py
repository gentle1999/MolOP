from __future__ import annotations

from pathlib import Path

import pytest

import molop.io as io_module
from molop import AutoFileParser
from molop.io.base_models.ChemFile import BaseChemFile
from molop.io.codec_exceptions import FormatMismatchError
from molop.io.FileBatchModelDisk import FileBatchModelDisk


def test_auto_file_parser_returns_a_file_model_with_automatic_detection(tmp_path: Path) -> None:
    source = tmp_path / "water.xyz"
    source.write_text("1\nwater\nH 0.0 0.0 0.0\n", encoding="utf-8")

    parsed = AutoFileParser(source)

    assert isinstance(parsed, BaseChemFile)
    assert not isinstance(parsed, FileBatchModelDisk)
    assert parsed.file_path == str(source.resolve())
    assert parsed.detected_format_id == "xyz"
    assert len(parsed) == 1


def test_auto_file_parser_does_not_construct_the_batch_parser(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    source = tmp_path / "water.xyz"
    source.write_text("1\nwater\nH 0.0 0.0 0.0\n", encoding="utf-8")

    class UnexpectedBatchParser:
        def __init__(self, *_args: object, **_kwargs: object) -> None:
            raise AssertionError("single-file parsing must not construct a batch parser")

    monkeypatch.setattr(io_module, "FileBatchParserDisk", UnexpectedBatchParser)

    parsed = io_module.AutoFileParser(source)
    assert parsed.detected_format_id == "xyz"


def test_auto_file_parser_reports_a_format_mismatch(tmp_path: Path) -> None:
    source = tmp_path / "invalid.xyz"
    source.write_text("this is not an xyz file\n", encoding="utf-8")

    with pytest.raises(FormatMismatchError):
        AutoFileParser(source)


def test_auto_file_parser_checks_content_for_ambiguous_output_suffix() -> None:
    source = Path("tests/test_files/g16irc/irc.out").resolve()

    parsed = AutoFileParser(source, only_last_frame=True)

    assert parsed.detected_format_id == "g16log"
    assert parsed.file_path == str(source)
