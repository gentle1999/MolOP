from pathlib import Path

import pytest

from molop.io.logic.qminput_parsers.GJFFileParser import GJFFileParserDisk, GJFFileParserMemory


def test_gjf_disk_parser_expands_at_include(tmp_path: Path) -> None:
    included = tmp_path / "included.gjf"
    included.write_text(
        "%chk=inc.chk\n#p hf/3-21g\n\nincluded title\n\n0 1\nH 0.0 0.0 0.0\nH 0.0 0.0 0.7\n"
    )

    main = tmp_path / "main.gjf"
    main.write_text("@included.gjf\n")

    parser = GJFFileParserDisk()
    parsed = parser.parse(str(main))
    assert len(parsed) == 1
    assert parsed[0].title_card.title_card.strip() == "included title"


def test_gjf_memory_parser_rejects_include_without_path_context() -> None:
    parser = GJFFileParserMemory()
    with pytest.raises(ValueError, match="without source file path context"):
        parser.parse("@other.gjf\n")


def test_gjf_disk_parser_detects_recursive_include_cycle(tmp_path: Path) -> None:
    a_file = tmp_path / "a.gjf"
    b_file = tmp_path / "b.gjf"
    a_file.write_text("@b.gjf\n")
    b_file.write_text("@a.gjf\n")

    parser = GJFFileParserDisk()
    with pytest.raises(ValueError, match="include cycle"):
        parser.parse(str(a_file))
