from pathlib import Path

import pytest

from molop.io.codec_exceptions import ParseError
from molop.io.logic.gaussian.input.parsers._gjf_file_extractors import (
    GJF_INCLUDE_PROVENANCE_DIAGNOSTIC,
)
from molop.io.logic.gaussian.input.parsers.GJFFileParser import (
    GJFFileParserDisk,
    GJFFileParserMemory,
)


def test_gjf_disk_parser_rejects_include_without_multi_artifact_spans(tmp_path: Path) -> None:
    included = tmp_path / "included.gjf"
    included.write_text(
        "%chk=inc.chk\n#p hf/3-21g\n\nincluded title\n\n0 1\nH 0.0 0.0 0.0\nH 0.0 0.0 0.7\n"
    )

    main = tmp_path / "main.gjf"
    main.write_text("@included.gjf\n")

    with pytest.raises(ParseError, match=GJF_INCLUDE_PROVENANCE_DIAGNOSTIC):
        GJFFileParserDisk().parse(str(main))


def test_gjf_memory_parser_rejects_include_without_multi_artifact_spans() -> None:
    with pytest.raises(ParseError, match=GJF_INCLUDE_PROVENANCE_DIAGNOSTIC):
        GJFFileParserMemory().parse("@other.gjf\n")


def test_gjf_disk_parser_rejects_include_before_following_recursive_sources(
    tmp_path: Path,
) -> None:
    a_file = tmp_path / "a.gjf"
    b_file = tmp_path / "b.gjf"
    a_file.write_text("@b.gjf\n")
    b_file.write_text("@a.gjf\n")

    with pytest.raises(ParseError, match=GJF_INCLUDE_PROVENANCE_DIAGNOSTIC):
        GJFFileParserDisk().parse(str(a_file))
