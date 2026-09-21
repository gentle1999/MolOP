from __future__ import annotations

import json
from pathlib import Path

from molop.io.logic.gaussian.log.parsers._g16log_archive_tail import _normalize_archive_text
from molop.io.logic.gaussian.log.parsers.G16LogFileParser import G16LogFileParserDisk


def test_unrecoverable_archive_byte_is_replaced_with_valid_unicode() -> None:
    assert _normalize_archive_text("bad byte: \udcff") == "bad byte: \ufffd"


def test_archive_comment_rejoins_split_utf8_and_is_json_safe() -> None:
    path = Path("tests/test_files/g16log/anie202116024-sup-0001-misc_information__71.log")
    parsed = G16LogFileParserDisk(
        source_decode_errors="surrogateescape",
        capture_source_evidence=True,
    ).parse(path)

    final_frame_comments = parsed.frames[-1].comments
    assert any(comment.text.endswith("TS6a-2\u2019") for comment in final_frame_comments)

    comments = [comment for frame in parsed.frames for comment in frame.comments]
    comments.extend(parsed.comments)
    assert all(
        not any(0xD800 <= ord(char) <= 0xDFFF for char in value)
        for comment in comments
        for value in (comment.text, comment.raw or "")
    )
    # PostgreSQL JSONB requires valid Unicode even when strings are JSON-escaped.
    json.dumps(
        [comment.model_dump(mode="json") for comment in comments], ensure_ascii=False
    ).encode("utf-8")

    assert parsed.source_complete is True
    assert parsed.file_content.encode("utf-8", errors="surrogateescape") == path.read_bytes()
