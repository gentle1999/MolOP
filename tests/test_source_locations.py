from __future__ import annotations

import codecs
from hashlib import sha256

import pytest

from molop.io.base_models.source import (
    DecodedSource,
    LocatedSourceSegment,
    LocatedTextBlock,
)


def test_utf8_bom_offsets_hashes_and_crlf_lines_are_exact() -> None:
    payload = b"alpha\r\nbeta"
    source = DecodedSource.from_bytes(codecs.BOM_UTF8 + payload, "utf-8-sig")
    complete = LocatedTextBlock(0, len(source.text))
    second_line = LocatedTextBlock(source.text.index("beta"), len(source.text))

    source.prepare_blocks((complete, second_line))
    complete_span = source.span(complete)
    second_line_span = source.span(second_line)

    assert complete_span.start_byte == len(codecs.BOM_UTF8)
    assert complete_span.end_byte == len(codecs.BOM_UTF8) + len(payload)
    assert complete_span.start_line == 1
    assert complete_span.end_line == 3
    assert source.block_sha256(complete_span) == sha256(payload).hexdigest()
    assert second_line_span.start_byte == len(codecs.BOM_UTF8) + len(b"alpha\r\n")
    assert second_line_span.end_byte == len(codecs.BOM_UTF8) + len(payload)
    assert second_line_span.start_line == 2
    assert second_line_span.end_line == 3


def test_decoded_source_builds_line_index_only_when_span_is_requested() -> None:
    source = DecodedSource.from_bytes(b"alpha\r\nbeta")
    block = LocatedTextBlock(0, len(source.text))

    assert source._line_starts_cache is None

    source.prepare_blocks((block,))
    assert source._line_starts_cache is None

    span = source.span(block)

    assert span.start_line == 1
    assert span.end_line == 3
    assert source._line_starts_cache == (0, len("alpha\r\n"))


@pytest.mark.parametrize(
    ("encoding", "bom", "payload_encoding"),
    [
        ("utf-16", codecs.BOM_UTF16_BE, "utf-16-be"),
        ("utf-32", codecs.BOM_UTF32_BE, "utf-32-be"),
    ],
)
def test_big_endian_bom_offsets_and_cr_only_lines_are_exact(
    encoding: str,
    bom: bytes,
    payload_encoding: str,
) -> None:
    text = "alpha\rbeta"
    payload = text.encode(payload_encoding)
    source = DecodedSource.from_bytes(bom + payload, encoding)
    complete = LocatedTextBlock(0, len(text))
    second_line = LocatedTextBlock(text.index("beta"), len(text))

    source.prepare_blocks((complete, second_line))
    complete_span = source.span(complete)
    second_line_span = source.span(second_line)

    assert complete_span.start_byte == len(bom)
    assert complete_span.end_byte == len(bom) + len(payload)
    assert source.block_sha256(complete_span) == sha256(payload).hexdigest()
    assert second_line_span.start_byte == len(bom) + len("alpha\r".encode(payload_encoding))
    assert second_line_span.end_byte == len(bom) + len(payload)
    assert second_line_span.start_line == 2
    assert second_line_span.end_line == 3


def test_decoded_source_uses_strict_decoding() -> None:
    with pytest.raises(UnicodeDecodeError):
        DecodedSource.from_bytes(b"\xff", "utf-8")


def test_located_source_segment_rejects_frames_outside_segment() -> None:
    with pytest.raises(ValueError, match="contained by their segment"):
        LocatedSourceSegment(
            segment=LocatedTextBlock(2, 8),
            frames=(LocatedTextBlock(1, 4),),
        )
