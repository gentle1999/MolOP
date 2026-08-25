"""Public automatic parsers for already-loaded text and byte sources."""

from __future__ import annotations

from typing import Any, BinaryIO, TextIO, TypeAlias, cast, overload

from molop.config import moloplogger
from molop.io import codec_registry
from molop.io.base_models.ChemFile import BaseChemFile
from molop.io.codec_exceptions import FormatMismatchError, ParseError, UnsupportedFormatError
from molop.io.codec_types import MemoryReaderCodec, ParseOptions, ParseResult


TextSource: TypeAlias = str | TextIO
BytesSource: TypeAlias = bytes | bytearray | memoryview | BinaryIO
MemorySource: TypeAlias = TextSource | BytesSource


def _resolved_options(
    *,
    total_charge: int | None,
    total_multiplicity: int | None,
    only_extract_structure: bool,
    only_last_frame: bool,
    capture_source_evidence: bool,
    source_encoding: str,
    release_file_content: bool,
    parse_options: ParseOptions | None,
) -> ParseOptions:
    if parse_options is not None and not isinstance(parse_options, ParseOptions):
        raise TypeError("parse_options must be a ParseOptions instance")
    return (
        parse_options
        or ParseOptions(
            total_charge=total_charge,
            total_multiplicity=total_multiplicity,
            only_extract_structure=only_extract_structure,
            only_last_frame=only_last_frame,
            capture_source_evidence=capture_source_evidence,
            source_encoding=source_encoding,
            release_file_content=release_file_content,
        )
    ).resolved()


def _coerce_text(source: TextSource) -> str:
    if isinstance(source, str):
        return source
    read = getattr(source, "read", None)
    if not callable(read):
        raise TypeError("text source must be str or a text stream")
    value = read()
    if not isinstance(value, str):
        raise TypeError("text stream must return str")
    return value


def _coerce_bytes(source: BytesSource) -> bytes:
    if isinstance(source, bytes):
        return source
    if isinstance(source, (bytearray, memoryview)):
        return bytes(source)
    read = getattr(source, "read", None)
    if not callable(read):
        raise TypeError("byte source must be bytes-like or a binary stream")
    value = read()
    if not isinstance(value, (bytes, bytearray, memoryview)):
        raise TypeError("binary stream must return bytes-like data")
    return bytes(value)


def _reader_name(reader: MemoryReaderCodec) -> str:
    return getattr(reader, "format_id", reader.__class__.__name__)


def _filter_readers_by_text(
    text: str,
    readers: tuple[MemoryReaderCodec, ...],
) -> tuple[MemoryReaderCodec, ...]:
    matching: list[MemoryReaderCodec] = []
    for reader in readers:
        probe = getattr(reader, "probe_text", None)
        if not callable(probe):
            matching.append(reader)
            continue
        try:
            if probe(text):
                matching.append(reader)
        except Exception as exc:
            moloplogger.debug(
                "Reader %s text probe failed; keeping it as a candidate. %s",
                _reader_name(reader),
                exc,
            )
            matching.append(reader)
    return tuple(matching) or readers


def _parse_with_readers(
    source: str | bytes,
    *,
    source_kind: str,
    readers: tuple[MemoryReaderCodec, ...],
    options: ParseOptions,
    allow_parse_fallback: bool,
) -> BaseChemFile[Any]:
    attempted = 0
    unsupported = 0
    mismatches: list[str] = []
    for reader in readers:
        reader_name = _reader_name(reader)
        read = getattr(reader, f"read_{source_kind}", None)
        if not callable(read):
            unsupported += 1
            continue
        attempted += 1
        try:
            result = cast(
                ParseResult[object],
                read(source, parse_options=options),
            )
        except UnsupportedFormatError:
            unsupported += 1
            continue
        except FormatMismatchError as exc:
            mismatches.append(f"{reader_name}: {exc}")
            continue
        except StopIteration as exc:
            # Open Babel-backed coordinate readers use StopIteration when a
            # text block is not a valid record for their format.
            mismatches.append(f"{reader_name}: {exc}")
            continue
        except ValueError as exc:
            # A simple coordinate locator can recognize a broad text block,
            # then reject it while validating frame coverage. Treat that
            # locator-only failure as a candidate mismatch so output readers
            # can still be tried when the source has no filename extension.
            if "Located source frames must cover" in str(exc):
                mismatches.append(f"{reader_name}: {exc}")
                continue
            if allow_parse_fallback:
                mismatches.append(f"{reader_name}: {exc}")
                continue
            raise ParseError(
                f"Failed to parse loaded {source_kind} with {reader_name}: {exc}"
            ) from exc
        except Exception as exc:
            if allow_parse_fallback:
                mismatches.append(f"{reader_name}: {exc}")
                continue
            raise ParseError(
                f"Failed to parse loaded {source_kind} with {reader_name}: {exc}"
            ) from exc

        value = result.value
        if not isinstance(value, BaseChemFile):
            raise ParseError(
                f"Reader {reader_name} returned unexpected value type: {type(value).__name__}"
            )
        if len(value) == 0:
            raise ParseError(f"Reader {reader_name} returned a file model with no frames.")
        if getattr(value, "source_format", None) is None and result.detected_format:
            value.source_format = result.detected_format
        return value

    if attempted == 0 and unsupported:
        raise UnsupportedFormatError(
            "No selected reader supports parsing the loaded source in memory."
        )
    details = "; ".join(mismatches)
    suffix = f" {details}" if details else ""
    raise FormatMismatchError(f"No reader accepted the loaded {source_kind} format.{suffix}")


def _parse_loaded_source(
    source: str | bytes,
    *,
    source_kind: str,
    parser_detection: str,
    options: ParseOptions,
) -> BaseChemFile[Any]:
    hint_format = None if parser_detection == "auto" else parser_detection
    readers = codec_registry.select_memory_reader(hint_format=hint_format)
    if parser_detection == "auto":
        probe_text = source if isinstance(source, str) else source.decode(options.source_encoding)
        readers = _filter_readers_by_text(probe_text, readers)
    return _parse_with_readers(
        source,
        source_kind=source_kind,
        readers=readers,
        options=options,
        allow_parse_fallback=parser_detection == "auto",
    )


@overload
def AutoMemoryParser(
    source: TextSource,
    *,
    total_charge: int | None = None,
    total_multiplicity: int | None = None,
    only_extract_structure: bool = False,
    only_last_frame: bool = False,
    capture_source_evidence: bool = False,
    source_encoding: str = "utf-8",
    release_file_content: bool = True,
    parser_detection: str = "auto",
    parse_options: ParseOptions | None = None,
) -> BaseChemFile[Any]: ...


@overload
def AutoMemoryParser(
    source: BytesSource,
    *,
    total_charge: int | None = None,
    total_multiplicity: int | None = None,
    only_extract_structure: bool = False,
    only_last_frame: bool = False,
    capture_source_evidence: bool = False,
    source_encoding: str = "utf-8",
    release_file_content: bool = True,
    parser_detection: str = "auto",
    parse_options: ParseOptions | None = None,
) -> BaseChemFile[Any]: ...


def AutoMemoryParser(
    source: MemorySource,
    *,
    total_charge: int | None = None,
    total_multiplicity: int | None = None,
    only_extract_structure: bool = False,
    only_last_frame: bool = False,
    capture_source_evidence: bool = False,
    source_encoding: str = "utf-8",
    release_file_content: bool = True,
    parser_detection: str = "auto",
    parse_options: ParseOptions | None = None,
) -> BaseChemFile[Any]:
    """Parse already-loaded text or bytes with automatic format detection.

    A string is treated as source text, never as a path. Use :func:`AutoFileParser`
    for path input. Binary and text streams are consumed from their current position.
    """

    options = _resolved_options(
        total_charge=total_charge,
        total_multiplicity=total_multiplicity,
        only_extract_structure=only_extract_structure,
        only_last_frame=only_last_frame,
        capture_source_evidence=capture_source_evidence,
        source_encoding=source_encoding,
        release_file_content=release_file_content,
        parse_options=parse_options,
    )
    if isinstance(source, str):
        return _parse_loaded_source(
            source,
            source_kind="text",
            parser_detection=parser_detection,
            options=options,
        )
    if isinstance(source, (bytes, bytearray, memoryview)):
        return _parse_loaded_source(
            bytes(source),
            source_kind="bytes",
            parser_detection=parser_detection,
            options=options,
        )
    read = getattr(source, "read", None)
    if callable(read):
        value = read()
        if isinstance(value, str):
            return _parse_loaded_source(
                value,
                source_kind="text",
                parser_detection=parser_detection,
                options=options,
            )
        if isinstance(value, (bytes, bytearray, memoryview)):
            return _parse_loaded_source(
                bytes(value),
                source_kind="bytes",
                parser_detection=parser_detection,
                options=options,
            )
        raise TypeError("source stream must return str or bytes-like data")
    return _parse_loaded_source(
        _coerce_bytes(cast(BytesSource, source)),
        source_kind="bytes",
        parser_detection=parser_detection,
        options=options,
    )


def AutoTextParser(
    text: TextSource,
    *,
    total_charge: int | None = None,
    total_multiplicity: int | None = None,
    only_extract_structure: bool = False,
    only_last_frame: bool = False,
    capture_source_evidence: bool = False,
    source_encoding: str = "utf-8",
    release_file_content: bool = True,
    parser_detection: str = "auto",
    parse_options: ParseOptions | None = None,
) -> BaseChemFile[Any]:
    """Parse text already held by the caller without creating a temporary file."""

    return AutoMemoryParser(
        _coerce_text(text),
        total_charge=total_charge,
        total_multiplicity=total_multiplicity,
        only_extract_structure=only_extract_structure,
        only_last_frame=only_last_frame,
        capture_source_evidence=capture_source_evidence,
        source_encoding=source_encoding,
        release_file_content=release_file_content,
        parser_detection=parser_detection,
        parse_options=parse_options,
    )


def AutoBytesParser(
    raw_bytes: BytesSource,
    *,
    total_charge: int | None = None,
    total_multiplicity: int | None = None,
    only_extract_structure: bool = False,
    only_last_frame: bool = False,
    capture_source_evidence: bool = False,
    source_encoding: str = "utf-8",
    release_file_content: bool = True,
    parser_detection: str = "auto",
    parse_options: ParseOptions | None = None,
) -> BaseChemFile[Any]:
    """Parse bytes already held by the caller without a filesystem round trip."""

    return AutoMemoryParser(
        _coerce_bytes(raw_bytes),
        total_charge=total_charge,
        total_multiplicity=total_multiplicity,
        only_extract_structure=only_extract_structure,
        only_last_frame=only_last_frame,
        capture_source_evidence=capture_source_evidence,
        source_encoding=source_encoding,
        release_file_content=release_file_content,
        parser_detection=parser_detection,
        parse_options=parse_options,
    )


AutoParserMemory = AutoMemoryParser


__all__ = [
    "AutoBytesParser",
    "AutoMemoryParser",
    "AutoParserMemory",
    "AutoTextParser",
    "BytesSource",
    "MemorySource",
    "TextSource",
]
