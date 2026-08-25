from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

from molop.io.base_models.FileParser import BaseFileParserDisk, BaseFileParserMemory
from molop.io.base_models.source import normalize_parser_line_endings
from molop.io.codec_exceptions import FormatMismatchError, UnsupportedFormatError
from molop.io.codec_types import ParseOptions, ParseResult, ReaderCodec, StructureLevel


def extensions_for_parser(
    parser_cls: type[BaseFileParserDisk[Any, Any, Any]],
) -> tuple[str, ...]:
    return parser_cls.allowed_formats


def _parse_options_from_kwargs(kwargs: dict[str, Any]) -> ParseOptions:
    parse_options = kwargs.pop("parse_options", None)
    if parse_options is not None and not isinstance(parse_options, ParseOptions):
        raise TypeError("parse_options must be a ParseOptions instance")
    return (
        parse_options
        or ParseOptions(
            total_charge=kwargs.pop("total_charge", None),
            total_multiplicity=kwargs.pop("total_multiplicity", None),
            only_extract_structure=kwargs.pop("only_extract_structure", False),
            only_last_frame=kwargs.pop("only_last_frame", False),
            capture_source_evidence=kwargs.pop("capture_source_evidence", False),
            source_encoding=kwargs.pop("source_encoding", "utf-8"),
            release_file_content=kwargs.pop("release_file_content", True),
        )
    ).resolved()


MemoryParserClass = type[BaseFileParserMemory[Any, Any, Any]]


@dataclass
class ParserMemoryReader:
    """Adapt a memory parser to the registry's loaded-source reader contract."""

    format_id: str
    extensions: frozenset[str]
    level: StructureLevel
    parser_cls: MemoryParserClass
    priority: int

    def _parser(self, options: ParseOptions) -> BaseFileParserMemory[Any, Any, Any]:
        return self.parser_cls(
            forced_charge=options.total_charge,
            forced_multiplicity=options.total_multiplicity,
            only_extract_structure=options.only_extract_structure,
            only_last_frame=options.only_last_frame,
            capture_source_evidence=options.capture_source_evidence,
            source_encoding=options.source_encoding,
            parse_options=options,
        )

    def probe_text(self, text: str) -> bool:
        quick_check = getattr(self.parser_cls, "_quick_check_file_format", None)
        if not callable(quick_check):
            return True
        try:
            quick_check(normalize_parser_line_endings(text))
        except FormatMismatchError:
            return False
        return True

    def read_text(self, text: str, **kwargs: Any) -> ParseResult[object]:
        options = _parse_options_from_kwargs(kwargs)
        try:
            value = self._parser(options).parse(
                text,
                total_charge=options.total_charge,
                total_multiplicity=options.total_multiplicity,
                release_file_content=options.release_file_content,
            )
        except ValueError as exc:
            if "no locator-provided source segments" in str(exc):
                raise FormatMismatchError(str(exc)) from exc
            raise
        return ParseResult(value=value, level=self.level, detected_format=self.format_id)

    def read_bytes(self, raw_bytes: bytes, **kwargs: Any) -> ParseResult[object]:
        options = _parse_options_from_kwargs(kwargs)
        try:
            value = self._parser(options).parse_bytes(
                raw_bytes,
                total_charge=options.total_charge,
                total_multiplicity=options.total_multiplicity,
                release_file_content=options.release_file_content,
            )
        except ValueError as exc:
            if "no locator-provided source segments" in str(exc):
                raise FormatMismatchError(str(exc)) from exc
            raise
        return ParseResult(value=value, level=self.level, detected_format=self.format_id)


@dataclass
class ParserDiskReader:
    format_id: str
    extensions: frozenset[str]
    level: StructureLevel
    parser_cls: type[BaseFileParserDisk[Any, Any, Any]]
    priority: int
    memory_parser_cls: MemoryParserClass | None = None

    def probe_file_format(self, path: str | Path) -> bool:
        return self.parser_cls.probe_file_format(path)

    def read(self, path: str | Path, **kwargs: Any) -> ParseResult[object]:
        parse_options = kwargs.pop("parse_options", None)
        if parse_options is not None and not isinstance(parse_options, ParseOptions):
            raise TypeError("parse_options must be a ParseOptions instance")
        options = (
            parse_options
            or ParseOptions(
                total_charge=kwargs.pop("total_charge", None),
                total_multiplicity=kwargs.pop("total_multiplicity", None),
                only_extract_structure=kwargs.pop("only_extract_structure", False),
                only_last_frame=kwargs.pop("only_last_frame", False),
                capture_source_evidence=kwargs.pop("capture_source_evidence", False),
                source_encoding=kwargs.pop("source_encoding", "utf-8"),
                release_file_content=kwargs.pop("release_file_content", True),
            )
        ).resolved()
        parser = self.parser_cls(
            forced_charge=options.total_charge,
            forced_multiplicity=options.total_multiplicity,
            only_extract_structure=options.only_extract_structure,
            only_last_frame=options.only_last_frame,
            capture_source_evidence=options.capture_source_evidence,
            source_encoding=options.source_encoding,
            parse_options=options,
        )
        value = parser.parse(
            str(path),
            total_charge=options.total_charge,
            total_multiplicity=options.total_multiplicity,
            release_file_content=options.release_file_content,
        )
        return ParseResult(value=value, level=self.level, detected_format=self.format_id)

    def _memory_reader(self) -> ParserMemoryReader:
        if self.memory_parser_cls is None:
            raise UnsupportedFormatError(
                f"Reader {self.format_id!r} has no in-memory parser implementation."
            )
        return ParserMemoryReader(
            format_id=self.format_id,
            extensions=self.extensions,
            level=self.level,
            parser_cls=self.memory_parser_cls,
            priority=self.priority,
        )

    def probe_text(self, text: str) -> bool:
        if self.memory_parser_cls is None:
            return True
        return self._memory_reader().probe_text(text)

    def read_text(self, text: str, **kwargs: Any) -> ParseResult[object]:
        return self._memory_reader().read_text(text, **kwargs)

    def read_bytes(self, raw_bytes: bytes, **kwargs: Any) -> ParseResult[object]:
        return self._memory_reader().read_bytes(raw_bytes, **kwargs)


__all__ = [
    "ParserDiskReader",
    "ParserMemoryReader",
    "ReaderCodec",
    "StructureLevel",
    "extensions_for_parser",
]
