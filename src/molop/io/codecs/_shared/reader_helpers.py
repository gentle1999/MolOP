from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

from molop.io.base_models.FileParser import BaseFileParserDisk
from molop.io.codec_types import ParseOptions, ParseResult, ReaderCodec, StructureLevel


def extensions_for_parser(
    parser_cls: type[BaseFileParserDisk[Any, Any, Any]],
) -> tuple[str, ...]:
    return parser_cls.allowed_formats


@dataclass
class ParserDiskReader:
    format_id: str
    extensions: frozenset[str]
    level: StructureLevel
    parser_cls: type[BaseFileParserDisk[Any, Any, Any]]
    priority: int

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


__all__ = ["ParserDiskReader", "ReaderCodec", "StructureLevel", "extensions_for_parser"]
