from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

from molop.io.base_models.FileParser import BaseFileParserDisk
from molop.io.codec_types import ParseResult, ReaderCodec, StructureLevel


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

    def read(self, path: str | Path, **kwargs: Any) -> ParseResult[object]:
        total_charge = kwargs.pop("total_charge", None)
        total_multiplicity = kwargs.pop("total_multiplicity", None)
        only_extract_structure = kwargs.pop("only_extract_structure", False)
        only_last_frame = kwargs.pop("only_last_frame", False)
        release_file_content = kwargs.pop("release_file_content", True)
        parser = self.parser_cls(
            forced_charge=total_charge,
            forced_multiplicity=total_multiplicity,
            only_extract_structure=only_extract_structure,
            only_last_frame=only_last_frame,
        )
        value = parser.parse(
            str(path),
            total_charge=total_charge,
            total_multiplicity=total_multiplicity,
            release_file_content=release_file_content,
        )
        return ParseResult(value=value, level=self.level, detected_format=self.format_id)


__all__ = ["ParserDiskReader", "ReaderCodec", "StructureLevel", "extensions_for_parser"]
