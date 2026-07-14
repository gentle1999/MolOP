"""Format-independent source decoding and location values."""

from __future__ import annotations

import codecs
import json
from bisect import bisect_right
from collections.abc import Mapping
from dataclasses import dataclass, field
from enum import Enum
from hashlib import sha256
from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field, model_validator
from typing_extensions import Self


class SourceSpan(BaseModel):
    """Half-open byte, character, and line offsets in one decoded source."""

    model_config = ConfigDict(extra="forbid", frozen=True)

    start_byte: int = Field(ge=0)
    end_byte: int = Field(ge=0)
    start_char: int = Field(ge=0)
    end_char: int = Field(ge=0)
    start_line: int = Field(ge=1)
    end_line: int = Field(ge=1)

    @model_validator(mode="after")
    def validate_half_open_ranges(self) -> Self:
        if self.end_byte <= self.start_byte:
            raise ValueError("end_byte must be greater than start_byte")
        if self.end_char <= self.start_char:
            raise ValueError("end_char must be greater than start_char")
        if self.end_line <= self.start_line:
            raise ValueError("end_line must be greater than start_line")
        return self


class ParsePresence(str, Enum):
    """Why an optional parsed value is present or missing."""

    NOT_REQUESTED = "not_requested"
    ABSENT_IN_SOURCE = "absent_in_source"
    PARSED = "parsed"
    PARSE_FAILED = "parse_failed"
    UNSUPPORTED = "unsupported"


class ParseCompleteness(str, Enum):
    """Whether the requested scientific parsing work was fully assessed."""

    NOT_ASSESSED = "not_assessed"
    COMPLETE = "complete"
    PARTIAL = "partial"


class ParseDiagnostic(BaseModel):
    """Stable structured diagnostic emitted by the normal parser lifecycle."""

    model_config = ConfigDict(extra="forbid", frozen=True)

    code: str = Field(min_length=1)
    severity: Literal["info", "warning", "error"]
    scope: str = Field(min_length=1)
    message: str = Field(min_length=1)
    segment_index: int | None = Field(default=None, ge=0)
    segment_frame_index: int | None = Field(default=None, ge=0)
    field: str | None = Field(default=None, min_length=1)


def canonical_json_sha256(value: Mapping[str, Any]) -> str:
    """Hash a JSON mapping with the calculation-export canonical encoding."""

    canonical_json = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    )
    return sha256(canonical_json.encode("utf-8")).hexdigest()


class ParserProvenance(BaseModel):
    """Version and effective-configuration snapshot for one file parse."""

    model_config = ConfigDict(extra="forbid", frozen=True)

    parser_id: str = Field(min_length=1)
    parser_version: str = Field(min_length=1)
    molop_version: str = Field(min_length=1)
    molgr_version: str = Field(min_length=1)
    rdkit_version: str = Field(min_length=1)
    effective_config: dict[str, Any]
    effective_config_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")

    @model_validator(mode="after")
    def validate_effective_config_sha256(self) -> Self:
        try:
            expected_hash = canonical_json_sha256(self.effective_config)
        except (TypeError, ValueError) as exc:
            raise ValueError(
                "effective_config must contain only finite canonical JSON values"
            ) from exc
        if self.effective_config_sha256 != expected_hash:
            raise ValueError("effective_config_sha256 does not match effective_config")
        return self


class SourceSegmentEvidence(BaseModel):
    """Located segment identity plus evidence whose scope is the whole segment."""

    model_config = ConfigDict(extra="forbid", frozen=True)

    segment_index: int = Field(ge=0)
    source_span: SourceSpan
    source_block_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    frame_count: int = Field(ge=0)
    captured_frame_indices: list[int] = Field(default_factory=list)
    qm_software: str | None = Field(default=None, min_length=1)
    qm_software_version: str | None = Field(default=None, min_length=1)
    protocol: dict[str, Any] | None = None
    task_requests: list[dict[str, Any]] = Field(default_factory=list)
    task_types: list[str] = Field(default_factory=list)
    termination_status: bool | None = None
    scf_status: bool | None = None
    parse_presence: dict[str, ParsePresence] = Field(default_factory=dict)
    parse_completeness: ParseCompleteness = ParseCompleteness.NOT_ASSESSED
    diagnostics: list[ParseDiagnostic] = Field(default_factory=list)

    @model_validator(mode="after")
    def validate_frame_indices(self) -> Self:
        if self.captured_frame_indices != sorted(set(self.captured_frame_indices)):
            raise ValueError("captured_frame_indices must be sorted and unique")
        if any(index >= self.frame_count for index in self.captured_frame_indices):
            raise ValueError("captured_frame_indices must be smaller than frame_count")
        return self

    def __getitem__(self, key: str) -> Any:
        """Keep legacy mapping-style reads while exposing a typed model."""

        return getattr(self, key)


@dataclass(frozen=True, slots=True)
class LocatedTextBlock:
    """A non-empty half-open character range in one decoded artifact."""

    start_char: int
    end_char: int

    def __post_init__(self) -> None:
        if self.start_char < 0 or self.end_char <= self.start_char:
            raise ValueError("located text blocks must be non-empty half-open ranges")

    def text(self, source_text: str) -> str:
        if self.end_char > len(source_text):
            raise ValueError("located text block exceeds its source text")
        return source_text[self.start_char : self.end_char]


@dataclass(frozen=True, slots=True)
class LocatedSourceSegment:
    """One exact source segment and the parser frame blocks it contains."""

    segment: LocatedTextBlock
    frames: tuple[LocatedTextBlock, ...]

    def __post_init__(self) -> None:
        for frame in self.frames:
            if frame.start_char < self.segment.start_char or frame.end_char > self.segment.end_char:
                raise ValueError("source frame blocks must be contained by their segment")


@dataclass(slots=True)
class DecodedSource:
    """One strictly decoded artifact with exact character-to-byte mapping."""

    raw_bytes: bytes
    text: str
    encoding: str
    offset_encoding: str
    content_start_byte: int
    byte_offsets: dict[int, int] = field(default_factory=dict)
    _line_starts_cache: tuple[int, ...] | None = field(default=None, init=False, repr=False)

    @staticmethod
    def _offset_codec(raw_bytes: bytes, encoding: str) -> tuple[str, int]:
        """Return a stateless codec and content offset for BOM-aware encodings."""

        if encoding == "utf-8-sig":
            offset = len(codecs.BOM_UTF8) if raw_bytes.startswith(codecs.BOM_UTF8) else 0
            return "utf-8", offset
        if encoding == "utf-16":
            if raw_bytes.startswith(codecs.BOM_UTF16_LE):
                return "utf-16-le", len(codecs.BOM_UTF16_LE)
            if raw_bytes.startswith(codecs.BOM_UTF16_BE):
                return "utf-16-be", len(codecs.BOM_UTF16_BE)
        if encoding == "utf-32":
            if raw_bytes.startswith(codecs.BOM_UTF32_LE):
                return "utf-32-le", len(codecs.BOM_UTF32_LE)
            if raw_bytes.startswith(codecs.BOM_UTF32_BE):
                return "utf-32-be", len(codecs.BOM_UTF32_BE)
        return encoding, 0

    @staticmethod
    def _line_starts(text: str) -> tuple[int, ...]:
        starts = [0]
        index = 0
        while index < len(text):
            char = text[index]
            if char == "\r":
                if index + 1 < len(text) and text[index + 1] == "\n":
                    index += 1
                starts.append(index + 1)
            elif char == "\n":
                starts.append(index + 1)
            index += 1
        return tuple(starts)

    @classmethod
    def from_bytes(cls, raw_bytes: bytes, encoding: str = "utf-8") -> DecodedSource:
        canonical_encoding = codecs.lookup(encoding).name
        text = raw_bytes.decode(canonical_encoding, errors="strict")
        offset_encoding, content_start_byte = cls._offset_codec(
            raw_bytes,
            canonical_encoding,
        )
        return cls(
            raw_bytes=raw_bytes,
            text=text,
            encoding=canonical_encoding,
            offset_encoding=offset_encoding,
            content_start_byte=content_start_byte,
            byte_offsets={0: content_start_byte},
        )

    @classmethod
    def from_text(cls, text: str, encoding: str = "utf-8") -> DecodedSource:
        canonical_encoding = codecs.lookup(encoding).name
        return cls.from_bytes(text.encode(canonical_encoding), canonical_encoding)

    def prepare_blocks(self, blocks: tuple[LocatedTextBlock, ...]) -> None:
        positions = sorted(
            {
                position
                for block in blocks
                for position in (block.start_char, block.end_char)
                if position not in self.byte_offsets
            }
        )
        if not positions:
            return

        encoder = codecs.getincrementalencoder(self.offset_encoding)(errors="strict")
        previous_char = 0
        previous_byte = self.content_start_byte
        for position in positions:
            if not 0 <= position <= len(self.text):
                raise ValueError("source block boundary is outside the decoded artifact")
            encoded_chunk = encoder.encode(self.text[previous_char:position], final=False)
            end_byte = previous_byte + len(encoded_chunk)
            if self.raw_bytes[previous_byte:end_byte] != encoded_chunk:
                raise ValueError(
                    f"encoding {self.encoding!r} cannot reproduce exact source byte offsets"
                )
            self.byte_offsets[position] = end_byte
            previous_char = position
            previous_byte = end_byte

    def _byte_offset(self, char_offset: int) -> int:
        if not 0 <= char_offset <= len(self.text):
            raise ValueError("character offset is outside the decoded artifact")
        if char_offset not in self.byte_offsets:
            self.prepare_blocks((LocatedTextBlock(0, char_offset),))
        return self.byte_offsets[char_offset]

    def _line_at_char(self, char_offset: int) -> int:
        line_starts = self._line_starts_cache
        if line_starts is None:
            line_starts = self._line_starts(self.text)
            self._line_starts_cache = line_starts
        return bisect_right(line_starts, char_offset)

    def span(self, block: LocatedTextBlock) -> SourceSpan:
        if not 0 <= block.start_char < block.end_char <= len(self.text):
            raise ValueError("source block must be a non-empty range within the decoded artifact")
        return SourceSpan(
            start_byte=self._byte_offset(block.start_char),
            end_byte=self._byte_offset(block.end_char),
            start_char=block.start_char,
            end_char=block.end_char,
            start_line=self._line_at_char(block.start_char),
            end_line=self._line_at_char(block.end_char - 1) + 1,
        )

    def block_sha256(self, span: SourceSpan) -> str:
        return sha256(self.raw_bytes[span.start_byte : span.end_byte]).hexdigest()


__all__ = [
    "DecodedSource",
    "LocatedSourceSegment",
    "LocatedTextBlock",
    "ParseCompleteness",
    "ParseDiagnostic",
    "ParsePresence",
    "ParserProvenance",
    "SourceSegmentEvidence",
    "SourceSpan",
    "canonical_json_sha256",
]
