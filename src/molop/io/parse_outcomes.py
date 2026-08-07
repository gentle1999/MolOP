"""Structured outcomes for batch parsing.

The regular ``AutoParser`` API intentionally remains a collection of successful
files.  These immutable records provide an opt-in diagnostic layer without
serializing exception objects or coupling the batch model to parser failures.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Generic, Literal, TypeVar

from molop.io.codec_types import ParseWarning


ParseStatus = Literal[
    "ok",
    "empty",
    "missing",
    "skipped",
    "unsupported",
    "mismatch",
    "error",
]
T = TypeVar("T")
B = TypeVar("B")


@dataclass(frozen=True, slots=True)
class ParseFailure:
    """Serializable description of one failed parse attempt."""

    kind: str
    message: str
    reader_format: str | None = None
    exception_type: str | None = None


@dataclass(frozen=True, slots=True)
class FileParseOutcome(Generic[T]):
    """Result of attempting to parse one input path."""

    file_path: str
    status: ParseStatus
    value: T | None = None
    warnings: tuple[ParseWarning, ...] = ()
    detected_format: str | None = None
    failure: ParseFailure | None = None
    input_index: int | None = None

    @property
    def succeeded(self) -> bool:
        return self.status == "ok"


@dataclass(frozen=True, slots=True)
class BatchParseResult(Generic[B, T]):
    """A successful batch plus one outcome for every supplied input path."""

    batch: B
    outcomes: tuple[FileParseOutcome[T], ...]

    @property
    def failures(self) -> tuple[FileParseOutcome[T], ...]:
        return tuple(outcome for outcome in self.outcomes if not outcome.succeeded)

    @property
    def succeeded(self) -> tuple[FileParseOutcome[T], ...]:
        return tuple(outcome for outcome in self.outcomes if outcome.succeeded)


__all__ = [
    "BatchParseResult",
    "FileParseOutcome",
    "ParseFailure",
    "ParseStatus",
]
