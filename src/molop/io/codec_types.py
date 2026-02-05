from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Generic, Literal, Protocol, TypeVar


class StructureLevel(str, Enum):
    GRAPH = "graph"
    COORDS = "coords"


GraphPolicy = Literal["prefer", "strict", "coords"]


@dataclass(frozen=True, slots=True)
class ParseWarning:
    code: str
    message: str


@dataclass(frozen=True, slots=True)
class ConversionWarning:
    code: str
    message: str


T = TypeVar("T")


@dataclass(frozen=True, slots=True)
class ParseResult(Generic[T]):
    value: T
    level: StructureLevel
    warnings: tuple[ParseWarning, ...] = ()
    detected_format: str | None = None


class ReaderCodec(Protocol):
    @property
    def format_id(self) -> str: ...

    @property
    def extensions(self) -> frozenset[str]: ...

    @property
    def priority(self) -> int: ...

    def read(self, path: str | Path, **kwargs) -> ParseResult[object]: ...


class WriterCodec(Protocol):
    format_id: str
    priority: int
    required_level: StructureLevel

    def write(self, value: object, **kwargs) -> object: ...
