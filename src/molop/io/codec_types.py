from __future__ import annotations

from dataclasses import dataclass, replace
from enum import Enum
from pathlib import Path
from typing import Generic, Literal, Protocol, TypeVar


class StructureLevel(str, Enum):
    GRAPH = "graph"
    COORDS = "coords"


GraphPolicy = Literal["prefer", "strict", "coords"]
WriterDomain = Literal["file", "frame"]


@dataclass(frozen=True, slots=True)
class ParseWarning:
    code: str
    message: str


@dataclass(frozen=True, slots=True)
class ConversionWarning:
    code: str
    message: str


@dataclass(frozen=True, slots=True)
class ParseOptions:
    """Immutable options shared by batch dispatch and reader codecs."""

    total_charge: int | None = None
    total_multiplicity: int | None = None
    only_extract_structure: bool = False
    only_last_frame: bool = False
    capture_source_evidence: bool = False
    source_encoding: str = "utf-8"
    release_file_content: bool = True
    force_unit_transform: bool | None = None
    graph_reconstruction_backend: Literal["cpp", "python"] | None = None
    make_dative_bonds: bool | None = None

    def resolved(self) -> ParseOptions:
        """Snapshot process-global defaults that affect parse-time model behavior."""

        from molop.config import molopconfig

        return replace(
            self,
            force_unit_transform=(
                molopconfig.force_unit_transform
                if self.force_unit_transform is None
                else self.force_unit_transform
            ),
            graph_reconstruction_backend=(
                molopconfig.graph_reconstruction_backend
                if self.graph_reconstruction_backend is None
                else self.graph_reconstruction_backend
            ),
            make_dative_bonds=(
                molopconfig.make_dative_bonds
                if self.make_dative_bonds is None
                else self.make_dative_bonds
            ),
        )


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
