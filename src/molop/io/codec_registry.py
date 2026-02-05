"""Deterministic registry for IO reader/writer codecs.

Design goals:
- Deterministic precedence: sort by priority (desc), then registration order.
- Lazy activation: built-in codecs are registered only when first needed.
- Lazy instantiation: codecs are registered as factories and instantiated only when used.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass
from itertools import count
from pathlib import Path
from typing import Any, cast

from molop.io.codec_exceptions import (
    ConversionError,
    MissingOptionalDependencyError,
    UnsupportedFormatError,
)
from molop.io.codec_types import (
    GraphPolicy,
    ParseResult,
    ReaderCodec,
    StructureLevel,
    WriterCodec,
)


ReaderFactory = Callable[[], ReaderCodec]
WriterFactory = Callable[[], WriterCodec]

MAX_FALLBACK_READERS = 3
_OPENBABEL_FALLBACK_FORMAT_ID = "openbabel"


@dataclass(frozen=True, slots=True)
class _ReaderSpec:
    format_id: str
    extensions: tuple[str, ...]
    priority: int
    factory: ReaderFactory
    order: int


@dataclass(frozen=True, slots=True)
class _WriterSpec:
    format_id: str
    required_level: StructureLevel
    priority: int
    factory: WriterFactory
    order: int


@dataclass(slots=True)
class _LazyReaderCodec:
    format_id: str
    extensions: frozenset[str]
    priority: int
    _factory: ReaderFactory

    def read(self, path: str | Path, **kwargs: Any) -> ParseResult[object]:
        return self._factory().read(path, **kwargs)


class Registry:
    def __init__(self, *, autoload_defaults: bool = True) -> None:
        self._registration_counter = count()
        self._readers_by_format: dict[str, list[_ReaderSpec]] = {}
        self._readers_by_extension: dict[str, list[_ReaderSpec]] = {}
        self._fallback_readers: list[_ReaderSpec] = []
        self._writers_by_format: dict[str, list[_WriterSpec]] = {}
        self._openbabel_fallback_factory: ReaderFactory | None = None
        self._default_codecs_registered = False
        self._autoload_defaults = autoload_defaults

    def ensure_default_codecs_registered(self) -> None:
        """Idempotently register default codecs.

        The registry is intentionally empty until this is called.
        """

        if self._default_codecs_registered:
            return
        # Do NOT flip the flag before successful activation.
        from molop.io.codecs import catalog

        catalog.load_builtin_codecs(self)
        catalog.load_plugin_codecs(self)
        self._default_codecs_registered = True

    def register_reader(self, codec: ReaderCodec, *, fallback: bool = False) -> None:
        self.register_reader_factory(
            lambda: codec,
            format_id=codec.format_id,
            extensions=codec.extensions,
            priority=codec.priority,
            fallback=fallback,
        )

    def register_reader_factory(
        self,
        factory: ReaderFactory,
        *,
        format_id: str,
        extensions: Iterable[str] = (),
        priority: int = 0,
        fallback: bool = False,
    ) -> None:
        """Register a reader factory keyed by format id and extension."""

        normalized_format_id = _normalize_format_id(format_id)
        normalized_extensions = _normalize_extensions(extensions)
        spec = _ReaderSpec(
            format_id=normalized_format_id,
            extensions=normalized_extensions,
            priority=priority,
            factory=factory,
            order=next(self._registration_counter),
        )
        self._readers_by_format.setdefault(normalized_format_id, []).append(spec)
        self._readers_by_format[normalized_format_id].sort(key=_reader_sort_key)
        for extension in normalized_extensions:
            self._readers_by_extension.setdefault(extension, []).append(spec)
            self._readers_by_extension[extension].sort(key=_reader_sort_key)
        if fallback:
            self._fallback_readers.append(spec)
            self._fallback_readers.sort(key=_reader_sort_key)

    def register_writer(self, codec: WriterCodec) -> None:
        self.register_writer_factory(
            lambda: codec,
            format_id=codec.format_id,
            required_level=codec.required_level,
            priority=codec.priority,
        )

    def register_writer_factory(
        self,
        factory: WriterFactory,
        *,
        format_id: str,
        required_level: StructureLevel,
        priority: int = 0,
    ) -> None:
        normalized_format_id = _normalize_format_id(format_id)
        spec = _WriterSpec(
            format_id=normalized_format_id,
            required_level=required_level,
            priority=priority,
            factory=factory,
            order=next(self._registration_counter),
        )
        self._writers_by_format.setdefault(normalized_format_id, []).append(spec)
        self._writers_by_format[normalized_format_id].sort(key=_writer_sort_key)

    def register_openbabel_fallback(self, factory: ReaderFactory) -> None:
        self._openbabel_fallback_factory = factory

    def reader_factory(
        self,
        *,
        format_id: str,
        extensions: Iterable[str] = (),
        priority: int = 0,
        fallback: bool = False,
    ):
        """Decorator that registers a reader factory into this registry."""

        def decorator(factory: ReaderFactory) -> ReaderFactory:
            self.register_reader_factory(
                factory,
                format_id=format_id,
                extensions=extensions,
                priority=priority,
                fallback=fallback,
            )
            return factory

        return decorator

    def writer_factory(
        self,
        *,
        format_id: str,
        required_level: StructureLevel,
        priority: int = 0,
    ):
        """Decorator that registers a writer factory into this registry."""

        def decorator(factory: WriterFactory) -> WriterFactory:
            self.register_writer_factory(
                factory,
                format_id=format_id,
                required_level=required_level,
                priority=priority,
            )
            return factory

        return decorator

    def select_reader(
        self, path: str | Path, hint_format: str | None = None
    ) -> tuple[ReaderCodec, ...]:
        """Return reader codecs ordered by deterministic precedence."""

        if self._autoload_defaults:
            self.ensure_default_codecs_registered()
        candidates: list[_ReaderSpec] = []
        seen: set[int] = set()

        if hint_format:
            normalized_hint = _normalize_format_id(hint_format)
            _extend_unique(candidates, self._readers_by_format.get(normalized_hint, ()), seen)

        extension = _normalize_extension(Path(path).suffix)
        if extension:
            _extend_unique(candidates, self._readers_by_extension.get(extension, ()), seen)

        if self._fallback_readers:
            _extend_unique(candidates, self._fallback_readers[:MAX_FALLBACK_READERS], seen)

        lazy_codecs = [
            _LazyReaderCodec(
                format_id=spec.format_id,
                extensions=frozenset(spec.extensions),
                priority=spec.priority,
                _factory=spec.factory,
            )
            for spec in candidates
        ]

        if not lazy_codecs and self._openbabel_fallback_factory is not None:
            lazy_codecs.append(
                _LazyReaderCodec(
                    format_id=_OPENBABEL_FALLBACK_FORMAT_ID,
                    extensions=frozenset(),
                    priority=-10_000,
                    _factory=self._openbabel_fallback_factory,
                )
            )

        if not lazy_codecs:
            raise UnsupportedFormatError(f"No reader codecs registered for {path}.")
        return tuple(lazy_codecs)

    def write(
        self,
        format_id: str,
        value: object,
        *,
        graph_policy: GraphPolicy = "prefer",
        **kwargs: Any,
    ) -> object:
        """Write a value using registered writer codecs."""

        if self._autoload_defaults:
            self.ensure_default_codecs_registered()
        normalized_format_id = _normalize_format_id(format_id)
        writer_specs = self._writers_by_format.get(normalized_format_id, [])
        if not writer_specs:
            raise UnsupportedFormatError(f"No writer codecs registered for {format_id}.")

        raw_value, structure_level = _unwrap_parse_result(value)

        graph_writers = [
            spec for spec in writer_specs if spec.required_level == StructureLevel.GRAPH
        ]
        coords_writers = [
            spec for spec in writer_specs if spec.required_level == StructureLevel.COORDS
        ]

        if graph_policy == "coords":
            if not coords_writers:
                raise ConversionError(f"No coords-level writers registered for {format_id}.")
            return _write_with_specs(coords_writers, raw_value, **kwargs)

        if graph_policy == "strict":
            if not graph_writers:
                raise ConversionError(f"No graph-level writers registered for {format_id}.")
            graph_value = _coerce_graph_value(raw_value, structure_level)
            return _write_with_specs(graph_writers, graph_value, **kwargs)

        if graph_writers:
            try:
                graph_value = _coerce_graph_value(raw_value, structure_level)
                return _write_with_specs(graph_writers, graph_value, **kwargs)
            except ConversionError:
                pass
        if coords_writers:
            return _write_with_specs(coords_writers, raw_value, **kwargs)
        if graph_writers:
            raise ConversionError(f"Unable to upgrade value for graph writers of {format_id}.")
        raise UnsupportedFormatError(f"No compatible writers registered for {format_id}.")


default_registry = Registry(autoload_defaults=True)


def ensure_default_codecs_registered() -> None:
    default_registry.ensure_default_codecs_registered()


def register_reader(codec: ReaderCodec, *, fallback: bool = False) -> None:
    default_registry.register_reader(codec, fallback=fallback)


def register_reader_factory(
    factory: ReaderFactory,
    *,
    format_id: str,
    extensions: Iterable[str] = (),
    priority: int = 0,
    fallback: bool = False,
) -> None:
    default_registry.register_reader_factory(
        factory,
        format_id=format_id,
        extensions=extensions,
        priority=priority,
        fallback=fallback,
    )


def register_writer(codec: WriterCodec) -> None:
    default_registry.register_writer(codec)


def register_writer_factory(
    factory: WriterFactory,
    *,
    format_id: str,
    required_level: StructureLevel,
    priority: int = 0,
) -> None:
    default_registry.register_writer_factory(
        factory,
        format_id=format_id,
        required_level=required_level,
        priority=priority,
    )


def register_openbabel_fallback(factory: ReaderFactory) -> None:
    default_registry.register_openbabel_fallback(factory)


def reader_factory(
    *,
    format_id: str,
    extensions: Iterable[str] = (),
    priority: int = 0,
    fallback: bool = False,
):
    return default_registry.reader_factory(
        format_id=format_id,
        extensions=extensions,
        priority=priority,
        fallback=fallback,
    )


def writer_factory(
    *,
    format_id: str,
    required_level: StructureLevel,
    priority: int = 0,
):
    return default_registry.writer_factory(
        format_id=format_id,
        required_level=required_level,
        priority=priority,
    )


def select_reader(path: str | Path, hint_format: str | None = None) -> tuple[ReaderCodec, ...]:
    return default_registry.select_reader(path, hint_format=hint_format)


def write(
    format_id: str,
    value: object,
    *,
    graph_policy: GraphPolicy = "prefer",
    **kwargs: Any,
) -> object:
    return default_registry.write(format_id, value, graph_policy=graph_policy, **kwargs)


def upgrade_coords_to_graph(
    value: object,
    *,
    total_charge: int | None = None,
    total_radical_electrons: int | None = None,
) -> object:
    """Upgrade a coordinate-level value to a graph-level representation."""

    raw_value, _ = _unwrap_parse_result(value)
    rdmol = _try_get_attr(raw_value, "rdmol")
    if rdmol is not None:
        return rdmol
    omol = _try_get_attr(raw_value, "omol")
    if omol is not None:
        return omol

    if isinstance(raw_value, str):
        xyz_block = raw_value
    else:
        xyz_block_getter = getattr(raw_value, "to_XYZ_block", None)
        if callable(xyz_block_getter):
            xyz_block = xyz_block_getter()
            if not isinstance(xyz_block, str):
                raise ConversionError("XYZ block must be a string for graph upgrade.")
        else:
            raise ConversionError("No coordinate data available for graph upgrade.")

    charge, radical_electrons = _resolve_charge_and_radicals(
        raw_value,
        total_charge=total_charge,
        total_radical_electrons=total_radical_electrons,
    )

    from molop.structure import xyz_to_rdmol

    graph_value = xyz_to_rdmol(
        xyz_block,
        total_charge=charge,
        total_radical_electrons=radical_electrons,
    )
    if graph_value is None:
        raise ConversionError("Graph reconstruction failed from coordinates.")
    return graph_value


def _coerce_graph_value(value: object, level: StructureLevel | None) -> object:
    if level == StructureLevel.GRAPH:
        return value
    if hasattr(value, "frames") and hasattr(value, "model_dump"):
        return value
    return upgrade_coords_to_graph(value)


def _write_with_specs(
    specs: Sequence[_WriterSpec],
    value: object,
    **kwargs: Any,
) -> object:
    last_missing: MissingOptionalDependencyError | None = None
    for spec in specs:
        try:
            codec = spec.factory()
        except MissingOptionalDependencyError as exc:
            last_missing = exc
            continue
        return codec.write(value, **kwargs)
    if last_missing is not None:
        raise last_missing
    raise UnsupportedFormatError("No usable writer codecs available.")


def _unwrap_parse_result(value: object) -> tuple[object, StructureLevel | None]:
    if isinstance(value, ParseResult):
        return value.value, value.level
    return value, None


def _resolve_charge_and_radicals(
    value: object,
    *,
    total_charge: int | None,
    total_radical_electrons: int | None,
) -> tuple[int, int]:
    charge = total_charge
    if charge is None:
        charge_attr = _try_get_attr(value, "charge")
        if charge_attr is not None:
            try:
                charge = int(cast(Any, charge_attr))
            except (TypeError, ValueError):
                charge = None
    radicals = total_radical_electrons
    if radicals is None:
        multiplicity_attr = _try_get_attr(value, "multiplicity")
        if multiplicity_attr is not None:
            try:
                radicals = int(cast(Any, multiplicity_attr)) - 1
            except (TypeError, ValueError):
                radicals = None
    return charge or 0, radicals or 0


def _try_get_attr(value: object, name: str) -> object | None:
    try:
        return getattr(value, name)
    except Exception:
        return None


def _reader_sort_key(spec: _ReaderSpec) -> tuple[int, int]:
    return -spec.priority, spec.order


def _writer_sort_key(spec: _WriterSpec) -> tuple[int, int]:
    return -spec.priority, spec.order


def _extend_unique(target: list[_ReaderSpec], specs: Sequence[_ReaderSpec], seen: set[int]) -> None:
    for spec in specs:
        spec_id = id(spec)
        if spec_id in seen:
            continue
        seen.add(spec_id)
        target.append(spec)


def _normalize_format_id(format_id: str) -> str:
    normalized = format_id.strip().lower()
    if not normalized:
        raise ValueError("format_id must be non-empty")
    return normalized


def _normalize_extension(extension: str) -> str:
    normalized = extension.strip().lower()
    if not normalized:
        return ""
    if not normalized.startswith("."):
        normalized = f".{normalized}"
    return normalized


def _normalize_extensions(extensions: Iterable[str]) -> tuple[str, ...]:
    normalized = {_normalize_extension(ext) for ext in extensions if ext}
    return tuple(sorted(ext for ext in normalized if ext))


__all__ = [
    "MAX_FALLBACK_READERS",
    "Registry",
    "default_registry",
    "ensure_default_codecs_registered",
    "reader_factory",
    "register_openbabel_fallback",
    "register_reader",
    "register_reader_factory",
    "register_writer",
    "register_writer_factory",
    "select_reader",
    "upgrade_coords_to_graph",
    "write",
    "writer_factory",
]
