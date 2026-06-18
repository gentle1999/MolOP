"""Deterministic registry for IO reader/writer codecs.

Design goals:
- Deterministic precedence: sort by priority (desc), then registration order.
- Lazy activation: built-in codecs are registered only when first needed.
- Lazy instantiation: codecs are registered as factories and instantiated only when used.
"""

from __future__ import annotations

import inspect
import re
from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass
from itertools import count
from pathlib import Path
from types import UnionType
from typing import Any, Literal, Union, cast, get_args, get_origin, get_type_hints

from molgr.interface import xyz_to_rdmol

from molop.config import molopconfig
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
    WriterDomain,
)


ReaderFactory = Callable[[], ReaderCodec]
WriterFactory = Callable[[], WriterCodec]

MAX_FALLBACK_READERS = 3
_OPENBABEL_FALLBACK_FORMAT_ID = "openbabel"
_PARAMETER_DOC_RE = re.compile(r"^([A-Za-z_]\w*)\s*(?:\([^)]*\))?\s*:\s*(.*)$")
_NUMPY_PARAMETER_DOC_RE = re.compile(r"^([A-Za-z_]\w*)\s*:\s*(.*)$")
_DEFAULT_SENTENCE_RE = re.compile(r"\s*Defaults? to .*$", re.IGNORECASE)
_OPTION_HELP_OVERRIDES = {
    "add_gjf_connectivity": "Append Gaussian connectivity data after the coordinate block.",
    "additional_sections": "Raw Gaussian extra-input sections appended after the molecule block.",
    "blocks": "Additional ORCA block syntax appended to the input preamble.",
    "chk": "Append a Gaussian %chk Link0 directive for this render.",
    "coords_type": "Preferred Gaussian coordinate representation.",
    "engine": "Rendering backend used to generate the target format.",
    "keywords": "Input keywords used in the generated quantum-chemistry input.",
    "link0_commands": "Gaussian Link0 resource and job-control directives.",
    "maxcore": "ORCA per-core memory limit in MB.",
    "molecule_specifications": "Gaussian charge, multiplicity, and coordinate block override.",
    "nprocs": "Number of CPU cores requested in the generated input.",
    "old_chk": "Append a Gaussian %oldchk Link0 directive for this render.",
    "parsed_additional_sections": "Structured Gaussian extra-input sections to render.",
    "resources_raw": "Raw ORCA resource/control block appended to the input preamble.",
    "route_section": "Gaussian route section keyword line.",
    "title_card": "Gaussian title section override.",
    "use_raw_geometry": "Preserve the original ORCA geometry block when possible.",
}


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
    domain: WriterDomain
    default_graph_policy: GraphPolicy
    priority: int
    factory: WriterFactory
    order: int


@dataclass(frozen=True, slots=True)
class WriterOptionSpec:
    name: str
    option: str
    value_candidates: tuple[str, ...] = ()
    takes_value: bool = True
    supports_no: bool = False
    help: str = ""


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
            default_graph_policy=getattr(codec, "default_graph_policy", None),
            priority=codec.priority,
        )

    def register_writer_factory(
        self,
        factory: WriterFactory,
        *,
        format_id: str,
        required_level: StructureLevel,
        domain: WriterDomain = "file",
        default_graph_policy: GraphPolicy | None = None,
        priority: int = 0,
    ) -> None:
        normalized_format_id = _normalize_format_id(format_id)
        spec = _WriterSpec(
            format_id=normalized_format_id,
            required_level=required_level,
            domain=domain,
            default_graph_policy=default_graph_policy
            or _default_graph_policy_for_level(required_level),
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
        domain: WriterDomain = "file",
        default_graph_policy: GraphPolicy | None = None,
        priority: int = 0,
    ):
        """Decorator that registers a writer factory into this registry."""

        def decorator(factory: WriterFactory) -> WriterFactory:
            self.register_writer_factory(
                factory,
                format_id=format_id,
                required_level=required_level,
                domain=domain,
                default_graph_policy=default_graph_policy,
                priority=priority,
            )
            return factory

        return decorator

    def get_supported_writer_formats(self) -> list[str]:
        """Return a sorted list of all registered writer format IDs."""
        if self._autoload_defaults:
            self.ensure_default_codecs_registered()
        return sorted(self._writers_by_format.keys())

    def get_writer_option_specs(self, format_id: str) -> list[WriterOptionSpec]:
        """Return completion metadata for writer-specific format options."""
        if self._autoload_defaults:
            self.ensure_default_codecs_registered()
        normalized_format_id = _normalize_format_id(format_id)
        writer_specs = self._writers_by_format.get(normalized_format_id, [])
        result: dict[str, WriterOptionSpec] = {}

        for spec in writer_specs:
            for option_spec in _writer_option_specs_from_factory(spec.factory):
                result.setdefault(option_spec.name, option_spec)

        return sorted(result.values(), key=lambda item: item.option)

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
        graph_policy: GraphPolicy | None = None,
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
        file_writer_specs = [spec for spec in writer_specs if spec.domain == "file"]
        if not file_writer_specs:
            raise UnsupportedFormatError(f"No file writer codecs registered for {format_id}.")

        return _write_with_domain_specs(
            file_writer_specs,
            raw_value,
            structure_level=structure_level,
            graph_policy=_resolve_graph_policy(file_writer_specs, graph_policy),
            **kwargs,
        )

    def write_frame(
        self,
        format_id: str,
        value: object,
        *,
        graph_policy: GraphPolicy | None = None,
        **kwargs: Any,
    ) -> object:
        if self._autoload_defaults:
            self.ensure_default_codecs_registered()
        normalized_format_id = _normalize_format_id(format_id)
        writer_specs = self._writers_by_format.get(normalized_format_id, [])
        if not writer_specs:
            raise UnsupportedFormatError(f"No writer codecs registered for {format_id}.")

        raw_value, structure_level = _unwrap_parse_result(value)
        frame_writer_specs = [spec for spec in writer_specs if spec.domain == "frame"]
        if not frame_writer_specs:
            raise UnsupportedFormatError(f"No frame writer codecs registered for {format_id}.")

        return _write_with_domain_specs(
            frame_writer_specs,
            raw_value,
            structure_level=structure_level,
            graph_policy=_resolve_graph_policy(frame_writer_specs, graph_policy),
            **kwargs,
        )


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
    domain: WriterDomain = "file",
    default_graph_policy: GraphPolicy | None = None,
    priority: int = 0,
) -> None:
    default_registry.register_writer_factory(
        factory,
        format_id=format_id,
        required_level=required_level,
        domain=domain,
        default_graph_policy=default_graph_policy,
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
    domain: WriterDomain = "file",
    default_graph_policy: GraphPolicy | None = None,
    priority: int = 0,
):
    return default_registry.writer_factory(
        format_id=format_id,
        required_level=required_level,
        domain=domain,
        default_graph_policy=default_graph_policy,
        priority=priority,
    )


def get_supported_writer_formats() -> list[str]:
    return default_registry.get_supported_writer_formats()


def get_writer_option_specs(format_id: str) -> list[WriterOptionSpec]:
    return default_registry.get_writer_option_specs(format_id)


def select_reader(path: str | Path, hint_format: str | None = None) -> tuple[ReaderCodec, ...]:
    return default_registry.select_reader(path, hint_format=hint_format)


def write(
    format_id: str,
    value: object,
    *,
    graph_policy: GraphPolicy | None = None,
    **kwargs: Any,
) -> object:
    return default_registry.write(format_id, value, graph_policy=graph_policy, **kwargs)


def write_frame(
    format_id: str,
    value: object,
    *,
    graph_policy: GraphPolicy | None = None,
    **kwargs: Any,
) -> object:
    return default_registry.write_frame(format_id, value, graph_policy=graph_policy, **kwargs)


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

    graph_value = xyz_to_rdmol(
        xyz_block,
        total_charge=charge,
        spin_multiplicity=radical_electrons + 1,
        backend=molopconfig.graph_reconstruction_backend,
        make_dative_bonds=molopconfig.make_dative_bonds,
    )
    if graph_value is None:
        raise ConversionError("Graph reconstruction failed from coordinates.")
    return graph_value


def _coerce_graph_value(value: object, level: StructureLevel | None) -> object:
    if level == StructureLevel.GRAPH:
        return value
    if _is_structured_render_value(value):
        return value
    return upgrade_coords_to_graph(value)


def _is_structured_render_value(value: object) -> bool:
    return hasattr(value, "model_dump") and (
        hasattr(value, "frames") or (hasattr(value, "frame_id") and hasattr(value, "_render"))
    )


def _write_with_domain_specs(
    specs: Sequence[_WriterSpec],
    value: object,
    *,
    structure_level: StructureLevel | None,
    graph_policy: GraphPolicy,
    **kwargs: Any,
) -> object:
    graph_writers = [spec for spec in specs if spec.required_level == StructureLevel.GRAPH]
    coords_writers = [spec for spec in specs if spec.required_level == StructureLevel.COORDS]

    if graph_policy == "coords":
        if not coords_writers:
            raise ConversionError("No coords-level writers registered.")
        return _write_with_specs(coords_writers, value, **kwargs)

    if graph_policy == "strict":
        if not graph_writers:
            raise ConversionError("No graph-level writers registered.")
        graph_value = _coerce_graph_value(value, structure_level)
        return _write_with_specs(graph_writers, graph_value, **kwargs)

    if graph_writers:
        try:
            graph_value = _coerce_graph_value(value, structure_level)
            return _write_with_specs(graph_writers, graph_value, **kwargs)
        except ConversionError:
            pass
    if coords_writers:
        return _write_with_specs(coords_writers, value, **kwargs)
    if graph_writers:
        raise ConversionError("Unable to upgrade value for graph writers.")
    raise UnsupportedFormatError("No compatible writers registered.")


def _default_graph_policy_for_level(required_level: StructureLevel) -> GraphPolicy:
    if required_level == StructureLevel.GRAPH:
        return "prefer"
    return "coords"


def _resolve_graph_policy(
    specs: Sequence[_WriterSpec], graph_policy: GraphPolicy | None
) -> GraphPolicy:
    if graph_policy is not None:
        return graph_policy
    return specs[0].default_graph_policy


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


def _writer_option_specs_from_factory(factory: WriterFactory) -> list[WriterOptionSpec]:
    try:
        writer = factory()
    except MissingOptionalDependencyError:
        return []

    specs: dict[str, WriterOptionSpec] = {
        "graph_policy": WriterOptionSpec(
            name="graph_policy",
            option="--graph-policy",
            value_candidates=("prefer", "strict", "coords"),
            help="Molecular graph policy used before rendering.",
        )
    }
    _collect_signature_options(specs, getattr(writer, "write", None))

    for attr_name in ("file_cls", "frame_cls"):
        cls = getattr(writer, attr_name, None)
        if cls is None:
            continue
        _collect_signature_options(specs, getattr(cls, "_render", None))

    return list(specs.values())


def _collect_signature_options(
    target: dict[str, WriterOptionSpec],
    callable_obj: Callable[..., object] | None,
) -> None:
    if callable_obj is None:
        return
    try:
        signature = inspect.signature(callable_obj)
    except (TypeError, ValueError):
        return
    try:
        type_hints = get_type_hints(callable_obj)
    except (NameError, TypeError):
        type_hints = {}

    parameter_docs = _extract_parameter_docs(callable_obj)
    for parameter in signature.parameters.values():
        if parameter.name in {
            "self",
            "value",
            "kwargs",
            "file_path",
            "frameID",
            "embed_in_one_file",
        }:
            continue
        if parameter.kind in {
            inspect.Parameter.VAR_KEYWORD,
            inspect.Parameter.VAR_POSITIONAL,
            inspect.Parameter.POSITIONAL_ONLY,
        }:
            continue
        annotation = type_hints.get(parameter.name, parameter.annotation)
        option = f"--{parameter.name.replace('_', '-')}"
        target.setdefault(
            parameter.name,
            WriterOptionSpec(
                name=parameter.name,
                option=option,
                value_candidates=_literal_value_candidates(annotation),
                takes_value=not _annotation_is_bool_only(annotation),
                supports_no=_annotation_accepts_bool(annotation),
                help=_writer_option_help(parameter, parameter_docs.get(parameter.name, "")),
            ),
        )


def _literal_value_candidates(annotation: object) -> tuple[str, ...]:
    origin = get_origin(annotation)
    if origin is Literal:
        return tuple(str(value) for value in get_args(annotation) if value is not None)
    if origin in {UnionType, Union}:
        values: list[str] = []
        for arg in get_args(annotation):
            values.extend(_literal_value_candidates(arg))
        return tuple(dict.fromkeys(values))
    return ()


def _annotation_is_bool_only(annotation: object) -> bool:
    if annotation is bool:
        return True
    origin = get_origin(annotation)
    if origin is Literal:
        return set(get_args(annotation)) <= {True, False}
    if origin in {UnionType, Union}:
        return set(get_args(annotation)) <= {bool, type(None)}
    return False


def _annotation_accepts_bool(annotation: object) -> bool:
    if annotation is bool:
        return True
    origin = get_origin(annotation)
    if origin is Literal:
        return bool(set(get_args(annotation)) & {True, False})
    if origin in {UnionType, Union}:
        return bool in get_args(annotation)
    return False


def _extract_parameter_docs(callable_obj: Callable[..., object]) -> dict[str, str]:
    doc = inspect.getdoc(callable_obj)
    if not doc:
        return {}

    docs: dict[str, str] = {}
    lines = doc.splitlines()
    for index, line in enumerate(lines):
        stripped = line.strip()
        match = _PARAMETER_DOC_RE.match(stripped) or _NUMPY_PARAMETER_DOC_RE.match(stripped)
        if match is None:
            continue
        name = match.group(1)
        description_parts = [match.group(2).strip()]
        for follow_line in lines[index + 1 :]:
            follow_stripped = follow_line.strip()
            if not follow_stripped:
                if any(description_parts):
                    break
                continue
            if (
                _PARAMETER_DOC_RE.match(follow_stripped)
                or _NUMPY_PARAMETER_DOC_RE.match(follow_stripped)
            ):
                break
            if follow_line.startswith((" ", "\t")):
                description_parts.append(follow_stripped)
                continue
            break
        description = _clean_option_help(" ".join(part for part in description_parts if part))
        if description:
            docs[name] = description
    return docs


def _writer_option_help(parameter: inspect.Parameter, doc_help: str = "") -> str:
    description = _OPTION_HELP_OVERRIDES.get(parameter.name, doc_help)
    default_help = _writer_option_default_help(parameter)
    if description and default_help:
        return f"{description} {default_help}"
    return description or default_help


def _writer_option_default_help(parameter: inspect.Parameter) -> str:
    if parameter.default is inspect.Parameter.empty:
        return ""
    return f"Default: {parameter.default!r}."


def _clean_option_help(value: str) -> str:
    compact = " ".join(value.split())
    compact = _DEFAULT_SENTENCE_RE.sub("", compact).strip()
    if not compact:
        return ""
    return compact[0].upper() + compact[1:]


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
    "get_supported_writer_formats",
    "get_writer_option_specs",
    "reader_factory",
    "register_openbabel_fallback",
    "register_reader",
    "register_reader_factory",
    "register_writer",
    "register_writer_factory",
    "select_reader",
    "upgrade_coords_to_graph",
    "write",
    "write_frame",
    "writer_factory",
    "WriterOptionSpec",
]
