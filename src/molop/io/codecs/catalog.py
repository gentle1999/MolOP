"""Codec discovery and activation.

This module is imported only during lazy activation
(`Registry.ensure_default_codecs_registered()`), not at package import time.

Design goals:
- Deterministic behavior across platforms and installs.
- Minimal "linkage" for adding builtin codecs: co-locate registration next to
  the parser/model implementation (no extra per-format codec glue modules).
- Preserve lazy import boundary: importing `molop.io.codecs` must not import
  codec implementations nor enumerate entry points.

Policies:
- Builtin scan scope: public `*File.py` and `*FileParser.py` modules under
  `molop.io.logic`, discovered recursively.
- Exclusions: module basenames starting with '_' are skipped.
- Ordering: deterministic sort by fully-qualified module name.
- Builtin failures: fail hard (raise) if import/register fails.
- Plugin failures: warn + skip per entry point (do not fail activation).
"""

from __future__ import annotations

import warnings
from importlib import import_module
from typing import Any


_BUILTIN_SCAN_ROOT_PACKAGE_NAME = "molop.io.logic"

_ENTRY_POINT_GROUP = "molop.codecs"


class CodecPluginWarning(RuntimeWarning):
    pass


_loaded_builtin_registries: set[int] = set()
_loaded_plugin_registries: set[int] = set()


def _discover_builtin_codec_module_names() -> list[str]:
    """Return candidate builtin codec module names (unsorted).

    Discovery is convention-based: concrete file models and file parsers live in
    public modules ending with `File` or `FileParser`.
    """

    import importlib
    import pkgutil

    discovered: list[str] = []
    root_pkg = importlib.import_module(_BUILTIN_SCAN_ROOT_PACKAGE_NAME)
    root_path = getattr(root_pkg, "__path__", None)
    if root_path is None:
        return discovered

    prefix = f"{root_pkg.__name__}."
    for mod in pkgutil.walk_packages(root_path, prefix):
        if mod.ispkg:
            continue
        parts = mod.name.split(".")
        base = parts[-1]
        if any(part.startswith("_") for part in parts):
            continue
        if not base.endswith(("File", "FileParser")):
            continue
        discovered.append(mod.name)
    return discovered


def _register_module_codecs(module_name: str, registry: Any) -> None:
    module = import_module(module_name)
    register = getattr(module, "register", None)
    if register is None:
        # Optional: allow non-registrar modules within scanned packages.
        return
    if not callable(register):
        raise AttributeError(
            f"Builtin module {module_name!r} has non-callable attribute 'register'."
        )
    register(registry)


def _register_openbabel_fallback(registry: Any) -> None:
    from molop.io.codecs.openbabel_reader import register_openbabel_fallback_reader

    register_openbabel_fallback_reader(registry)


def _register_cml_writer(registry: Any) -> None:
    from molop.io.codecs.cml_codec import register as register_cml

    register_cml(registry)


def _iter_entry_points(group: str):
    """Return EntryPoint objects for a group (handles API differences)."""

    from importlib.metadata import entry_points

    eps = entry_points()
    select = getattr(eps, "select", None)
    if callable(select):
        return eps.select(group=group)
    # Older API: dict-like mapping group -> list[EntryPoint]
    try:
        return eps.get(group, ())
    except Exception:
        return ()


def _load_plugins_from_entry_points(registry: Any) -> None:
    entry_points = list(_iter_entry_points(_ENTRY_POINT_GROUP))
    entry_points.sort(key=lambda ep: (getattr(ep, "name", ""), getattr(ep, "value", "")))

    for ep in entry_points:
        name = getattr(ep, "name", "<unknown>")
        value = getattr(ep, "value", "<unknown>")
        try:
            loaded = ep.load()
        except Exception as exc:
            warnings.warn(
                f"Codec plugin entry point {name!r} ({value}) failed to load: {exc}",
                category=CodecPluginWarning,
                stacklevel=2,
            )
            continue

        register = loaded if callable(loaded) else getattr(loaded, "register", None)
        if not callable(register):
            warnings.warn(
                f"Codec plugin entry point {name!r} ({value}) has no callable register(registry).",
                category=CodecPluginWarning,
                stacklevel=2,
            )
            continue
        try:
            register(registry)
        except Exception as exc:
            warnings.warn(
                f"Codec plugin entry point {name!r} ({value}) register() failed: {exc}",
                category=CodecPluginWarning,
                stacklevel=2,
            )


def load_plugin_codecs(registry) -> None:
    """Load and register third-party codecs discovered via entry points."""

    reg_id = id(registry)
    if reg_id in _loaded_plugin_registries:
        return
    _load_plugins_from_entry_points(registry)
    _loaded_plugin_registries.add(reg_id)


def load_builtin_codecs(registry) -> None:
    """Load and register builtin codecs into the provided registry.

    This function is idempotent per-registry instance.
    """

    reg_id = id(registry)
    if reg_id in _loaded_builtin_registries:
        return

    module_names = sorted(set(_discover_builtin_codec_module_names()))
    for module_name in module_names:
        _register_module_codecs(module_name, registry)
    _register_cml_writer(registry)
    _register_openbabel_fallback(registry)

    _loaded_builtin_registries.add(reg_id)


__all__ = [
    "CodecPluginWarning",
    "load_builtin_codecs",
    "load_plugin_codecs",
]
