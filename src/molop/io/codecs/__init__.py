"""Codec package.

Importing this package MUST NOT implicitly populate the global codec registry.

Builtin registration is co-located with original implementations and loaded via
`molop.io.codecs.catalog.load_builtin_codecs()` as part of lazy activation:
- `molop.io.logic.coords_parsers` and `molop.io.logic.QM_parsers` register reader codecs.
- `molop.io.logic.coords_models` and `molop.io.logic.QM_models` register writer codecs.

Third-party codecs may be discovered via Python entry points (group:
`molop.codecs`) and loaded via `molop.io.codecs.catalog.load_plugin_codecs()`.
"""


def load_builtin_codecs(registry=None) -> None:
    from molop.io.codec_registry import default_registry
    from molop.io.codecs.catalog import load_builtin_codecs as _load

    _load(default_registry if registry is None else registry)


def load_plugin_codecs(registry=None) -> None:
    from molop.io.codec_registry import default_registry
    from molop.io.codecs.catalog import load_plugin_codecs as _load

    _load(default_registry if registry is None else registry)


__all__ = [
    "load_builtin_codecs",
    "load_plugin_codecs",
]
