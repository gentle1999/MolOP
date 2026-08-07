# Reader/writer plugin development

MolOP has two IO extension modes: codecs builtin to the repository and independently distributed third-party plugins. Before adding a format, confirm that an existing dedicated codec or the OpenBabel fallback cannot meet the requirement.

See the [Parser contract](parser-contract.md) for the full File/Frame hierarchy, source spans, parse diagnostics, and renderer constraints. This page focuses on connecting a codec to the registry.

## Choose a registration mode

| Scenario | Registration mode |
| --- | --- |
| Format builtin to the MolOP repository | Define `register(registry)` in a public `*FileParser.py` or `*File.py` module under `src/molop/io/logic/`; the catalog discovers it automatically. |
| Independently installed third-party package | Declare a `molop.codecs` entry point in the package's own `pyproject.toml`. |
| Private host-application integration | Obtain the target `Registry` during application startup and call the plugin's `register(registry)` explicitly. |

Builtin formats do not modify MolOP's `pyproject.toml` and do not declare entry points. Third-party packages are outside the builtin scan and must use an entry point or explicit host registration.

## Minimal third-party package

Keep parsing, models, and registration separate:

```text
molop-myfmt/
  pyproject.toml
  src/
    molop_myfmt/
      __init__.py
      plugin.py
      reader.py
      writer.py
```

Declare the dependency and entry point in the plugin package's `pyproject.toml`:

```toml
[project]
name = "molop-myfmt"
version = "0.1.0"
dependencies = ["molop"]

[project.entry-points."molop.codecs"]
myfmt = "molop_myfmt.plugin:register"
```

The unconstrained dependency in this skeleton only illustrates the minimum package structure. Before publishing, replace it with the MolOP version range actually covered by plugin tests. Development documentation may be ahead of PyPI. If the plugin depends on an interface marked as unreleased in the banner, pin the displayed commit while developing, then publish the plugin with an updated lower bound after that interface enters a release. See [Documentation versions](../versioning.md).

The entry-point target must be a callable `register(registry)`. It may also resolve to a module object with a callable method of that name. Loading occurs when the registry is first activated; module import must not parse files, load large models, or mutate the global registry.

```python
from __future__ import annotations

from typing import TYPE_CHECKING

from .reader import MyFmtReader

if TYPE_CHECKING:
    from molop.io.codec_registry import Registry


def register(registry: Registry) -> None:
    @registry.reader_factory(
        format_id="myfmt",
        extensions={".myfmt"},
        priority=100,
    )
    def make_reader() -> MyFmtReader:
        return MyFmtReader()
```

One entry point may register multiple readers and writers. Keep factories lightweight; construct expensive resources only when the codec is actually selected.

## Reader contract

A reader provides at least `format_id`, `extensions`, `priority`, and `read(path, **kwargs)`, and returns a `ParseResult`. The value returned by `read()` must be a file model that can join a MolOP disk batch, not a bare dictionary, frame, or external-tool object.

```python
from __future__ import annotations

from pathlib import Path
from typing import Any

from molop.io.codec_types import ParseOptions, ParseResult, StructureLevel

from .parser import MyFmtFileParserDisk


class MyFmtReader:
    format_id = "myfmt"
    extensions = frozenset({".myfmt"})
    priority = 100

    def probe_file_format(self, path: str | Path) -> bool:
        return MyFmtFileParserDisk.probe_file_format(path)

    def read(self, path: str | Path, **kwargs: Any) -> ParseResult[object]:
        supplied = kwargs.get("parse_options")
        if supplied is not None and not isinstance(supplied, ParseOptions):
            raise TypeError("parse_options must be a ParseOptions instance")
        options = (supplied or ParseOptions()).resolved()

        parser = MyFmtFileParserDisk(
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
        return ParseResult(
            value=value,
            level=StructureLevel.COORDS,
            detected_format=self.format_id,
        )
```

`ParseOptions` is the immutable configuration snapshot shared by the batch, reader, file parser, and frame parser. A plugin must prefer the supplied `parse_options` and must not reread global configuration for every frame. Resolve defaults locally only when the reader is called independently without that argument.

`StructureLevel.COORDS` means the minimum guarantee is atoms plus coordinates. Use `StructureLevel.GRAPH` only when the reader guarantees a valid molecular graph.

### Probes and exceptions

`probe_file_format(path)` is recommended but is not required by `ReaderCodec`. If it is absent, the registry treats the reader as a candidate and proceeds directly to full parsing. Readers for shared extensions and fallback readers should therefore provide a fast, side-effect-free probe.

- Raise `FormatMismatchError` when the content does not belong to this format so the batch can try the next reader.
- Raise a specific parse exception when the file belongs to the format but is corrupt or the parser fails; do not disguise that failure as a format mismatch.
- Put recoverable issues in `ParseResult.warnings` instead of only logging them.

```python
from molop.io.codec_types import ParseResult, ParseWarning, StructureLevel

return ParseResult(
    value=value,
    level=StructureLevel.COORDS,
    warnings=(
        ParseWarning(
            code="MYFMT.MISSING_OPTIONAL_BLOCK",
            message="Optional property block was not present.",
        ),
    ),
    detected_format="myfmt",
)
```

These warnings survive in the per-file outcomes returned by `AutoParser(..., return_report=True)`. The default `AutoParser()` API still returns only the successfully parsed batch.

### Parse state and source loading

A frame parser must read its current block, file/segment metadata, and `ParseOptions` from the immutable `FrameParseContext`. Do not store the current block, cursor, or metadata on the parser instance. Context metadata overrides same-named fields from the format payload during frame assembly.

Parsers based on `BaseFileParserDisk` normally call `parse(path)`. When a caller already owns the source content, use `parse_bytes(..., file_path=...)` or `parse_decoded_source(..., file_path=...)` to retain disk-source identity and exact source offsets without reopening the file. Do not read an entire file merely for a probe; probe/full-read reuse belongs to a caller that explicitly owns the preloaded source.

## Writer contract

A writer provides at least `format_id`, `priority`, `required_level`, and `write(value, **kwargs)`. Factory registration also declares the output domain and graph policy:

```python
def register(registry: Registry) -> None:
    @registry.writer_factory(
        format_id="myfmt",
        required_level=StructureLevel.COORDS,
        domain="file",
        default_graph_policy="coords",
        priority=100,
    )
    def make_writer() -> MyFmtWriter:
        return MyFmtWriter()
```

- Use `domain="file"` for file composition and multi-frame embedding or splitting.
- Use `domain="frame"` only when the format has a valid independent single-frame representation.
- Set `default_graph_policy` deliberately to match the format's requirements.
- Preserve important format-specific directives or source blocks where possible to support round trips.

Dynamic CLI writer options are inferred from named parameters on `write()` and the associated file/frame `_render()` methods. Give those parameters accurate types, defaults, and docstrings: `Literal` supplies value candidates, while `bool` supplies a boolean switch. Do not hide user-facing options inside `**kwargs`.

## Loading and precedence

Readers and writers are selected by descending `priority`; equal priorities retain deterministic registration order. Use a priority above a builtin codec only when overriding it is intentional, and add tests for extension conflicts.

Third-party entry points are sorted by name and value before loading. If one plugin fails to import, has no callable `register()`, or fails during registration, MolOP emits `CodecPluginWarning` and skips that plugin without blocking other codecs. Plugin tests should treat these warnings as installation or compatibility errors.

## Verification checklist

Cover at least these behaviors:

1. The installed entry point is discovered by the default registry and registers only once.
2. Both explicit `parser_detection="myfmt"` and extension-based automatic selection parse a fixture.
3. A probe mismatch falls through to the next reader, while a real parse error is not swallowed.
4. `ParseOptions` reaches the file/frame parser, and reader warnings appear in the structured report.
5. The returned disk file model joins a batch and supports summaries, frame selection, and `format_transform()`.
6. Writer tests cover string rendering, disk output, naming, file/frame domains, and failure behavior.
7. Run the contract suite against both the minimum supported MolOP release and current mainline.

Builtin codecs must also update paired format pages and user examples, then run the [quality gates](quality.md). Do not maintain a second format list in prose; capability status comes from `tests/format_feature_coverage/support_matrix.py` and generated sections on format pages.
