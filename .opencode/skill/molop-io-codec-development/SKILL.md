---
name: MolOP IO Codec Development
description: Use this skill when changing MolOP parsers, frame/file models, codec registration, format transformation, or render/write behavior. It captures the repository's current IO architecture and the constraints that keep frame/file conversion logic consistent.
---

# MolOP IO Codec Development

Use this skill for work in:

- `src/molop/io/codec_registry.py`
- `src/molop/io/codecs/`
- `src/molop/io/base_models/`
- `src/molop/io/logic/*_models/`
- `src/molop/io/logic/*_frame_models/`
- CLI or batch conversion paths that depend on format rendering

## Core Architecture Rules

- File-level rendering goes through `codec_registry.write(...)`.
- Frame-level rendering goes through `codec_registry.write_frame(...)`.
- Do not wrap a frame as a fake one-frame file just to reuse file writers.
- Writer registration must reflect the correct domain:
  - `domain="file"`
  - `domain="frame"`
- Keep format behavior owned by the format codec whenever practical.

## File vs Frame Responsibilities

### File-level

File objects are responsible for:

- holding multiple frames
- frame selection (`frameID`)
- `embed_in_one_file` behavior
- file-path derivation for multi-frame outputs

### Frame-level

Frame objects are responsible for:

- rendering a single structure/input/result
- exposing the format-specific `_render(...)` logic
- carrying the canonical docstring that explains format-specific parameters

If you need single-frame rendering semantics, implement them at the frame writer layer, not by faking file containers.

## Codec Registry Rules

- Add new writer registrations through `@registry.writer_factory(...)`.
- Always decide explicitly whether a writer is for `file` or `frame`.
- If a format supports both file and frame transforms, register both.
- Preserve `graph_policy` behavior when adding or changing writers.
- If the format requires graph-level rendering, make the registry path perform the correct upgrade or fail explicitly.

## Render Logic Placement

- Put format-specific behavior in the relevant codec or frame model, not in unrelated helpers.
- If a codec can own its rendering semantics directly, prefer that over helper indirection.
- Shared helper code should remain generic; avoid adding one-off format wrappers there unless the abstraction is genuinely reusable.

## Docstring Rules for Render Methods

- `_render(...)` docstrings are the source of truth for format-specific documentation.
- Use NumPy-style sections:

```text
Parameters
----------
```

- Explain not only parameter type but actual behavior and fallback precedence.
- If `_render(...)` can use instance defaults when parameters are `None`, say that explicitly.

## Sequence and Iterator Rules

MolOP has had multiple bugs from stateful container iteration.

For file/frame containers, batch models, orbital collections, and vibration collections:

- Containers should be re-iterable.
- `__iter__` must return a fresh iterator.
- Avoid `__next__` on plain containers.
- Do not use shared cursor state like `_index_` unless the object is a real iterator.
- If filtering/grouping depends on positional alignment, snapshot the source sequence first.

## Validation Rules

After IO/codec changes, prefer this order:

```bash
make -f makefile format-check
make -f makefile lint-check
make -f makefile type-check
make -f makefile pyright
```

If transform or registry behavior changed, also run:

```bash
make -f makefile check-format-transform-stubs
make -f makefile check-cli-transform-stubs
make -f makefile check-typing-stubs
```

Then run the most relevant targeted test or a focused runtime check.

## High-Value Regression Checks

When touching parser/render/iterator logic, explicitly verify at least one of:

- nested iteration on a container still works
- repeated traversal returns the same values twice
- frame-level `format_transform(...)` routes through `write_frame(...)`
- file-level `format_transform(...)` routes through `write(...)`
- generated text and file outputs are both correct
- generator-based parallel logic preserves alignment with the source snapshot

## Common Mistakes to Avoid

- registering only a file writer when frame transform should exist too
- putting format-specific behavior in a shared helper instead of the codec
- reintroducing fake frame wrappers for file-oriented code paths
- mixing container semantics with shared iteration state
- changing `_render(...)` parameters without regenerating the stubs
