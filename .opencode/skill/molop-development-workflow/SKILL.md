---
name: MolOP Development Workflow
description: Use this skill when developing inside the MolOP repository. It captures the project's actual uv-based workflow, lowercase makefile conventions, stub generator rules, and validation sequence for parser/model/codec changes.
---

# MolOP Development Workflow

This skill defines the expected development workflow for the MolOP repository.

Use it whenever you modify parsing, rendering, codecs, CLI transforms, data models, typing stubs, or project tooling in this repo.

## Core Rules

- Use `uv` for Python execution in this repo.
- The project task wrapper is `makefile`, not `Makefile`.
- Prefer `make -f makefile ...` for standard project workflows.
- Treat `make check` as non-mutating verification.
- If you change generated artifacts or their source logic, regenerate them immediately.
- Keep sequence-like containers stateless; do not mix container semantics with shared iterator state.

## Canonical Commands

### Environment

```bash
uv sync --all-groups
```

If a command should run inside the project environment, prefer:

```bash
uv run python ...
uv run pytest ...
uv run ruff ...
uv run mypy ...
uv run pyright
```

### Make Targets

Always invoke the project wrapper explicitly:

```bash
make -f makefile <target>
```

Most important targets:

```bash
make -f makefile format
make -f makefile format-check
make -f makefile lint
make -f makefile lint-check
make -f makefile type-check
make -f makefile pyright
make -f makefile check
make -f makefile test
make -f makefile check-typing-stubs
make -f makefile check-format-transform-stubs
make -f makefile check-cli-transform-stubs
```

### Quality Contract

- `format` and `lint` may modify files.
- `format-check`, `lint-check`, and `check` must not modify files.
- `check` currently means:
  - formatting is already correct
  - lint passes without fixes
  - mypy passes
  - pyright passes
  - generated stubs are up to date

## Generated Files and Their Sources

If you touch the generator inputs, regenerate the outputs in the same work session.

### IO typing catalog

Source:

```bash
scripts/generate_io_typing_catalog.py
```

Generate/check:

```bash
uv run python scripts/generate_io_typing_catalog.py
uv run python scripts/generate_io_typing_catalog.py --check
```

Output:

```bash
src/molop/io/_typing_catalog.pyi
```

### Chemfile format_transform stubs

Source:

```bash
scripts/generate_chemfile_format_transform_stubs.py
```

Generate/check:

```bash
uv run python scripts/generate_chemfile_format_transform_stubs.py
uv run python scripts/generate_chemfile_format_transform_stubs.py --check
```

Outputs:

```bash
src/molop/io/base_models/_format_transform.pyi
src/molop/io/_batch_format_transform.pyi
```

### CLI transform stubs

Source:

```bash
scripts/generate_cli_transform_stubs.py
```

Generate/check:

```bash
uv run python scripts/generate_cli_transform_stubs.py
uv run python scripts/generate_cli_transform_stubs.py --check
```

Output:

```bash
src/molop/cli/app.pyi
```

## Format Transform Rules

- File-level transforms go through `codec_registry.write(...)`.
- Frame-level transforms go through `codec_registry.write_frame(...)`.
- Do not fake frame rendering by wrapping a frame as a one-frame file.
- Stub docs for format transforms should preserve the original `_render` docstring content.
- Generated stub section headings should use NumPy-style headings such as:

```text
Parameters
----------
```

## Iterator and Container Rules

MolOP has already had real bugs from stateful iterator implementations.

When a class is conceptually a container or `Sequence`, follow these rules:

- `__iter__` should return a fresh iterator.
- Do not store shared iteration cursors like `_index_` unless the object is truly an iterator.
- Avoid implementing `__next__` on plain containers.
- Nested iteration and repeated traversal must produce stable results.
- If parallel work depends on positional alignment, snapshot inputs first.

This rule is especially important for:

- batch models
- chemfile/frame containers
- orbital/vibration collections
- mask/filter/groupby logic paired with parallel execution

## Testing and Verification Strategy

Use the smallest correct verification set for the files you changed, then escalate if needed.

### First line

```bash
make -f makefile format-check
make -f makefile lint-check
```

### Typing

```bash
make -f makefile type-check
make -f makefile pyright
```

### Stub consistency

```bash
make -f makefile check-typing-stubs
make -f makefile check-format-transform-stubs
make -f makefile check-cli-transform-stubs
```

### Tests

Prefer targeted tests first.

```bash
uv run pytest path/to/test_file.py
```

If pytest is blocked by unrelated environment/plugin issues, document the blocker and run the narrowest direct runtime check you can justify.

## Repo-Specific Pitfalls

- The file is `makefile`, not `Makefile`.
- `pyproject.toml` sets pytest `addopts` with coverage; if you intentionally bypass that for a narrow diagnostic run, do it explicitly and say why.
- Importing the full package may pull in optional external dependencies such as `rdkit_dof` or `molgr`; for isolated runtime checks, stub only the missing external boundary, not project code.
- If `uv run ... --check` disagrees with a file you generated via another path, trust the `uv` path and regenerate with `uv run`.

## Recommended End-of-Task Checklist

For most MolOP code changes:

```bash
make -f makefile format-check
make -f makefile lint-check
make -f makefile check-format-transform-stubs
make -f makefile check-cli-transform-stubs
make -f makefile check-typing-stubs
```

Then run the most relevant targeted tests or direct runtime checks for the modified subsystem.

If you changed public rendering, parsing, or iteration semantics, include one explicit behavior check that proves the bug is actually fixed.
