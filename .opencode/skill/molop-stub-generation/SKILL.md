---
name: MolOP Stub Generation
description: Use this skill when modifying MolOP stub generators, generated .pyi files, format_transform overloads, CLI transform stubs, or _render docstrings that should flow into generated typing/help output.
---

# MolOP Stub Generation

Use this skill when changing:

- `scripts/generate_io_typing_catalog.py`
- `scripts/generate_chemfile_format_transform_stubs.py`
- `scripts/generate_cli_transform_stubs.py`
- generated `.pyi` files under `src/molop/io/` or `src/molop/cli/`
- `_render(...)` docstrings that should appear in generated stubs

## Source-of-Truth Rules

- Never hand-edit generated `.pyi` files as the final fix.
- Fix the generator or the generator input instead.
- For transform stubs, `_render(...)` signatures and docstrings are the authoritative source.
- If a generated file differs from `--check`, regenerate it using the same `uv run` path the check uses.

## Canonical Commands

### IO typing catalog

```bash
uv run python scripts/generate_io_typing_catalog.py
uv run python scripts/generate_io_typing_catalog.py --check
```

### Chemfile format transform stubs

```bash
uv run python scripts/generate_chemfile_format_transform_stubs.py
uv run python scripts/generate_chemfile_format_transform_stubs.py --check
```

### CLI transform stubs

```bash
uv run python scripts/generate_cli_transform_stubs.py
uv run python scripts/generate_cli_transform_stubs.py --check
```

### Wrapper targets

```bash
make -f makefile check-typing-stubs
make -f makefile check-format-transform-stubs
make -f makefile check-cli-transform-stubs
```

## Documentation Rules for Generated Stubs

- Generated stub docs should use NumPy-style section headers.
- Prefer:

```text
Parameters
----------
Returns
-------
Source
------
```

- Preserve full `_render(...)` docstring content when the goal is user understanding.
- Keep raw docstring content in the generated output without inventing misleading summaries.
- If generator-added wrapper sections exist, place raw render doc text where it best preserves readability and VS Code rendering.

## File-Level vs Frame-Level Stub Rules

- `src/molop/io/base_models/_format_transform.pyi` should describe both:
  - `FrameFormatTransformMixin`
  - `FormatTransformMixin`
- Frame and file overload sets must be generated separately.
- Batch transform stubs belong in:

```bash
src/molop/io/_batch_format_transform.pyi
```

- CLI transform overloads belong in:

```bash
src/molop/cli/app.pyi
```

## Import and Annotation Rules

- If a generated stub mentions types such as `GJFLink0Commands`, ensure the generator also emits the necessary imports.
- Annotation-dependent imports should be derived from the real source signatures, not hard-coded if avoidable.
- Keep generated imports deterministic and sorted.

## Determinism Rules

- Generator output must be stable across repeated runs.
- `--check` must pass immediately after generation.
- If formatting differs between `python` and `uv run python`, trust the `uv run` output and align generation to that path.
- Sorted overload order matters; do not emit format literals in a nondeterministic order.

## Stub Test Expectations

The invariant tests should continue to verify:

- literal overloads are sorted
- generated docstrings contain the expected section headings
- no duplicate parameter entries appear in generated wrapper docs
- generator `--check` mode is deterministic

When generator output format changes, update the invariant tests in the same change set.

## Recommended Workflow

When you modify a render docstring or transform signature:

1. Update the real source (`_render(...)`, writer registration, or generator logic)
2. Regenerate the relevant stubs
3. Run:

```bash
uv run python scripts/generate_chemfile_format_transform_stubs.py --check
uv run python scripts/generate_cli_transform_stubs.py --check
python3 -m pytest -o addopts='' tests/test_stub_generators_invariants.py
```

4. Read the generated `.pyi` sections directly to confirm the user-facing output is actually good

## Common Mistakes to Avoid

- editing generated `.pyi` files without fixing the generator
- keeping only a summary when users need the full `_render(...)` docstring
- using Google-style headings when the project has standardized on NumPy-style docs for stub rendering
- forgetting to update invariant tests after changing stub layout
- regenerating via the wrong command path and then wondering why `--check` fails
