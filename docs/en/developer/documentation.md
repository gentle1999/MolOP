# Documentation contributions

User documentation starts from tasks, and core Chinese/English pages are maintained in pairs.

## Page order

1. State the task and input.
2. Provide a 5-20 line shortest runnable example.
3. Follow it with real output, return shape, or generated file tree.
4. Explain common variants and actual limitations.
5. Link to the format matrix or exact contract.

End-user examples use `molop`; only development environment commands use `uv run molop`. The quick
start must not depend on `tests/test_files/`, a source root, or test environment variables.

## Sources of truth

- API names and defaults: current Python source and type hints.
- CLI: `molop --help` and operation-specific help.
- Format capability: `tests/format_feature_coverage/support_matrix.py`.
- Core examples: `docs/assets/examples/` and `tests/test_documentation_examples.py`.

Do not hand-maintain a second supported/unsupported matrix in prose.

## New pages

One change should include:

- the `docs/zh/...` Chinese page;
- the semantically localized `docs/en/...` page;
- `mkdocs.yml` navigation;
- tests for executable code.

## Validate

```bash
uv run rumdl check docs README.md README.zh.md
uv run pytest -q tests/test_documentation_examples.py --no-cov
NO_MKDOCS_2_WARNING=1 uv run mkdocs build --strict
```
