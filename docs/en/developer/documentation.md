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

## Output rendering

Treat output as a separate, collapsible part of every Markdown example. In the
documentation site, use a Material details block:

````markdown
??? example "Output"
    ```text
    ...
    ```
````

Keep commands, input code, and explanations outside the block. This applies to
printed values, returned text, command output, generated file trees, and other
artifacts that are useful for verification. README pages use native HTML
`<details>` because GitHub does not render Material admonitions.

Notebook output has a different ownership boundary: edit code and Markdown
cells, then let CI execute the notebooks with `nbconvert --execute` before the
MkDocs build. `mkdocs-jupyter` renders the saved notebook as HTML with
`execute: false`; do not hand-edit a notebook's `outputs`, execution counts,
or generated HTML.

To show one notebook cell's HTML output directly in a Markdown page, reference its cell `id`:

```markdown
<!-- notebook-output: examples/02-batch-summary-filter-select.ipynb#batch-summary -->
```

`scripts/mkdocs_hooks.py` reads the executed notebook's `text/html` output and embeds it in the
matching locale page. The cell must have a stable `id` and produce HTML after CI execution. Do not
copy DataFrame HTML into Markdown.

## New pages

One change should include:

- the `docs/zh/...` Chinese page;
- the semantically localized `docs/en/...` page;
- `mkdocs.yml` navigation;
- tests for executable code.

## Validate

```bash
uv run rumdl check README.md README.zh.md $(rg --files docs -g '*.md')
uv run pytest -q tests/test_documentation_examples.py --no-cov
NO_MKDOCS_2_WARNING=1 uv run mkdocs build --strict
```
