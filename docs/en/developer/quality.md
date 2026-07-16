# Development environment and quality gates

This page is for source contributors, not end-user installation.

## Create the environment

```bash
git clone https://github.com/gentle1999/MolOP.git
cd MolOP
uv sync
```

## Common checks

```bash
uv run ruff format . --check
uv run ruff check .
uv run mypy src
uv run pyright
uv run pytest --no-cov -q
```

Run the repository aggregate gate with:

```bash
make check
```

## Documentation checks

```bash
uv run rumdl check docs README.md README.zh.md
uv run pytest -q tests/test_documentation_examples.py --no-cov
uv run pytest -q tests/format_feature_coverage --no-cov
NO_MKDOCS_2_WARNING=1 uv run mkdocs build --strict
```

Hooks and tests maintain format capability pages. Reader/writer changes require format coverage tests,
not only a site build.

## Before committing

```bash
git diff --check
git status --short
```

Confirm that local logs, build directories, and unrelated user changes are not staged.
