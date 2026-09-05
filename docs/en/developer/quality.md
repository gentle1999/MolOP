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
uv run rumdl check README.md README.zh.md $(rg --files docs -g '*.md')
uv run pytest -q tests/test_documentation_examples.py --no-cov
uv run pytest -q tests/format_feature_coverage --no-cov
NO_MKDOCS_2_WARNING=1 uv run mkdocs build --strict
```

Hooks and tests maintain format capability pages. Reader/writer changes require format coverage tests,
not only a site build.

## Release workflow

Releases are created manually on GitHub. After the release changes are merged into `main`, open the
[new release page](https://github.com/gentle1999/MolOP/releases/new), draft a release, create a `v*` tag
such as `v0.2.13` from `main`, and publish it.

The resulting `v*` tag event triggers `.github/workflows/ci.yaml` and
`.github/workflows/docs-deploy.yml`. After all quality, build, documentation, and test gates pass, CI
uploads the distributions to the existing GitHub Release and publishes them to PyPI through trusted
publishing. The documentation workflow publishes the versioned documentation and updates the `latest`
alias. Verify both workflow runs before announcing the version.

## Before committing

```bash
git diff --check
git status --short
```

Confirm that local logs, build directories, and unrelated user changes are not staged.
