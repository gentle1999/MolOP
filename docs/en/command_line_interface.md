# Command Line Interface

MolOP exposes a single business command: `molop parse`.

The command first parses files into a `FileBatchModelDisk` batch state, then
executes a checked chain of operations. Operations that return another
`FileBatchModelDisk` can be followed by more operations. Operations that return
another result must be the last operation, and this is checked before files are
parsed.

Examples use `uv run molop ...` for in-repo reproducibility. If MolOP is
installed in your environment, use `molop ...` directly.

## Global Options

- `-v`, `--verbose`: Enable verbose output.
- `-q`, `--quiet`: Enable quiet mode.
- `--version`: Show the version and exit.

## Discover Commands

```bash
uv run molop --help
uv run molop parse --help
```

## Shell Completion

Install shell completion for the current shell:

```bash
molop completion install
```

Use `--shell bash`, `--shell zsh`, or `--shell fish` to install for a specific
shell. Restart the shell after installation, or source the updated rc file.

To inspect the generated completion script without installing it:

```bash
molop completion show --shell bash
```

## Parse And Chain Operations

List parsed files after filtering by detected codec:

```bash
uv run molop -q parse "tests/test_files/orca/h2_grad_orca.inp" \
  --parser-detection orcainp \
  --n-jobs 1 \
  --output-format json \
  filter-by-codec --codec-id orcainp
```

Convert parsed files to another format:

```bash
uv run molop -q parse "tests/test_files/orca/h2_grad_orca.inp" \
  --parser-detection orcainp \
  --n-jobs 1 \
  format-transform --format xyz --output-dir .tmp/molop_xyz
```

When `format-transform` writes files through `--output-dir`, it does not print
the rendered file contents to stdout. The generated files are the operation
result.

Generate a summary table:

```bash
uv run molop -q parse "tests/test_files/orca/h2_grad_orca.inp" \
  --parser-detection orcainp \
  --n-jobs 1 \
  to-summary-df --out .tmp/molop_summary.csv
```

## Operation Reference

Chainable operations:

| Operation | Meaning |
| --- | --- |
| `filter-state` | Filter by calculation state. |
| `filter-value` | Filter by charge, multiplicity, or file format. |
| `filter-by-codec` | Filter by detected reader codec id. |
| `sample` | Randomly sample files from the current batch. |

Terminal operations:

| Operation | Meaning |
| --- | --- |
| `format-transform` | Transform files to another format. With `--output-dir`, files are written and stdout stays quiet. |
| `to-summary-df` | Build a summary table. |
| `draw-grid-image` | Render a molecule grid image. |
| `groupby` | Group files and print grouped paths. |
| `copy-to` | Copy current batch files to a directory. |
| `move-to` | Move current batch files to a directory. |

A terminal operation must be the last operation in the chain.

## Format-Specific Options

`format-transform` accepts additional writer-specific options after its normal
options:

```bash
uv run molop parse "input.log" \
  format-transform --format gjf --output-dir out \
  --route-section "#p B3LYP/6-31G(d) opt" \
  --link0-commands "%nprocshared=8"
```

These dynamic options come from the registered writer metadata. Shell
completion can suggest available option names for the selected `--format`, and
can also show short descriptions for supported options.

Use the command help for the static surface:

```bash
uv run molop parse PATTERN format-transform --help
```
