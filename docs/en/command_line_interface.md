# CLI command reference

MolOP's business entry point is:

```text
molop [global options] parse PATTERN [parse options] operation [operation options] ...
```

Start with [CLI first steps](getting_started/cli.md) or
[CLI task recipes](guides/cli-recipes.md). Use this page for options and chaining rules.

## Global options

| Option | Purpose |
| --- | --- |
| `-v, --verbose` | Show more detailed logs |
| `-q, --quiet` | Suppress progress and non-result output for scripts |
| `--version` | Show the version |
| `-h, --help` | Show help |

## `parse` options

| Option | Default | Purpose |
| --- | --- | --- |
| `PATTERN` | required | One path or glob |
| `--parser-detection` | `auto` | Automatic detection or an ID such as `g16log`, `orcaout`, `xtbout` |
| `-j, --n-jobs` | `-1` | Parser process count |
| `--output-format` | `text` | Final terminal output as `text` or `json` |

## Chainable operations

| Operation | Required options | Result state |
| --- | --- | --- |
| `filter-state` | `--state ts|error|opt|normal|thermal|no-img` | New batch |
| `filter-value` | `--target charge|multiplicity|format --value VALUE` | New batch |
| `filter-by-codec` | `--codec-id FORMAT_ID` | New batch |
| `sample` | `--n N [--seed SEED]` | New batch |

Another filter or one terminal operation can follow these operations.

## Terminal operations

| Operation | Main options | Output |
| --- | --- | --- |
| `to-summary-df` | `--frame`, `--full`, `--out`, `--format` | CSV/JSON or terminal table |
| `format-transform` | `--format`, `--output-dir`, `--frame`, `--write/--no-write` | Rendered text or files |
| `draw-grid-image` | `--out` | PNG/SVG image |
| `groupby` | grouping options | Grouped paths |
| `copy-to` | destination | Copied source files |
| `move-to` | destination | Moved source files |

A terminal operation must end the operation chain.

## Summary options

```bash
molop parse "results/*.log" to-summary-df --help
```

Common options:

- `--mode frame|file`: frame- or file-level summary.
- `--frame -1|all|0,2`: final, all, or selected frames.
- `--brief/--full`: base columns or extended scientific results.
- `--flatten-columns/--multi-index-columns`: flat or three-level columns.
- `--on-missing-frame skip|error`: missing frame behavior.
- `--out PATH --format csv|json`: output path and format.

Run:

```bash
molop -q parse "water_mp2.out" \
  to-summary-df --full --format json
```

The output object includes:

```json
{
  "Calc Parameter.Software": "ORCA",
  "Status.IsNormal": true,
  "Energy.total_energy.hartree": -74.9993745981
}
```

## Transform options

```bash
molop parse "input.out" format-transform --help
```

- `--format FORMAT_ID` is required.
- `--frame -1|all|0,2` selects frames.
- `--embed/--no-embed` controls multi-frame embedding.
- `--output-dir DIR` selects a directory and enables writing by default.
- `--write` writes beside the source when no output directory is given.
- `--no-write` returns rendered output even when an output directory is present.

Writer-specific options follow static options:

```bash
molop parse "input.out" \
  format-transform --format orcainp --output-dir next \
  --keywords "B3LYP def2-SVP Opt" --nprocs 8 --maxcore 2000
```

## Inspect exact help

```bash
molop --help
molop parse --help
molop parse PATTERN filter-state --help
molop parse PATTERN to-summary-df --help
molop parse PATTERN format-transform --help
```

See the [Advanced CLI contract](advanced/cli-contract.md) for internal plan validation, terminal
operation constraints, and dynamic completion behavior.
