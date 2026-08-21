# CLI command reference

MolOP's business entry point is:

```text
molop [global options] parse PATTERN [parse options] operation [operation options] ...
```

Start with [CLI first steps](getting_started/cli.md) or
[CLI task recipes](guides/cli-recipes.md). Use this page for options and chaining rules.

The commands below use the [shared ORCA water example](../assets/examples/water_mp2.out). Download
it into the current directory before running them; replace `water_mp2.out` with your own input when
needed.

## Global options

| Option | Purpose |
| --- | --- |
| `-v, --verbose` | Show more detailed logs |
| `-q, --quiet` | Suppress progress and non-result output for scripts |
| `--max-jobs N` | Set a process-wide worker ceiling for parsing and operations |
| `--version` | Show the version |
| `-h, --help` | Show help |

## `parse` options

| Option | Default | Purpose |
| --- | --- | --- |
| `PATTERN` | required | One path or glob |
| `--input PATH_OR_GLOB` | none | Add another path or glob; repeat the option for more inputs |
| `--parser-detection` | `auto` | Automatic detection or an ID such as `g16log`, `orcaout`, `xtbout` |
| `-j, --n-jobs` | `-1` | Parser process count |
| `--output-format` | `text` | Final terminal output as `text` or `json` |
| `--total-charge` | source value | Override molecular charge for every input |
| `--total-multiplicity` | source value | Override spin multiplicity for every input |
| `--only-extract-structure` | off | Skip energies, frequencies, and other non-structural results |
| `--only-last-frame` | off | Retain only the final frame from each file |
| `--capture-source-evidence` | off | Capture supported source spans, hashes, and parser provenance |
| `--source-encoding` | `utf-8` | Strict text decoding and source-offset encoding |
| `--release-file-content / --keep-file-content` | release | Release or retain raw source text after parsing |
| `--force-unit-transform / --no-force-unit-transform` | configured value | Override unit conversion for this parse |
| `--graph-reconstruction-backend` | configured value | Use the `cpp` or `python` graph backend for this parse |
| `--reconstruction-failure-policy` | configured value | Use `raise` or retain a `return_suspicious` fallback graph |
| `--make-dative-bonds / --no-make-dative-bonds` | configured value | Override dative-bond reconstruction for this parse |
| `--report` | off | Return one parse outcome per supplied input instead of a batch |

`--n-jobs -1` is the default automatic setting. MolOP uses the process-aware CPU limit and the
optional `molopconfig.max_jobs` ceiling; pass `--n-jobs 1` when diagnosing a parser or native-library
issue. Every parse and operation command defaults independently to `-1`. Use global `--max-jobs N` to
cap the entire command, or pass a positive `--n-jobs` on one operation when it needs a lower limit.

Parsing keeps the full process-aware limit. Operations that may enter MolGR, including custom graph
callbacks, automatically use the separate two-thirds CPU budget.

The charge and multiplicity options apply to every matched input, so split heterogeneous inputs into
separate commands. `--only-extract-structure` is a fast path and deliberately omits non-structural
scientific results. The unit and topology options are immutable per-command overrides; when omitted,
they inherit `molopconfig` without changing process-global configuration.

Add multiple input specifications without relying on one broad glob:

```bash
molop -q parse "gaussian/*.log" \
  --input "orca/*.out" \
  --input selected.xyz \
  --only-last-frame \
  to-summary-df --full --out summary.csv
```

??? example "Created file"

    ```text
    Summary written to summary.csv
    ```

## Parse outcome report

Use `--report` to audit successful, missing, unsupported, empty, and failed inputs. It is a terminal
parse mode and cannot be combined with operation commands.

```bash
molop -q parse "results/*.log" \
  --input required.out \
  --report --output-format json > parse-report.json
```

??? example "Created file"

    ```text
    parse-report.json
    ```

The JSON object contains `summary` counts and an ordered `outcomes` array. Each outcome includes the
input index, absolute path, status, detected format, parser warnings, and a serializable failure when
applicable. The default text output is a tab-separated diagnostic table.

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
molop parse "water_mp2.out" to-summary-df --help
```

??? example "Help summary"

    ```text
    Usage: molop parse PATTERN to-summary-df [OPTIONS]
    --mode [file|frame]       Summary mode.  [default: frame]
    --frame TEXT              Frame selection: all, int, or csv ints.
    --brief / --full          Use compact or full summary fields.
    --out FILE                Output file.
    --format [csv|json]       Output file format.  [default: csv]
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
molop -q parse "water_mp2.out" --parser-detection orcaout --n-jobs 1 \
  to-summary-df --full --format json --out summary.json
```

The output object includes the following fields (the full JSON also contains other parsed fields):

??? example "JSON output"

    ```json
    {
      "Calc Parameter.Software": "ORCA",
      "Status.IsNormal": true,
      "Energy.total_energy.hartree": -74.9993745981
    }
    ```

    Created file: `summary.json`.

## Transform options

```bash
molop parse "input.out" format-transform --help
```

??? example "Help summary"

    ```text
    Usage: molop parse PATTERN format-transform [OPTIONS] [EXTRA_ARGS]...
    --format TEXT           Target writer format id.  [required]
    --output-dir DIRECTORY  Directory for generated files.
    --frame TEXT            Frame selection: all, int, or csv ints.
    --write / --no-write    Write generated files.
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

??? example "Created file"

    ```text
    next/input.inp
    ```

## Inspect exact help

```bash
molop --help
molop parse --help
molop parse PATTERN filter-state --help
molop parse PATTERN to-summary-df --help
molop parse PATTERN format-transform --help
```

??? example "Common help shape"

    ```text
    molop [OPTIONS] parse PATTERN OPERATION [OPTIONS]
    ```

See the [Advanced CLI contract](advanced/cli-contract.md) for internal plan validation, terminal
operation constraints, and dynamic completion behavior.
