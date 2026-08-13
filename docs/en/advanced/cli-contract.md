# CLI Chain Contract

This page explains the `molop parse` chain rules and terminal-operation boundary. Start with [CLI
first steps](../getting_started/cli.md) for normal use; see the [CLI command reference](../command_line_interface.md)
for the complete option surface.

## Command shape

```text
molop parse PATTERN [parse options] OPERATION [operation options] ...
                                                     -> terminal operation
```

`PATTERN` is a path or glob; repeat `--input PATH_OR_GLOB` for additional input groups. After parsing,
operations that return a batch can be followed by more operations. Operations that return a summary,
rendered text, or path mapping must be last. The CLI checks this rule before reading files. The
`--report` parse mode is also terminal and cannot be combined with operations.

The examples below use the bundled [water_mp2.out](../../assets/examples/water_mp2.out).

## Common chains

=== "Filter by codec"

    ```bash
    molop -q parse "water_mp2.out" \
      --parser-detection orcaout \
      --n-jobs 1 \
      --output-format json \
      filter-by-codec --codec-id orcaout
    ```

    ??? example "Output shape"

        ```json
        [
          "<absolute path to water_mp2.out>"
        ]
        ```

=== "Summary"

    ```bash
    molop -q parse "water_mp2.out" --n-jobs 1 \
      filter-state --state normal \
      to-summary-df --full --out summary.csv
    ```

    ??? example "Terminal output and created file"

        ```text
        Summary written to summary.csv
        summary.csv
        ```

=== "Conversion"

    ```bash
    molop -q parse "water_mp2.out" --n-jobs 1 \
      format-transform --format xyz --output-dir converted
    ```

    ??? example "Created file"

        ```text
        converted/water_mp2.xyz
        ```

## Operation classes

| Class | Operations | Returned state |
| --- | --- | --- |
| Chainable | `filter-state`, `filter-value`, `filter-by-codec`, `sample` | New batch |
| Terminal | `to-summary-df` | CSV, JSON, or table result |
| Terminal | `format-transform` | Rendered text or generated files |
| Terminal | `draw-grid-image` | SVG/PNG image |
| Terminal | `groupby`, `copy-to`, `move-to` | Grouped paths or file-operation result |

A terminal operation cannot be followed by another operation. The following chain is rejected before
parsing:

```bash
molop parse "water_mp2.out" \
  to-summary-df --out summary.csv \
  filter-state --state normal
```

??? example "Error shape"

    ```text
    Error: to-summary-df returns non-FileBatchModelDisk and must be the last operation.
    ```

The exact error depends on the validation branch; move the filter before `to-summary-df` to correct it.

## Writer-specific options

`format-transform` accepts writer-specific options after its common options. The registered writer
provides the option names:

```bash
molop parse "water_mp2.out" --n-jobs 1 \
  format-transform --format gjf --output-dir gaussian_inputs \
  --route-section "#p B3LYP/6-31G(d) opt" \
  --link0-commands "%nprocshared=8"
```

??? example "Created file"

    ```text
    gaussian_inputs/water_mp2.gjf
    ```

Inspect the static and writer-specific options in the installed version:

```bash
molop parse PATTERN format-transform --help
```

??? example "Help output shape"

    ```text
    Usage: molop parse PATTERN format-transform [OPTIONS] [EXTRA_ARGS]...
    --format TEXT           Target writer format id.  [required]
    --output-dir DIRECTORY  Directory for generated files.
    --frame TEXT            Frame selection: all, int, or csv ints.
    ```

## Shell completion

Inspect the completion script without changing shell configuration:

```bash
molop completion show --shell bash
```

??? example "Output shape"

    ```text
    # bash completion script
    <script generated from the current Click command tree>
    ```

Install completion with `molop completion install --shell bash`, then reload the shell. Installation
modifies the current user's shell configuration; automation should use `completion show` and manage
the script explicitly.

## Related pages

- [CLI first steps](../getting_started/cli.md)
- [CLI command reference](../command_line_interface.md)
- [API Contracts](../reference/api_contracts.md)
