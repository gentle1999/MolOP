# CLI first steps

Summarize, filter, and convert a batch without writing a Python script.

## Inspect commands

```bash
molop --help
molop parse --help
```

??? example "Help summary"

    ```text
    Usage: molop [OPTIONS] COMMAND [ARGS]...
    Commands: completion, parse
    Usage: molop parse [OPTIONS] PATTERN COMMAND1 [ARGS]...
    Commands: copy-to, draw-grid-image, filter-by-codec, filter-state,
              filter-value, format-transform, groupby, move-to, sample,
              to-summary-df
    ```

Every data operation starts with `molop parse PATTERN`. `PATTERN` can be a path or a shell glob; add
repeatable `--input PATH_OR_GLOB` options when one command must cover several input groups.

## Export a summary CSV

```bash
molop -q parse "results/*.out" \
  to-summary-df --full --out summary.csv
```

The default selects the final frame of each file. Add `--frame all` to export every frame.

For the bundled example:

```bash
molop -q parse "water_mp2.out" --n-jobs 1 \
  to-summary-df --full --out summary.csv
```

??? example "Output"

    ```text
    Summary written to summary.csv
    ```

## Filter before exporting

```bash
molop -q parse "results/*.log" \
  filter-state --state normal \
  filter-state --state opt \
  to-summary-df --full --out optimized.csv
```

Operations run from left to right. `filter-state` keeps a batch state, so another filter or one
terminal output operation can follow it.

??? example "Created file"

    ```text
    optimized.csv
    ```

## Export final structures

```bash
mkdir -p structures
molop -q parse "results/*.out" \
  filter-state --state normal \
  format-transform --format xyz --output-dir structures
```

Providing `--output-dir` creates the directory and writes files. Use `--no-write` to preview output
in the terminal without writing files.

With `water_mp2.out`, the output directory contains:

??? example "Created file"

    ```text
    structures/water_mp2.xyz
    ```

## Resolve an ambiguous extension

```bash
molop parse "results/*.out" --parser-detection orcaout \
  to-summary-df --out summary.csv
```

??? example "Created file"

    ```text
    summary.csv
    ```

## Next steps

- [CLI task recipes](../guides/cli-recipes.md)
- [Filter and select](../guides/filtering.md)
- [Full CLI reference](../command_line_interface.md)
