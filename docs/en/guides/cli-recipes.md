# CLI task recipes

These commands target an installed MolOP environment and all start with `molop parse PATTERN`.

## Summarize final frames

```bash
molop -q parse "results/*.log" \
  to-summary-df --full --out summary.csv
```

Result: `summary.csv` is created in the current directory with one row per successfully parsed file
and columns such as `Status.IsNormal` and `Energy.total_energy.hartree`.

With the bundled `water_mp2.out` example, the terminal reports:

??? example "Terminal output"

    ```text
    Summary written to summary.csv
    ```

    Created file:

    ```text
    summary.csv
    ```

## Summarize every frame

```bash
molop -q parse "results/*.log" \
  to-summary-df --frame all --full --out trajectory.csv
```

Result: `trajectory.csv` contains one row per frame and `General.FrameID` identifies each frame.

??? example "Created file"

    ```text
    trajectory.csv
    ```

## Keep normal optimized results

```bash
molop -q parse "results/*.log" \
  filter-state --state normal \
  filter-state --state opt \
  to-summary-df --full --out optimized.csv
```

Result: `optimized.csv` contains only files that pass both filters.

??? example "Created file"

    ```text
    optimized.csv
    ```

## Select ORCA outputs

```bash
molop -q parse "results/*" \
  filter-by-codec --codec-id orcaout \
  to-summary-df --out orca.csv
```

`filter-by-codec` uses the actual reader ID rather than the extension string.

??? example "Created file"

    ```text
    orca.csv
    ```

## Export final XYZ files

```bash
molop -q parse "results/*.out" \
  filter-state --state normal \
  format-transform --format xyz --output-dir structures
```

Result: `.xyz` files named from their sources are written under `structures/`. The CLI creates a
missing output directory.

??? example "Created files"

    ```text
    structures/
      water_mp2.xyz
    ```

## Preview without writing

```bash
molop -q parse "water_mp2.out" \
  format-transform --format xyz --no-write
```

The terminal displays the source path and rendered XYZ text. No `.xyz` file is created.

For the shared example, the rendered content is:

??? example "Rendered XYZ"

    ```text
    3
    comment charge 0 multiplicity 1
    O               1.7849140000      1.2624220000      0.5119850000
    H               2.6482370000      1.0729290000      0.1316310000
    H               1.1831680000      1.2568160000     -0.2388350000
    ```

## Produce JSON

```bash
molop -q parse "water_mp2.out" --output-format json \
  filter-state --state normal \
  sample --n 1 --seed 1
```

Output is a path array. The path is absolute and depends on the environment; the stable shape is
shown below rather than a machine-specific prefix:

??? example "JSON output shape"

    ```json
    [
      "<absolute-path-to-working-directory>/water_mp2.out"
    ]
    ```

## Rules

- Operations run in written order.
- `filter-state`, `filter-value`, `filter-by-codec`, and `sample` can be chained further.
- `to-summary-df`, `format-transform`, `copy-to`, `move-to`, `groupby`, and `draw-grid-image` are
  terminal operations; no business operation may follow them.
- Use `molop parse --help` for the operation list and a specific operation's `--help` for options.

See the [CLI command reference](../reference/cli.md) for exact constraints.
