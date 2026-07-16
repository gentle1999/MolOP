# CLI task recipes

These commands target an installed MolOP environment and all start with `molop parse PATTERN`.

## Summarize final frames

```bash
molop -q parse "results/*.log" \
  to-summary-df --full --out summary.csv
```

Result: `summary.csv` is created in the current directory with one row per successfully parsed file
and columns such as `Status.IsNormal` and `Energy.total_energy.hartree`.

## Summarize every frame

```bash
molop -q parse "results/*.log" \
  to-summary-df --frame all --full --out trajectory.csv
```

Result: `trajectory.csv` contains one row per frame and `General.FrameID` identifies each frame.

## Keep normal optimized results

```bash
molop -q parse "results/*.log" \
  filter-state --state normal \
  filter-state --state opt \
  to-summary-df --full --out optimized.csv
```

Result: `optimized.csv` contains only files that pass both filters.

## Select ORCA outputs

```bash
molop -q parse "results/*" \
  filter-by-codec --codec-id orcaout \
  to-summary-df --out orca.csv
```

`filter-by-codec` uses the actual reader ID rather than the extension string.

## Export final XYZ files

```bash
molop -q parse "results/*.out" \
  filter-state --state normal \
  format-transform --format xyz --output-dir structures
```

Result: `.xyz` files named from their sources are written under `structures/`. The CLI creates a
missing output directory.

## Preview without writing

```bash
molop -q parse "water_mp2.out" \
  format-transform --format xyz --no-write
```

The terminal displays the source path and rendered XYZ text. No `.xyz` file is created.

## Produce JSON

```bash
molop -q parse "water_mp2.out" --output-format json \
  filter-state --state normal \
  sample --n 1 --seed 1
```

Output is a path array:

```json
[
  "/absolute/path/to/water_mp2.out"
]
```

## Rules

- Operations run in written order.
- `filter-state`, `filter-value`, `filter-by-codec`, and `sample` can be chained further.
- `to-summary-df`, `format-transform`, `copy-to`, `move-to`, `groupby`, and `draw-grid-image` are
  terminal operations; no business operation may follow them.
- Use `molop parse --help` for the operation list and a specific operation's `--help` for options.

See the [CLI command reference](../command_line_interface.md) for exact constraints.
