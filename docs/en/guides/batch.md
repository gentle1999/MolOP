# Batch summaries

Turn a batch of calculation outputs into a pandas DataFrame and export CSV or JSON.

## Shortest example

```python
from molop import AutoParser

batch = AutoParser("water_mp2.out", n_jobs=1)
summary = batch.to_summary_df(
    frame=-1,
    brief=False,
    flatten_columns=True,
)
summary.to_csv("summary.csv", index=False)
print(summary.shape)
```

??? example "Output"

    ```text
    (1, 23)
    ```

Replace `water_mp2.out` with `results/*.out` to process a batch. Each successfully parsed file
produces one row; inputs that do not parse produce no row.

## Parallel batch operations

Parsing and batch operations default to automatic parallelism with `n_jobs=-1`:

```python
from molop import AutoParser

batch = AutoParser("results/*.out")
normal = batch.filter_state("normal", n_jobs=-1)
summary = normal.to_summary_df(
    frame=-1,
    brief=False,
    flatten_columns=True,
    n_jobs=-1,
)
summary.to_csv("normal.csv", index=False)
```

`n_jobs=-1` follows the process-aware CPU limit. With the default `max_jobs=None`, it respects
observable scheduler/container quotas and CPU affinity where the platform exposes them; set
`molopconfig.max_jobs` before creating the batch when a lower ceiling is required. Use `n_jobs=1`
while diagnosing a parser or a native-library failure. CLI parse and operation commands have the same
independent default.

For a custom per-file operation, use `parallel_execute` and return explicit values:

```python
def file_format(parsed_file):
    return parsed_file.filename, parsed_file.detected_format_id

rows = batch.parallel_execute(file_format, n_jobs=-1, return_results=True)
print(list(rows))
```

For the bundled sample, the result is:

??? example "Output"

    ```text
    [('water_mp2.out', 'orcaout')]
    ```

## Check the shared example

```python
import pandas as pd
from IPython.display import display
from molop import AutoParser

batch = AutoParser("water_mp2.out", n_jobs=1)
summary = batch.to_summary_df(brief=False, flatten_columns=True)
with pd.option_context("display.max_columns", None, "display.max_rows", None, "display.width", 240):
    display(summary)
```

### Complete table output

<!-- notebook-output: examples/02-batch-summary-filter-select.ipynb#batch-summary -->

The table comes from the executed output of [Notebook 02](../examples/02-batch-summary-filter-select.ipynb)
and contains all 23 columns from `summary` without a temporary slice. CI executes the notebook
before MkDocs embeds the saved HTML.

## Final and all frames

```python
final_rows = batch.to_summary_df(frame=-1)
all_rows = batch.to_summary_df(frame="all")
selected_rows = batch.to_summary_df(frame=[0, -1])
print(final_rows.shape, all_rows.shape, selected_rows.shape)
```

??? example "Output"

    ```text
    (1, 19) (1, 19) (2, 19)
    ```

The default `on_missing_frame="skip"` omits unavailable indices. Use
`on_missing_frame="error"` when a missing frame should fail the operation.

## Brief and full fields

```python
brief = batch.to_summary_df(brief=True, flatten_columns=True)
full = batch.to_summary_df(brief=False, flatten_columns=True)
print(brief.columns[:7].tolist())
```

The brief summary includes stable file, molecule, calculation, and status columns. The first seven
columns for the shared sample are:

??? example "Columns"

    ```text
    ['DiskStorage.FilePath', 'DiskStorage.FileFormat', 'General.Charge',
     'General.Multiplicity', 'General.CanonicalSMILES', 'General.NumAtoms',
     'General.FrameID']
    ```

The full summary adds available `Energy.*`, `Thermal.*`, and `Vibration.*` columns.

## Flat and MultiIndex columns

```python
flat = batch.to_summary_df(flatten_columns=True)
print(flat["General.Charge"].tolist())

multi = batch.to_summary_df(flatten_columns=False)
print(multi[("General", "Charge", "")].tolist())
```

??? example "Output"

    ```text
    [0]
    [0]
    ```

Flat columns are easier for CSV, JSON, and general analysis. Keep MultiIndex columns when the
three-level group, field, and unit meaning matters.

## Group, sample, and preview

These operations keep the batch model and are useful before committing to a full export:

```python
groups = batch.groupby(lambda parsed_file: parsed_file.detected_format_id, n_jobs=-1)
print({key: len(group) for key, group in groups.items()})

sample = batch.sample(n=1, seed=1)
print(sample.file_names)

grid = batch.draw_grid_image(maxMols=16, useSVG=True, n_jobs=1)
print(type(grid).__name__, grid.lstrip()[:4])
```

For the bundled sample, the textual output is:

??? example "Output"

    ```text
    {'orcaout': 1}
    ['water_mp2.out']
    str <svg
    ```

With `useSVG=True`, `grid` is an SVG string that can be embedded in a notebook or written to an
`.svg` file. MolOP uses `rdkit-dof` by default for depth-aware drawing; set
`molopconfig.use_dof_effect_drawer = False` to use the standard RDKit drawer. Set `useSVG=False`
for a raster image. The CLI equivalent is `draw-grid-image --out structures.svg` or
`draw-grid-image --out structures.png`.

## Export JSON

```python
summary.to_json("summary.json", orient="records", indent=2)
print("summary.json")
```

??? example "Created file"

    ```text
    summary.json
    ```

## Large-batch guidance

- Start with 5-20 representative files and `n_jobs=1` to inspect fields.
- Increase `n_jobs` after confirming the format; disks and large files may limit scaling.
- Use `only_last_frame=True` only when optimization trajectories are unnecessary.
- Compare input paths with `batch.file_paths` to record omitted files.
- Group calculations by theory level before comparing `total_energy`.

## CLI equivalent

```bash
molop -q parse "results/*.out" \
  to-summary-df --full --flatten-columns --out summary.csv
```

??? example "Terminal output"

    ```text
    Summary written to summary.csv
    ```

## Next steps

- [Filter and select](filtering.md)
- [Export energy and thermochemistry CSV](../tutorials/energy-csv.md)
- [Find fields by scientific property](../reference/model_fields.md)
- [CLI task recipes](cli-recipes.md)
