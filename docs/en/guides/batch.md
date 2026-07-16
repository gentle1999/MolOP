# Batch summaries

Turn a batch of calculation outputs into a pandas DataFrame and export CSV or JSON.

## Shortest example

```python
from molop import AutoParser

batch = AutoParser("results/*.out", n_jobs=4)
summary = batch.to_summary_df(
    frame=-1,
    brief=False,
    flatten_columns=True,
)
summary.to_csv("summary.csv", index=False)
print(summary.shape)
```

If `results/` contains 100 successfully parsed single-job outputs, the shape is normally
`(100, column count)`. Inputs that do not parse produce no row.

## Check the shared example

```python
from molop import AutoParser

batch = AutoParser("water_mp2.out", n_jobs=1)
summary = batch.to_summary_df(brief=False, flatten_columns=True)
print(summary[["Status.IsNormal", "Energy.total_energy.hartree"]])
```

Output:

```text
   Status.IsNormal  Energy.total_energy.hartree
0             True                    -74.999375
```

## Final and all frames

```python
final_rows = batch.to_summary_df(frame=-1)
all_rows = batch.to_summary_df(frame="all")
selected_rows = batch.to_summary_df(frame=[0, -1])
```

The default `on_missing_frame="skip"` omits unavailable indices. Use
`on_missing_frame="error"` when a missing frame should fail the operation.

## Brief and full fields

```python
brief = batch.to_summary_df(brief=True, flatten_columns=True)
full = batch.to_summary_df(brief=False, flatten_columns=True)
print(brief.columns.tolist())
```

The brief summary includes stable file, molecule, calculation, and status columns such as:

```text
['DiskStorage.FilePath', ..., 'Status.IsNormal', 'Status.IsTS',
 'Status.IsOptimized']
```

The full summary adds available `Energy.*`, `Thermal.*`, and `Vibration.*` columns.

## Flat and MultiIndex columns

```python
flat = batch.to_summary_df(flatten_columns=True)
print(flat["General.Charge"])

multi = batch.to_summary_df(flatten_columns=False)
print(multi[("General", "Charge", "")])
```

Flat columns are easier for CSV, JSON, and general analysis. Keep MultiIndex columns when the
three-level group, field, and unit meaning matters.

## Export JSON

```python
summary.to_json("summary.json", orient="records", indent=2)
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

## Next steps

- [Filter and select](filtering.md)
- [Export energy and thermochemistry CSV](../tutorials/energy-csv.md)
- [Find fields by scientific property](../reference/model_fields.md)
