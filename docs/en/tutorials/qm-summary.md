# Summarize Gaussian, ORCA, and xTB

Put outputs from different programs in one batch and export final status and energy through common
fields.

## Input

```text
results/
  gaussian_opt.log
  orca_sp.out
  xtb.out
```

The `.out` extension can match both ORCA and xTB. Automatic mode probes the file contents.

## Python

```python
from molop import AutoParser

batch = AutoParser("results/*", n_jobs=3)
summary = batch.to_summary_df(
    frame=-1,
    brief=False,
    flatten_columns=True,
)

columns = [
    "DiskStorage.FilePath",
    "Calc Parameter.Software",
    "Calc Parameter.Method",
    "Status.IsNormal",
    "Energy.total_energy.hartree",
]
result = summary.reindex(columns=columns)
print(result.to_string(index=False))
result.to_csv("qm_summary.csv", index=False)
```

## Output

`qm_summary.csv` contains one row per successfully parsed file with this schema; values come from
each input:

```text
DiskStorage.FilePath | Calc Parameter.Software | Calc Parameter.Method |
Status.IsNormal | Energy.total_energy.hartree
```

Unavailable fields remain empty. Do not compare `total_energy` values across unrelated methods or
Hamiltonians without scientific normalization.

## Inspect detected programs

```python
for parsed_file in batch:
    print(parsed_file.filename, parsed_file.detected_format_id)
```

Typical output:

```text
gaussian_opt.log g16log
orca_sp.out orcaout
xtb.out xtbout
```

See [Gaussian log](../reference/formats/g16log.md),
[ORCA output](../reference/formats/orcaout.md), and
[xTB output](../reference/formats/xtbout.md) for exact capabilities.
