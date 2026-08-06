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

`qm_summary.csv` contains one row per successfully parsed file. The stable output below uses the
bundled single-file ORCA sample:

??? example "Output"

    ```text
    Calc Parameter.Software Calc Parameter.Method Status.IsNormal Energy.total_energy.hartree
                       ORCA                   MP2            True                  -74.999375
    ```

For a mixed directory, the result has one row per successfully parsed input and unavailable fields
remain empty. Environment-dependent columns such as absolute file paths should be handled separately;
do not treat them as fixed output. Do not compare `total_energy` values across unrelated methods or
Hamiltonians without scientific normalization.

To reproduce the displayed row without a mixed directory, use:

```python
sample = AutoParser("water_mp2.out", n_jobs=1)
sample_summary = sample.to_summary_df(
    frame=-1,
    brief=False,
    flatten_columns=True,
)
sample_columns = [
    "Calc Parameter.Software",
    "Calc Parameter.Method",
    "Status.IsNormal",
    "Energy.total_energy.hartree",
]
print(sample_summary.reindex(columns=sample_columns).to_string(index=False))
```

??? example "Output"

    ```text
    Calc Parameter.Software Calc Parameter.Method Status.IsNormal Energy.total_energy.hartree
                       ORCA                   MP2            True                  -74.999375
    ```

## Inspect detected programs

```python
for parsed_file in batch:
    print(parsed_file.filename, parsed_file.detected_format_id)
```

Typical output:

??? example "Output"

    ```text
    gaussian_opt.log g16log
    orca_sp.out orcaout
    xtb.out xtbout
    ```

See [Gaussian log](../reference/formats/g16log.md),
[ORCA output](../reference/formats/orcaout.md), and
[xTB output](../reference/formats/xtbout.md) for exact capabilities.
