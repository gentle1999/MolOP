# 5-minute start

Use a downloadable ORCA output to parse a file, read its final structure and energy, build a
summary, and export XYZ.

## Prepare

1. Install MolOP using the [installation guide](installation.md).
2. Download [water_mp2.out](../../assets/examples/water_mp2.out) into an empty directory.
3. Start Python in that directory.

The example comes from cclib's ORCA regression data. See the
[example data notes](../../assets/examples/SOURCE.txt) for its source and license. This tutorial
parses an existing output; it does not run ORCA.

## Parse and read results

```python
from molop import AutoParser

batch = AutoParser("water_mp2.out", n_jobs=1)
parsed_file = batch[0]
frame = parsed_file[-1]

print(parsed_file.detected_format_id)
print(len(frame.atoms), frame.coords.shape)
print(frame.energies.total_energy.m_as("hartree"))
```

Expected output:

```text
orcaout
3 (3, 3)
-74.999374598107
```

`AutoParser` always returns a batch. The access sequence here is:

```text
batch -> first input file -> final frame
          batch[0]           [-1]
```

## Build a full summary

```python
summary = batch.to_summary_df(
    frame=-1,
    brief=False,
    flatten_columns=True,
)
print(summary[[
    "DiskStorage.FilePath",
    "Status.IsNormal",
    "Energy.total_energy.hartree",
]])
summary.to_csv("summary.csv", index=False)
```

`summary.csv` contains one final-frame row. Use `brief=False` to include extended energy,
thermochemistry, and vibration columns.

## Export the final structure

```python
rendered = batch.format_transform("xyz", frame=-1, write_to_disk=False)
print(rendered[parsed_file.file_path])

batch.format_transform(
    "xyz",
    output_dir=".",
    frame=-1,
    write_to_disk=True,
)
```

The first call returns XYZ text. The second writes `water_mp2.xyz` in the current directory.

## Use your own files

Only the path changes:

```python
gaussian = AutoParser("calculation.log")
orca = AutoParser("job.out")
xtb = AutoParser("xtb.out")
many_files = AutoParser("results/*.log")
```

When an extension is ambiguous, specify the format explicitly:

```python
batch = AutoParser("calculation.out", parser_detection="orcaout")
```

## Next steps

- [Python API first steps](python-api.md)
- [CLI first steps](cli.md)
- [Read calculation results](../guides/results.md)
- [Format overview](../reference/format_support.md)
