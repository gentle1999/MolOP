# Convert and export

Render selected frames as text or write batches as structure files and next-step calculation inputs.

## Preview without writing

```python
from molop import AutoParser

batch = AutoParser("water_mp2.out", n_jobs=1)
rendered = batch.format_transform("xyz", frame=-1, write_to_disk=False)
print(rendered[batch[0].file_path])
```

??? example "Output"

    ```text
    3
    comment charge 0 multiplicity 1
    O               1.7849140000      1.2624220000      0.5119850000
    H               2.6482370000      1.0729290000      0.1316310000
    H               1.1831680000      1.2568160000     -0.2388350000
    ```

The result is `{absolute source path: rendered value}`. `write_to_disk=False` creates no file.

## Write a batch

```python
from pathlib import Path

Path("structures").mkdir(exist_ok=True)
batch.format_transform(
    "xyz",
    output_dir="structures",
    frame=-1,
    write_to_disk=True,
)
```

??? example "Created files"

    ```text
    structures/
      water_mp2.xyz
    ```

`output_dir` must already exist. Output names preserve the source stem and replace the final suffix.

## Select frames

```python
final_only = batch.format_transform("xyz", frame=-1)
all_in_one = batch.format_transform("xyz", frame="all", embed_in_one_file=True)
one_per_frame = batch.format_transform("xyz", frame="all", embed_in_one_file=False)
print(len(final_only), len(all_in_one), len(one_per_frame))
```

??? example "Output"

    ```text
    1 1 1
    ```

`embed_in_one_file=False` returns a list of strings. When writing, the writer produces multiple
frame-specific files.

## Available targets

| Format ID | Typical use | Main information level |
| --- | --- | --- |
| `xyz` | Coordinates and visualization | Elements, coordinates, charge/multiplicity comment |
| `sdf` | Molecular graph and property exchange | Atoms, bonds, coordinates |
| `smi` | SMILES datasets | Molecular graph, no 3D coordinates |
| `gjf` | Next Gaussian input | Coordinates, charge/multiplicity, route/link0 |
| `orcainp` | Next ORCA input | Coordinates, keywords, resources, `%block` data |
| `cml` | XML chemical structure exchange | Molecular graph and coordinates |
| `fakeg` | Gaussian-like text | Compatibility display, not an original Gaussian output |

See the [format overview](../reference/format_support.md) for exact reader and writer status.

## Generate next-step inputs

```python
from pathlib import Path

Path("gaussian_inputs").mkdir(exist_ok=True)
batch.format_transform(
    "gjf",
    output_dir="gaussian_inputs",
    write_to_disk=True,
    route_section="#p B3LYP/6-31G(d) opt",
    link0_commands={"nprocshared": "8"},
)
```

```python
Path("orca_inputs").mkdir(exist_ok=True)
batch.format_transform(
    "orcainp",
    output_dir="orca_inputs",
    write_to_disk=True,
    keywords="B3LYP def2-SVP Opt",
    maxcore=2000,
)
print("gaussian_inputs/water_mp2.gjf")
print("orca_inputs/water_mp2.inp")
```

??? example "Created files"

    ```text
    gaussian_inputs/water_mp2.gjf
    orca_inputs/water_mp2.inp
    ```

Writer-specific options vary by format. Use CLI completion or the individual format page.

## CLI

```bash
mkdir -p structures
molop -q parse "results/*.out" \
  format-transform --format xyz --output-dir structures --frame -1
```

??? example "Created file"

    ```text
    structures/water_mp2.xyz
    ```

## Information boundary

XYZ and SMILES cannot carry full energy, thermochemistry, vibration, or NMR containers. SDF can
embed selected atom-level QM properties, but ordinary `format_transform("sdf")` is not a lossless
serialization of every calculation result. Use model serialization for complete data exchange.

## Next steps

- [Export structures and next-step inputs](../tutorials/export-inputs.md)
- [Structure recovery](structure-recovery.md)
- [Exact transform behavior](../reference/transform_behavior.md)
