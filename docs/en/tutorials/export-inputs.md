# Export structures and next-step inputs

Generate XYZ, Gaussian input, and ORCA input from normally terminated final frames.

## Python

```python
from pathlib import Path
from molop import AutoParser

batch = AutoParser("water_mp2.out", n_jobs=1)
selected = batch.filter_state("normal")

for directory in ("xyz", "gjf", "orca"):
    Path(directory).mkdir(exist_ok=True)

selected.format_transform(
    "xyz", output_dir="xyz", write_to_disk=True
)
selected.format_transform(
    "gjf",
    output_dir="gjf",
    write_to_disk=True,
    route_section="#p B3LYP/6-31G(d) opt",
)
selected.format_transform(
    "orcainp",
    output_dir="orca",
    write_to_disk=True,
    keywords="B3LYP def2-SVP Opt",
)
```

## Output

??? example "Created files"

    ```text
    xyz/
      water_mp2.xyz
    gjf/
      water_mp2.gjf
    orca/
      water_mp2.inp
    ```

Each source produces one final-frame file named from its stem. For a batch, replace
`water_mp2.out` with a path or glob. Before submission, check
route/keywords, resources, charge, multiplicity, and solvent. Conversion does not choose scientific
settings for you.

## Export XYZ with the CLI

```bash
molop -q parse "results/*" \
  filter-state --state normal \
  format-transform --format xyz --output-dir xyz
```

??? example "Created file"

    ```text
    xyz/<source stem>.xyz
    ```

Gaussian and ORCA writers have many dynamic options. Python calls with explicit arguments are often
easier to audit in batch research workflows.
