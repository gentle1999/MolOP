# MolOP

MolOP provides one Python API and one command-line tool for reading computational chemistry
outputs, extracting results, filtering batches, and exporting structures.

## Read a final energy in three steps

Install the current version:

```bash
pip install git+https://github.com/gentle1999/MolOP.git
```

Download the documentation's [ORCA water example](../assets/examples/water_mp2.out), save it as
`water_mp2.out`, and run:

```python
from molop import AutoParser

batch = AutoParser("water_mp2.out", n_jobs=1)
frame = batch[0][-1]

print(frame.atoms)
print(frame.energies.total_energy.m_as("hartree"))
```

The last line is approximately `-74.999374598107`. `batch[0]` selects the first file and `[-1]`
selects its final frame.

## Choose a task

| Task | Start here |
| --- | --- |
| Install and verify the environment | [Installation](getting_started/installation.md) |
| Parse a real file for the first time | [5-minute start](getting_started/quickstart.md) |
| Read energies, frequencies, populations, or NMR data | [Read calculation results](guides/results.md) |
| Export a batch summary to CSV | [Batch summaries](guides/batch.md) |
| Select normal, optimized, or transition-state jobs | [Filter and select](guides/filtering.md) |
| Export XYZ, SDF, Gaussian, or ORCA inputs | [Convert and export](guides/conversion.md) |
| Check what a format can provide | [Format overview](reference/format_support.md) |

## Common inputs

MolOP provides dedicated readers or writers for Gaussian log/fchk/input, ORCA output/input, xTB
output, XYZ, SDF/MOL, SMILES, CML, and related formats. Use the
[format overview](reference/format_support.md) and individual format pages for exact fields and
limitations.

MolOP does not run Gaussian, ORCA, or xTB calculations. It reads existing files and turns their
contents into queryable and exportable objects.

## Next step

[Install MolOP](getting_started/installation.md), then follow the
[5-minute start](getting_started/quickstart.md).
