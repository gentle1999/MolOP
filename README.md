# MolOP (Molecule OPerator)

[中文](README.zh.md) | [English](README.md)

MolOP is a Python 3.10+ library and command-line tool for computational
chemistry files. It reads Gaussian, ORCA, xTB, and common structure formats,
extracts scientific results, filters batches, and exports structures or
next-step calculation inputs.

## Installation

MolOP is not currently published on PyPI or Conda. Install it from GitHub:

```bash
python -m pip install git+https://github.com/gentle1999/MolOP.git
```

Verify the installation:

```bash
python -c "import molop; print(molop.__version__)"
molop --help
```

## Read one calculation result

Download the
[ORCA water example](https://gentle1999.github.io/MolOP/assets/examples/water_mp2.out)
and save it as `water_mp2.out`:

```python
from molop import AutoParser

batch = AutoParser("water_mp2.out", n_jobs=1)
frame = batch[0][-1]

print(batch[0].detected_format_id)
print(len(frame.atoms), frame.coords.shape)
print(frame.energies.total_energy.m_as("hartree"))
```

Output:

```text
orcaout
3 (3, 3)
-74.999374598107
```

## Common tasks

### Export a batch CSV

```python
batch = AutoParser("results/*.log")
summary = batch.to_summary_df(
    frame=-1,
    brief=False,
    flatten_columns=True,
)
summary.to_csv("summary.csv", index=False)
```

The result is `summary.csv` with one row per successfully parsed file and
unit-bearing columns such as `Energy.total_energy.hartree`.

### Filter and export structures

```bash
molop -q parse "results/*.out" \
  filter-state --state normal \
  format-transform --format xyz --output-dir structures
```

The result is one `.xyz` file per selected source under `structures/`.

## Supported scope

- QM outputs: Gaussian log/fchk, ORCA output, and xTB output.
- QM inputs: Gaussian and ORCA input reading and canonical writing.
- Structure formats: XYZ, SDF/MOL, SMILES, and a CML writer.
- Common results: structures, energies, thermochemistry, vibrations, orbitals,
  atomic populations, dipole/polarizability, NMR, and calculation status,
  depending on the format and printed source content.

See the
[format overview](https://gentle1999.github.io/MolOP/en/reference/format_support/)
for exact reader/writer status and field boundaries.

MolOP does not run quantum chemistry calculations and is not a dedicated
molecular viewer or molecular dynamics engine.

## Documentation

- [5-minute start](https://gentle1999.github.io/MolOP/en/getting_started/quickstart/)
- [Read calculation results](https://gentle1999.github.io/MolOP/en/guides/results/)
- [Batch summaries](https://gentle1999.github.io/MolOP/en/guides/batch/)
- [Filter and select](https://gentle1999.github.io/MolOP/en/guides/filtering/)
- [Convert and export](https://gentle1999.github.io/MolOP/en/guides/conversion/)
- [Contributing](https://gentle1999.github.io/MolOP/en/contributing/)

## Development

```bash
git clone https://github.com/gentle1999/MolOP.git
cd MolOP
uv sync
make check
```

See the documentation site's Developer section for implementation contracts
and quality gates.

## Citation and license

If MolOP helps your research, please cite:

> MolOP (Molecule OPerator), <https://github.com/gentle1999/MolOP>

This project is licensed under the [MIT License](LICENSE).
