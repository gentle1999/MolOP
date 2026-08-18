# MolOP (Molecule OPerator)

[中文](https://github.com/gentle1999/MolOP/blob/main/README.zh.md) |
[English](https://github.com/gentle1999/MolOP/blob/main/README.md)

[![PyPI](https://img.shields.io/pypi/v/molop.svg)](https://pypi.org/project/molop/)
[![Python](https://img.shields.io/pypi/pyversions/molop.svg)](https://pypi.org/project/molop/)
[![Typing: Typed](https://img.shields.io/badge/typing-typed-blue.svg)](https://typing.python.org/en/latest/spec/distributing.html#packaging-type-information)
[![Status](https://img.shields.io/pypi/status/molop.svg)](https://pypi.org/project/molop/)
[![Wheel](https://img.shields.io/pypi/wheel/molop.svg)](https://pypi.org/project/molop/)
[![Downloads](https://img.shields.io/pypi/dm/molop.svg)](https://pypi.org/project/molop/)
[![CI](https://github.com/gentle1999/MolOP/actions/workflows/ci.yaml/badge.svg?branch=main)](https://github.com/gentle1999/MolOP/actions/workflows/ci.yaml)
[![Docs](https://github.com/gentle1999/MolOP/actions/workflows/docs-deploy.yml/badge.svg?branch=main)](https://gentle1999.github.io/MolOP/)
[![License](https://img.shields.io/github/license/gentle1999/MolOP.svg)](https://github.com/gentle1999/MolOP/blob/main/LICENSE)
[![Last commit](https://img.shields.io/github/last-commit/gentle1999/MolOP.svg)](https://github.com/gentle1999/MolOP/commits/main/)
[![Issues](https://img.shields.io/github/issues/gentle1999/MolOP.svg)](https://github.com/gentle1999/MolOP/issues)
[![Stars](https://img.shields.io/github/stars/gentle1999/MolOP.svg)](https://github.com/gentle1999/MolOP/stargazers)
[![Forks](https://img.shields.io/github/forks/gentle1999/MolOP.svg)](https://github.com/gentle1999/MolOP/network/members)

MolOP is a Python 3.10+ library and command-line tool for computational chemistry files. It selects
a registered reader from file content and maps different quantum-chemistry programs and structure
formats into common batch, file, frame, and scientific-result containers. Downstream processing does
not need a separate extraction path for every program.

## Installation

```bash
pip install molop
```

Verify the installation:

```bash
python -c "import molop; print(molop.__version__)"
molop --help
```

`molop --help` should list the `parse` command.

<details>
<summary>Verification output shape</summary>

```text
<version>
Usage: molop [OPTIONS] COMMAND [ARGS]...
...
  parse       Parse files into a FileBatchModelDisk state, then run...
```

</details>

## Read different quantum-chemistry files through one model

This workflow uses the [Gaussian 16 example](https://gentle1999.github.io/MolOP/assets/examples/mn_complex_sp.log)
and [ORCA 6 example](https://gentle1999.github.io/MolOP/assets/examples/water_mp2.out) as two
representative inputs. Save both files in the current directory, then run:

```python
from molop import AutoParser

batch = AutoParser(["mn_complex_sp.log", "water_mp2.out"], n_jobs=1)
print(type(batch).__name__, len(batch))

for parsed_file in batch:
    frame = parsed_file[-1]
    energy = frame.energies.total_energy.m_as("hartree")
    print(parsed_file.detected_format_id, frame.qm_software, frame.method, energy)
```

<details>
<summary>Output</summary>

```text
FileBatchModelDisk 2
g16log Gaussian DFT -2182.472195
orcaout ORCA MP2 -74.999374598107
```

</details>

Both samples enter one `FileBatchModelDisk` and expose results through the same public fields. New
readers follow the same container contract. Missing fields remain `None`; MolOP does not synthesize
scientific results that the source did not provide.

## Common tasks

### Export a batch CSV

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

For the bundled `water_mp2.out` example, the complete table has one row and 23 columns. Replace the
input with a path or glob for a batch:

<details>
<summary>Output</summary>

```text
(1, 23)
```

</details>

The result is `summary.csv` with one row per successfully parsed file and unit-bearing columns such
as `Energy.total_energy.hartree`. [Notebook 02](https://gentle1999.github.io/MolOP/en/examples/02-batch-summary-filter-select/)
renders the complete DataFrame from this call without selecting a temporary subset of columns.

### Filter and export structures

```bash
molop -q parse "water_mp2.out" --n-jobs 1 \
  filter-state --state normal \
  format-transform --format xyz --output-dir structures
```

The command creates:

<details>
<summary>Created files</summary>

```text
structures/water_mp2.xyz
```

</details>

Replace `water_mp2.out` with a path or glob for your own batch.

### Render trajectories and transition-state candidates

Parsed multi-frame files can be rendered as GIF or animated SVG, and calculation files with one
imaginary mode can export geometry-based pre- and post-TS candidates:

```python
from molop import AutoParser

trajectory = AutoParser("trajectory.log", n_jobs=1)[0]
ts_frame = next(frame for frame in trajectory if frame.is_TS)
trajectory.draw_animation(file_path="trajectory.gif")
pre_path, post_path = ts_frame.save_pre_post_ts("ts-endpoints", format="sdf")
```

See the [transition-state analysis tutorial](https://gentle1999.github.io/MolOP/en/tutorials/transition-states/)
for frame selection, vibration animation, batch export, and the non-optimized endpoint boundary.

## Supported scope

- QM outputs: Gaussian log/fchk, ORCA output, and xTB output.
- QM inputs: Gaussian and ORCA input reading and canonical writing.
- Structure formats: XYZ, SDF/MOL, SMILES, and a CML writer.
- Common results: structures, energies, thermochemistry, vibrations, orbitals,
  atomic populations, dipole/polarizability, NMR, and calculation status,
  depending on the format and printed source content.
- Visualization: trajectory and vibration GIF/SVG animations, plus geometry-based transition-state
  endpoint candidates.

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
- [Optional structure recovery and graph visualization](https://gentle1999.github.io/MolOP/en/guides/structure-recovery/)
- [Transition-state analysis and animations](https://gentle1999.github.io/MolOP/en/tutorials/transition-states/)
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

This project is licensed under the
[MIT License](https://github.com/gentle1999/MolOP/blob/main/LICENSE).
