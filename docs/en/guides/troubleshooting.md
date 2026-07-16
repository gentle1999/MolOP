# Troubleshooting

Diagnose common installation, format detection, missing-field, and export problems.

## `AutoParser` returns an empty batch

First confirm that the path matches files:

```python
from pathlib import Path

paths = list(Path("results").glob("*.out"))
print(len(paths), paths[:3])
```

If the file exists, retry in one process with an explicit format:

```python
from molop import AutoParser

batch = AutoParser(
    "results/job.out",
    parser_detection="orcaout",
    n_jobs=1,
)
print(len(batch))
```

A successful single-file parse prints `1`. If it remains `0`, inspect the terminal and `molop.log`
for unsupported format, probe, or parser messages.

## Extension and content disagree

An `.out` file can come from ORCA or xTB, and a `.log` file may not be Gaussian. Use the matching
format ID:

```python
AutoParser("orca.out", parser_detection="orcaout")
AutoParser("xtb.out", parser_detection="xtbout")
AutoParser("gaussian.log", parser_detection="g16log")
```

If explicit detection still fails, check for a truncated output, scheduler log, input file, or empty
file.

## A result field is `None`

`None` usually means that the source did not contain the result or that its printed form is not
structured. It does not mean zero.

```python
frame = batch[0][-1]
if frame.vibrations is None:
    print("this frame has no structured frequency data")
```

Check whether:

- you selected the frame that actually contains the result;
- the calculation requested and printed the property;
- the format page declares support for the field;
- `only_extract_structure=True` was enabled accidentally.

## `is_normal` or `is_optimized` is `None`

MolOP preserves unknown status. Inputs, truncated outputs, and segments without termination evidence
may not support a true/false decision. Record unknown separately instead of silently converting it
with `bool(None)`.

## Output directory error

The Python API requires an existing `output_dir`:

```python
from pathlib import Path

Path("output").mkdir(parents=True, exist_ok=True)
batch.format_transform(
    "xyz", output_dir="output", write_to_disk=True
)
```

CLI `--output-dir` creates the directory. The Python API writes only with `write_to_disk=True`.

## SDF or SMILES conversion fails

Graph-level formats require usable topology:

```python
frame = batch[0][-1]
print(frame.rdmol)
print(frame.topology_reconstruction_status)
```

If these show `None` and `failed`, check charge, multiplicity, and geometry, then see
[Structure recovery](structure-recovery.md).

## Convert units

Do not read a bare magnitude and guess its unit:

```python
energy = frame.energies.total_energy
print(energy.m_as("hartree"))
print(energy.m_as("eV"))
```

Hartree-to-eV is a per-particle energy conversion. Converting to `kcal/mol` requires explicit
per-mole semantics and cannot be requested directly from an ordinary electronic-energy quantity.

## Report a reproducible issue

Include:

- `molop --version` and `python --version`;
- the smallest input that reproduces the issue;
- the exact Python code or CLI command;
- `parser_detection`, `n_jobs`, and the error text;
- the expected field and corresponding raw output lines.

Open an issue in [GitHub Issues](https://github.com/gentle1999/MolOP/issues).
