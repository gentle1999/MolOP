# Core Concepts

MolOP reads existing computational-chemistry files and exposes their contents through one object
model. The central access pattern is:

```text
input path or glob
        |
   AutoParser
        v
      batch  ->  parsed file  ->  frame
                  batch[0]       [-1]
```

## Batch, file, and frame

`AutoParser(...)` returns a batch. A batch is a collection of successfully parsed files and the
place for collection-level operations such as summaries, filtering, grouping, and conversion.

```python
from molop import AutoParser

batch = AutoParser("water_mp2.out", n_jobs=1)
parsed_file = batch[0]
frame = parsed_file[-1]

print(len(batch), len(parsed_file), parsed_file.detected_format_id)
print(frame.charge, frame.multiplicity)
print(len(frame.atoms), frame.coords.shape)
```

??? example "Output"

    ```text
    1 1 orcaout
    0 1
    3 (3, 3)
    ```

The levels have different meanings:

- A **batch** contains files and supports operations across them. Filters return a new batch.
- A **file** represents one source artifact. It may contain one frame, an optimization trajectory,
  multiple calculation segments, or a combination of these.
- A **frame** is one structure and result snapshot. `[-1]` selects the final retained frame; it does
  not guarantee that the frame is the converged or scientifically preferred structure.

MolOP parses files; it does not run Gaussian, ORCA, or xTB calculations.

## Detection and source facts

Automatic detection uses the extension as a candidate and then inspects file content. Extensions
are not unique: `.out` can be ORCA or xTB, and `.log` can be Gaussian or another program. Use an
explicit `parser_detection` when the source is ambiguous:

```python
batch = AutoParser("calculation.out", parser_detection="orcaout")
```

The parsed model retains facts that are present in the source. It does not fill missing scientific
properties with defaults. A field may be absent because the calculation did not request it, the
program did not print it, or the selected format does not expose it.

## Optional result containers

Check a container before reading its fields:

```python
frame = AutoParser("water_mp2.out", n_jobs=1)[0][-1]

if frame.energies and frame.energies.total_energy is not None:
    print(frame.energies.total_energy.m_as("hartree"))

if frame.vibrations:
    print(frame.vibrations.num_imaginary)
```

??? example "Output"

    ```text
    -74.999374598107
    ```

The water sample has an energy section but no frequency section, so the second branch produces no
line. `None` means that the source did not provide a structured value; it is not a numeric zero.

## Structure and topology

Coordinates and elements are available independently of a molecular graph. Accessing `frame.rdmol`
returns a provided graph when the source contains one; otherwise MolOP may reconstruct topology
from coordinates. The reconstruction is lazy and its status is recorded on the frame:

```python
mol = frame.rdmol
print(mol.GetNumAtoms(), mol.GetNumBonds())
print(frame.topology_reconstruction_status)
```

??? example "Output"

    ```text
    3 2
    succeeded
    ```

`failed` means that no RDKit molecule was built. `suspicious_fallback` is usable as a candidate but
requires scientific review before high-confidence use, especially for metals, radicals, ion pairs,
or distorted geometries. See [Structure recovery](guides/structure-recovery.md).

## Summary and conversion are different operations

Use `to_summary_df(...)` to make a tabular view of parsed results. Its default frame selector is
`-1`; use `frame="all"` for every frame. Use `flatten_columns=True` when writing CSV or working with
ordinary column names such as `Energy.total_energy.hartree`.

Use `format_transform(...)` to render a selected structure into a target file format. Python calls
return rendered content by default; they write files only when `write_to_disk=True`. The target
format defines the information boundary: XYZ cannot carry the full energy and vibration containers,
and an SDF conversion is not a lossless serialization of every calculation result.

## A typical workflow

```text
parse one or more source files
        |
select the frame or filter the batch
        |
read optional scientific containers
        |
review status and topology
        |
summarize, render, or export
```

Start with [5-minute start](getting_started/quickstart.md), then choose the [Python API](getting_started/python-api.md),
[CLI](getting_started/cli.md), or a task tutorial.
