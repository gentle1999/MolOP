# Troubleshooting

Start with one representative file and `n_jobs=1`. This keeps parser diagnostics
ordered and separates format problems from parallel execution problems.

## `AutoParser` returns an empty batch

Confirm that the path or glob matches a regular file, then force the expected
format ID for one run:

```python
from pathlib import Path

from molop import AutoParser

path = Path("results/job.out")
print(path.is_file())

batch = AutoParser(path, parser_detection="orcaout", n_jobs=1)
print(len(batch), batch[0].detected_format_id if batch else None)
```

??? example "Expected output for a valid ORCA output"
    ```text
    True
    1 orcaout
    ```

An unmatched glob, missing file, unsupported format ID, or failed content probe
produces no parsed file and emits a diagnostic. If explicit detection still
fails, check whether the file is empty, truncated, an input file, or a scheduler
log rather than a calculation output.

## Extension and content disagree

Extensions such as `.out` and `.log` are ambiguous. Select the reader by its
canonical format ID while diagnosing:

```python
orca = AutoParser("orca.out", parser_detection="orcaout", n_jobs=1)
xtb = AutoParser("xtb.out", parser_detection="xtbout", n_jobs=1)
gaussian = AutoParser("gaussian.log", parser_detection="g16log", n_jobs=1)
```

The [format overview](../reference/format_support.md) lists the registered IDs
and extensions. Explicit detection narrows candidate selection; it does not make
malformed content valid.

## A result field is `None`

`None` means the selected frame has no structured value for that property. It
does not mean zero. Check the frame, the calculation request, the source text,
and the corresponding format page:

```python
frame = batch[0][-1]
if frame.vibrations is None:
    print("no structured frequency data on this frame")
```

??? example "Output when frequency data is absent"
    ```text
    no structured frequency data on this frame
    ```

`only_extract_structure=True` deliberately skips many non-structural fields.
Inputs, truncated outputs, and jobs without termination evidence can also leave
`is_normal` or `is_optimized` as `None`; preserve that unknown state instead of
coercing it with `bool(...)`.

## Writing to an output directory fails

The batch Python API requires an existing directory and writes only when
`write_to_disk=True`:

```python
from pathlib import Path

Path("structures").mkdir(parents=True, exist_ok=True)
batch.format_transform(
    "xyz",
    output_dir="structures",
    write_to_disk=True,
)
```

??? example "Generated files"
    ```text
    structures/job.xyz
    ```

The CLI creates `--output-dir` when needed. With `write_to_disk=False`, Python
returns rendered text and ignores `output_dir`.

## SDF, SMILES, or CML conversion fails

Graph-level formats require a usable molecular graph. Inspect graph recovery
before retrying the writer:

```python
frame = batch[0][-1]
print(frame.rdmol is not None)
print(frame.topology_reconstruction_status)
```

??? example "Output after failed topology recovery"
    ```text
    False
    failed
    ```

Check charge, multiplicity, and geometry before changing reconstruction policy.
See [Structure recovery](structure-recovery.md) for the status contract and
configuration options.

## Native crash when importing RDKit and Open Babel

Importing Open Babel before RDKit can cause a native segmentation fault. Import
MolOP before using the native packages, or load RDKit first in a process that
imports both directly:

```python
import molop
from rdkit import Chem
from openbabel import pybel
```

MolOP initializes RDKit before loading Open Babel-dependent modules. Restart the
Python process after changing the import order; native library state cannot be
repaired after the first conflicting import.

## Report a reproducible issue

Include the MolOP and Python versions, the smallest input that reproduces the
problem, the exact command or Python code, `parser_detection`, `n_jobs`, the
complete error text, and the source lines corresponding to the expected field.
Use [GitHub Issues](https://github.com/gentle1999/MolOP/issues).
