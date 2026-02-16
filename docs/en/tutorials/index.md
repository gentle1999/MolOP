# Tutorials

This section provides navigation and key takeaways for using MolOP. **The notebooks under `docs/en/examples/` are the canonical baseline**: if any Markdown page conflicts with the notebook behavior/API/output, the notebook wins.

MkDocs is configured with `mkdocs-jupyter.execute: false`, so notebooks are not executed during docs builds; notebooks in the repository should include reproducible outputs.

Recommended reading order:

- [01-Gaussian Parse and Inspect (Notebook)](../examples/01-gaussian-parse-and-inspect.ipynb)
- [02-Batch Summary, Filter and Select (Notebook)](../examples/02-batch-summary-filter-select.ipynb)
- [03-Transform and Export (Notebook)](../examples/03-transform-and-export.ipynb)

## 1. IO Parsing

### Prerequisites

- MolOP installed.
- Computational chemistry output files (e.g., `.log`) or coordinate files (e.g., `.xyz`, `.sdf`).

### Code Snippet

```python
from molop.io import AutoParser
# Parse a single file
files = AutoParser("example.log")
file = files[0]
# Get the last frame
frame = file[-1]
print(f"Energy: {frame.energies.total_energy}")
```

### Expected Output

- `files` is a `FileBatchModelDisk` object.
- `frame` contains fields like `energies`, `coords`, `atoms`.

## 2. Batch Processing

### Prerequisites

- Multiple files of the same type in a directory.

### Code Snippet

```python
from molop.io import AutoParser
# Parallel parsing using wildcards
batch = AutoParser("*.log", n_jobs=-1)
# Filter successfully optimized structures
opt_batch = batch.filter_state("opt")
print(f"Parsed {len(opt_batch)} optimized structures")
```

### Expected Output

- A filtered `FileBatchModelDisk` object.

## 3. Structure Reconstruction

### Prerequisites

- Atomic coordinates only (e.g., from an XYZ file), needing bond inference.

### Code Snippet

```python
from molop.io import AutoParser
batch = AutoParser("molecule.xyz")
# Automatically triggers the reconstruction algorithm
rdmol = batch[0][0].to_rdmol()
```

### Expected Output

- An RDKit `Mol` object with inferred bonds and bond orders.

## 4. Transforms / Format Conversion

### Prerequisites

- A parsed batch object.

### Code Snippet

```python
from molop.io import AutoParser
batch = AutoParser("*.log")
# Batch convert to SDF and save (see Notebook 03 for canonical parameters)
batch.format_transform(format="sdf", output_dir="./output", frameID=-1)
```

### Expected Output

- `.sdf` files generated in the `./output` directory.

## 5. CLI Usage

### Prerequisites

- `molop` command configured in your terminal.

### Command Examples

```bash
# Generate a summary CSV
uv run molop -q summary "tests/test_files/g16log/2-TS1-Opt.log" --out tutorial_ts_summary.csv --format csv --mode frame --frame -1 --n-jobs 1

# Transform molecular files to another format
uv run molop -q transform "tests/test_files/orca/h2_grad_orca.inp" --to sdf --output-dir ./tutorial_transform_out --frame -1 --embed --parser-detection orcainp --n-jobs 1
```

## 6. Plugin/Codec Extension

### Prerequisites

- Need to support a custom private format.

### Code Snippet

Define in `src/molop/io/logic/coords_parsers/MyParser.py`:

```python
def register(registry):
    @registry.reader_factory(format_id="myfmt", extensions={".myfmt"}, priority=0)
    def my_reader_factory():
        return MyReader()

class MyReader:
    format_id = "myfmt"
    extensions = frozenset({".myfmt"})
    priority = 0
    def read(self, path, **kwargs):
        # Your parsing logic here
        ...
```

### Expected Output

- `AutoParser` can now recognize the `.myfmt` extension.

### Notes

- The registration function must be named `register`; see [Core Concepts](../concepts.md).
