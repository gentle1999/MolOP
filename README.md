# MolOP (Molecule OPerator)

[中文](README.zh.md) | [English](README.md)

MolOP is a Python 3.10+ library and command-line tool designed for computational chemistry workflows. It aims to bridge the gap between raw computational output and structured, analysis-ready molecular data.

## Core Features

- **Unified Parser**: Read various formats like Gaussian log, GJF, XYZ, SDF through a single interface (`AutoParser`).
- **Structure Recovery**: Advanced molecular graph reconstruction algorithm capable of recovering bonding information from coordinates, with excellent support for radicals and metal complexes.
- **Data Modeling**: Pydantic-based models providing type-safe access to data such as energies, vibrations, orbitals, and population analysis.
- **Batch Processor**: Supports parallel processing of thousands of files with built-in flexible filtering and format conversion.

## Use Cases

- Need to extract thermodynamic data or molecular properties from hundreds of Gaussian log files.
- Need to convert between different chemical file formats while preserving or recovering bond information.
- Building machine learning pipelines and need to reliably extract molecular features from quantum chemistry output.
- Prefer using "chained" command-line tools to quickly inspect and process data without writing Python scripts.

## Non-goals

- **Quantum Chemistry Solver**: MolOP does not perform quantum chemistry calculations; it only parses and processes the results.
- **Visualization Tool**: While integrated with RDKit, it is not a dedicated molecular viewer.
- **Force Field Engine**: It is not designed for running molecular dynamics simulations.

## Installation

### For End Users

Install the latest version directly from the GitHub repository:

```bash
pip install git+https://github.com/gentle1999/MolOP.git
```

### For Developers

Clone the repository and sync the environment using `uv`:

```bash
git clone https://github.com/gentle1999/MolOP.git
cd MolOP
uv sync
```

Common development commands:

- `make test`: Run the test suite
- `make format`: Format the code
- `make check`: Run all quality checks (lint, type-check, etc.)

Contributions via [GitHub Issues](https://github.com/gentle1999/MolOP/issues) or Pull Requests are welcome.

_Note: MolOP is currently not published on PyPI or Conda._

## Quick Start

### Python API

Use `AutoParser` to batch parse files and generate summaries:

```python
from molop.io import AutoParser

# Parse all Gaussian log files in the specified path
batch = AutoParser("path/to/*.log", n_jobs=-1)

# Get the summary DataFrame of the first file
df = batch[0].to_summary_df()
print(df)
```

### Command Line Interface (CLI)

MolOP provides a modern Typer CLI and powerful chained operations:

```bash
# Generate a summary CSV
molop summary "path/to/*.log" -o summary.csv

# Chained operation: Filter transition states and convert to SDF (using legacy chain mode)
molop chain read "path/to/*.log" - filter_state ts - transform sdf --output_dir=. --embed_in_one_file=True - end
```

## Documentation & Tutorials

- **Official Documentation**: [https://gentle1999.github.io/MolOP/](https://gentle1999.github.io/MolOP/)
- **Quick Start**: [English Documentation](https://gentle1999.github.io/MolOP/en/getting_started/quickstart/)
- **Command Line Interface**: [CLI Documentation](https://gentle1999.github.io/MolOP/en/command_line_interface/)
- **Example Notebooks**:
  - 01-Gaussian Parse and Inspect: [Docs](https://gentle1999.github.io/MolOP/en/examples/01-gaussian-parse-and-inspect/) | [Source](https://github.com/gentle1999/MolOP/blob/main/docs/en/examples/01-gaussian-parse-and-inspect.ipynb)
  - 02-Batch Summary, Filter and Select: [Docs](https://gentle1999.github.io/MolOP/en/examples/02-batch-summary-filter-select/) | [Source](https://github.com/gentle1999/MolOP/blob/main/docs/en/examples/02-batch-summary-filter-select.ipynb)
  - 03-Transform and Export: [Docs](https://gentle1999.github.io/MolOP/en/examples/03-transform-and-export/) | [Source](https://github.com/gentle1999/MolOP/blob/main/docs/en/examples/03-transform-and-export.ipynb)

## Citation

If MolOP helps your research, please cite:

> MolOP (Molecule OPerator), <https://github.com/gentle1999/MolOP>

## License

This project is licensed under the [MIT License](LICENSE).
