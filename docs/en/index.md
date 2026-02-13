# MolOP (Molecule OPerator)

MolOP is a Python 3.10+ library and CLI designed to streamline computational chemistry workflows. It bridges the gap between raw calculation outputs and structured, analysis-ready molecular data.

## What is MolOP?

- **Unified Parser**: A single interface (`AutoParser`) to read Gaussian logs, GJF, XYZ, SDF, and more.
- **Structure Recovery**: Advanced algorithms to reconstruct molecular graphs (bonds) from coordinates, with superior support for radicals and metal complexes.
- **Data Modeling**: Pydantic-based models providing type-safe access to energies, vibrations, orbitals, and population analysis.
- **Batch Processor**: Parallelized processing of thousands of files with built-in filtering and transformation capabilities.

## When to use

- You need to extract thermodynamic data or properties from hundreds of Gaussian log files.
- You want to convert between chemical file formats while preserving or recovering bonding information.
- You are building a machine learning pipeline and need a reliable way to featurize molecular structures from QM outputs.
- You prefer a "chained" CLI for quick data inspection and manipulation without writing Python scripts.

## Non-goals

- **Quantum Solver**: MolOP does not perform QM calculations; it parses and processes their results.
- **Visualization**: While it integrates with RDKit, it is not a dedicated molecular viewer.
- **Force Field Engine**: It is not designed for running molecular dynamics simulations.

## Getting Started

- [Installation](getting_started/installation.md)
- [Quickstart](getting_started/quickstart.md)
- [Command Line Interface](command_line_interface.md)

## Citation

If MolOP helps your research, please cite it as:
> MolOP (Molecule OPerator), <https://github.com/gentle1999/MolOP>
