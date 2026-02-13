# Installation

## Goal

Install MolOP and its dependencies to start your computational chemistry workflow.

## Prerequisites

- **Python 3.10+**: Ensure you have a compatible Python version installed.
- **RDKit & OpenBabel**: These are required at runtime for full functionality (molecular graph reconstruction, format conversion, etc.).

## Steps

### Option 1: For End-Users (via pip)

Install the latest version directly from the GitHub repository:

```bash
pip install git+https://github.com/gentle1999/MolOP.git
```

### Option 2: For Developers (via uv)

If you want to contribute to MolOP or run the test suite, we recommend using `uv`:

```bash
git clone https://github.com/gentle1999/MolOP.git
cd MolOP
uv sync
```

## Expected Output

After installation, you should be able to run the following verification commands without errors.

### Verify Version

```bash
python -c "import molop; print(molop.__version__)"
```

*Expected: A version string like `0.1.0`.*

### Verify Core Components

```bash
python -c "from molop.io import AutoParser; print(AutoParser)"
```

*Expected: `<function AutoParser at ...>`.*

## Related Links

- [Quickstart](quickstart.md)
- [GitHub Repository](https://github.com/gentle1999/MolOP)
