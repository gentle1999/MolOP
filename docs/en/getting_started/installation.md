# Installation

Install MolOP and verify both the Python API and command-line entry point.

## Prerequisites

- Python 3.10 or newer.
- Access to install a Python package from GitHub.
- A dedicated virtual environment is recommended to avoid conflicts with an existing RDKit setup.

MolOP is not currently published on PyPI or Conda. End users install it from GitHub.

## Install

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install git+https://github.com/gentle1999/MolOP.git
```

On Windows PowerShell, activate the environment with `.venv\Scripts\Activate.ps1`.

## Verify

```bash
python -c "import molop; print(molop.__version__)"
molop --version
molop --help
```

`molop --help` should list the `parse` command. The version is derived from Git tags; running an
untagged source checkout can show a development version.

## Troubleshooting

### Incompatible Python version

Check `python --version`. `python` and `python -m pip` in the same terminal must point to the active
Python 3.10+ environment.

### RDKit installation fails

RDKit is a MolOP runtime dependency. Use an environment supported by RDKit. If pip has no suitable
wheel, install RDKit in a Conda environment first, then install MolOP from GitHub.

### Is OpenBabel required?

The dedicated Gaussian, ORCA, xTB, XYZ, SDF, and SMILES readers do not require the OpenBabel
fallback. OpenBabel Python bindings are needed only for unknown-extension fallback or when you
explicitly select the OpenBabel rendering engine.

### How do I install a development environment?

Source checkout, `uv sync`, and test commands belong to the contributor workflow. See
[Development environment and quality gates](../developer/quality.md).

## Next step

Use the shared example in the [5-minute start](quickstart.md).
