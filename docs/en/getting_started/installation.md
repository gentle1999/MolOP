# Installation

```bash
pip install molop
```

## Verify

```bash
python -c "import molop; print(molop.__version__)"
molop --version
molop --help
```

`molop --help` should list the `parse` command. An untagged source checkout can show a development
version.

The stable output shape is:

??? example "Verification output"

    ```text
    <version>
    molop, version <version>
    Usage: molop [OPTIONS] COMMAND [ARGS]...
    ...
      parse       Parse files into a FileBatchModelDisk state, then run...
    ```

The two `<version>` values should match. The exact version and full help text are release-dependent.

## Troubleshooting

???+ note "Install from a source checkout"
    Contributors who need an editable checkout should use the [development environment and quality
    gates](../developer/contributing/quality.md). End-user examples assume the `molop` command is on `PATH`.

### Native crash when importing RDKit and Open Babel

Importing Open Babel before RDKit can cause a native segmentation fault. When using MolOP, import it
first. If you import the native packages directly, load RDKit before Open Babel:

```python
import molop
from rdkit import Chem
from openbabel import pybel
```

MolOP initializes RDKit before Open Babel to avoid this native-library conflict.

## Development environment

Source checkout, `uv sync`, and test commands belong to the contributor workflow. See
[Development environment and quality gates](../developer/contributing/quality.md).

## Next step

Use the shared example in the [5-minute start](quickstart.md).
