# SMILES

<!-- format-support:smi -->

| Item | Value |
| ---- | ----- |
| Format ID | `smi` |
| Extensions | `.smi`, `.txt` |
| Read | Yes |
| Write | Yes |
| Registry role | Reader, file writer, frame writer |
| Data level | Coordinates on read; graph required on write |

SMILES record IO with graph-derived charge/multiplicity and canonical SMILES writing.

```python
from molop import AutoParser

frame = AutoParser("molecules.smi")[0][0]
print(frame.to_canonical_SMILES(), frame.coords.shape)
```

??? example "Output"

    For `CCO ethanol`, output is similar to `CCO (3, 3)`.

The trailing name is outside the parsed contract and three-dimensional coordinates are not preserved.

| Feature | Support | Scope | Limits |
| ------- | ------- | ----- | ------ |
| <!-- feature-area:Reader -->Reader | Supported | Non-empty SMILES records; first whitespace token parsed; 2D coordinates, graph, formal charge/radical, charge, and multiplicity derived by RDKit. | Trailing record names or columns are not part of the advertised parsed-data contract; input 3D coordinates cannot be preserved. |
| <!-- feature-area:Writer graph policy -->Writer graph policy | Supported | Canonical graph-derived SMILES writer with strict graph semantics. | Coords-only rendering is not the default writer contract. |
| <!-- feature-area:Registry conversion -->Registry conversion | Supported | Gaussian output structures can be converted to SMILES through the registry. | Conversion depends on successful graph recovery or transformation. |
