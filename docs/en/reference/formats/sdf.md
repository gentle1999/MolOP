# SDF/MOL

<!-- format-support:sdf -->

| Item | Value |
| ---- | ----- |
| Format ID | `sdf` |
| Extensions | `.sdf`, `.sd`, `.mol` |
| Read | Yes |
| Write | Yes |
| Registry role | Reader, file writer, frame writer |
| Data level | Coordinates on read; graph required on write |

SDF/MOL structure IO backed by RDKit graph extraction and strict graph writer semantics.

```python
from molop import AutoParser

frame = AutoParser("molecule.sdf")[0][0]
print(frame.rdmol.GetNumAtoms(), frame.rdmol.GetNumBonds())
```

??? example "Output"

    Output is atom and bond counts, for example `3 2` for water.

Arbitrary SD data fields are not guaranteed to round-trip.

| Feature | Support | Scope | Limits |
| ------- | ------- | ----- | ------ |
| <!-- feature-area:Reader -->Reader | Supported | SDF/MOL blocks parsed into atom, coordinate, bond, formal charge, radical, total charge, and multiplicity fields. | Arbitrary SD data-field round-trip is not advertised by the tested contract. |
| <!-- feature-area:Writer graph policy -->Writer graph policy | Supported | SDF file/frame writer requires graph-capable input, defaults to strict graph semantics, and supports RDKit/OpenBabel rendering engines. | Coords-only override is rejected when no coords-only SDF writer exists. |
| <!-- feature-area:Registry conversion -->Registry conversion | Supported | Parsed structures can be converted to SDF through the codec registry. | Conversion quality depends on recovered or transformed graph data. |
