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
| <!-- feature-area:Reader -->Reading | Supported | SDF/MOL blocks parsed into atom, coordinate, bond, formal charge, radical, total charge, and multiplicity fields. | Arbitrary SD properties are not guaranteed to survive a round trip. |
| <!-- feature-area:Writer graph policy -->Molecular graph and writing | Supported | SDF file/frame writer requires graph-capable input, defaults to strict graph semantics, and supports RDKit/OpenBabel rendering engines. | Coords-only override is rejected when no coords-only SDF writer exists. |
| <!-- feature-area:Registry conversion -->Format conversion | Supported | Parsed structures can be converted to SDF through the codec registry. | Conversion quality depends on recovered or transformed graph data. |

## Property sources and reading conventions

<!-- parsing-basis:start -->

The table maps file contents to MolOP data. See [data containers and cross-QM
semantics](../data_model.md) for shared units, indices and missing values.

References: [RDKit file API](https://www.rdkit.org/docs/source/rdkit.Chem.rdmolfiles.html).

| Property | File marker or data source | Reading or writing convention |
| -------- | -------------------------- | ----------------------------- |
| <!-- parsing-basis:Reader -->Reading | `MOL counts/atom/bond blocks`; `M CHG/M RAD`; `$$$$`; `RDKit SDMolSupplier` | `atoms`, Å `coords`, bonds, formal charges/radicals and title. Charge is summed; multiplicity = radical count + 1 is a MolOP projection, not a measured spin state. Arbitrary SD properties are not mapped. |
| <!-- parsing-basis:Writer graph policy -->Molecular graph and writing | `RDKit MolToMolBlock / Open Babel sdf writer`; `recovered graph` | Serialize a graph-capable molecule; RDKit/Open Babel backends determine atom/bond encoding. A coordinates-only graph reconstruction is a separate MolOP operation. |
| <!-- parsing-basis:Registry conversion -->Format conversion | `Registry graph-level requirements and selected-frame writer` | Conversion reads existing graph/coordinate fields and writer parameters; it does not constitute new source evidence for QM properties. |

<!-- parsing-basis:end -->
