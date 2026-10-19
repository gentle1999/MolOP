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
| <!-- feature-area:Reader -->Reading | Supported | Non-empty SMILES records; first whitespace token parsed; 2D coordinates, graph, formal charge/radical, charge, and multiplicity derived by RDKit. | Trailing record names or columns are not part of the advertised parsed-data contract; input 3D coordinates cannot be preserved. |
| <!-- feature-area:Writer graph policy -->Molecular graph and writing | Supported | Canonical graph-derived SMILES writer with strict graph semantics. | Coords-only rendering is not the default writer contract. |
| <!-- feature-area:Registry conversion -->Format conversion | Supported | Gaussian output structures can be converted to SMILES through the registry. | Conversion depends on successful graph recovery or transformation. |

## Property sources and reading conventions

<!-- parsing-basis:start -->

The table maps file contents to MolOP data. See [data containers and cross-QM
semantics](../data_model.md) for shared units, indices and missing values.

References: [RDKit file API](https://www.rdkit.org/docs/source/rdkit.Chem.rdmolfiles.html).

| Property | File marker or data source | Reading or writing convention |
| -------- | -------------------------- | ----------------------------- |
| <!-- parsing-basis:Reader -->Reading | `First whitespace-delimited SMILES token`; `RDKit MolFromSmiles and Compute2DCoords` | Atoms/bonds/formal charges/radicals come from the decoded graph. Coordinates are generated 2D depiction, not input geometry. Trailing tokens are not preserved as molecular properties. |
| <!-- parsing-basis:Writer graph policy -->Molecular graph and writing | `Canonical SMILES from the recovered RDKit/Open Babel graph` | Canonicalization is a backend algorithm, not preservation of input token order. Graph recovery and stereochemistry depend on the available molecule. |
| <!-- parsing-basis:Registry conversion -->Format conversion | `Graph-capable writer contract`; `source molecule` | Graph fields feed canonical SMILES output; QM energies, NMR, spectra and source coordinates are not encoded by this writer. |

<!-- parsing-basis:end -->
