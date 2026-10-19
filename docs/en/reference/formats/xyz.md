# XYZ

<!-- format-support:xyz -->

| Item | Value |
| ---- | ----- |
| Format ID | `xyz` |
| Extensions | `.xyz` |
| Read | Yes |
| Write | Yes |
| Registry role | Reader, file writer, frame writer |
| Data level | Coordinates |

Standard XYZ coordinate IO with charge/multiplicity comment support.

```python
from molop import AutoParser

frame = AutoParser("molecule.xyz")[0][0]
print(len(frame.atoms), frame.coords.shape, frame.charge, frame.multiplicity)
```

??? example "Output"

    A three-atom neutral singlet prints `3 (3, 3) 0 1`.

XYZ stores neither energy nor a complete molecular graph.

| Feature | Support | Scope | Limits |
| ------- | ------- | ----- | ------ |
| <!-- feature-area:Reader -->Reading | Supported | Standard multi-frame XYZ files; charge and multiplicity in comments; file-level charge/multiplicity finalized from the first frame. | No graph recovery; malformed atom counts or incomplete frames fail during parse-phase mismatch handling. |
| <!-- feature-area:Writer and conversion -->Writing and conversion | Supported | File/frame XYZ rendering, explicit comment override, and registry conversion from coordinate-bearing parsed files. | Writer uses coordinate semantics; graph-only metadata is not encoded by XYZ. |
| <!-- feature-area:Format mismatch handling -->Frame boundaries and invalid input | Supported | Checks atom counts, coordinate columns and frame completeness while reading. | Incomplete frames raise a format mismatch. |

## Property sources and reading conventions

<!-- parsing-basis:start -->

The table maps file contents to MolOP data. See [data containers and cross-QM
semantics](../data_model.md) for shared units, indices and missing values.

References: [XYZ conventions](https://openbabel.org/docs/FileFormats/XYZ_cartesian_coordinates_format.html).

| Property | File marker or data source | Reading or writing convention |
| -------- | -------------------------- | ----------------------------- |
| <!-- parsing-basis:Reader -->Reading | `Atom-count line`; `comment`; `element + X Y Z`; `Open Babel XYZ reference` | `atoms`, `coords` (Å), `comment`; charge/multiplicity in comments are a MolOP extension, defaulting to 0/1. No bond graph is observed. |
| <!-- parsing-basis:Writer and conversion -->Writing and conversion | `Selected source atoms/coordinates`; `optional comment` | Write atom count, comment and Å coordinates; explicit/stored/generated comment precedence is MolOP policy. Graph-only properties are not encoded. |
| <!-- parsing-basis:Format mismatch handling -->Frame boundaries and invalid input | `Positive integer atom count and exactly N coordinate lines` | Frame boundaries follow the count; incomplete frames raise format mismatch during parsing rather than a separate pre-scan. |

<!-- parsing-basis:end -->
