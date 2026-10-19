# OpenBabel Fallback

<!-- format-support:openbabel-fallback -->

| Item | Value |
| ---- | ----- |
| Format ID | `openbabel-fallback` in the support matrix; runtime reader reports `openbabel` |
| Extensions | Unknown extension path, then OpenBabel-supported candidate formats |
| Read | Yes |
| Write | No |
| Registry role | Special fallback reader |
| Data level | Coordinates |

Fallback reader for unknown extensions when OpenBabel can parse the source file.

```python
from molop import AutoParser

parsed = AutoParser("molecule.unknown")[0]
print(parsed.detected_format_id, len(parsed[-1].atoms))
```

??? example "Output contract"

    On success, the reader ID is `openbabel` and the second value is the first molecule's atom
    count.

Available formats depend on the local OpenBabel installation.

| Feature | Support | Scope | Limits |
| ------- | ------- | ----- | ------ |
| <!-- feature-area:Unknown-extension fallback -->Reading unknown extensions | Partial | Unknown extensions can be read through OpenBabel-compatible formats and converted to an XYZFile/XYZFileFrame model with detected format ID persistence and first-molecule coordinate preservation. | Only the first parsed molecule is converted; output is coordinate-level and OpenBabel availability controls format reach. |

## Property sources and reading conventions

<!-- parsing-basis:start -->

The table maps file contents to MolOP data. See [data containers and cross-QM
semantics](../data_model.md) for shared units, indices and missing values.

References: [Open Babel Pybel API](https://openbabel.org/docs/UseTheLibrary/Python_Pybel.html).

| Property | File marker or data source | Reading or writing convention |
| -------- | -------------------------- | ----------------------------- |
| <!-- parsing-basis:Unknown-extension fallback -->Reading unknown extensions | `Open Babel detected format and first pybel.readfile molecule`; `XYZ-level model projection` | Retain detected format ID and the first molecule’s XYZ geometry/charge/spin projection. Only the coordinate-level model fields are transferred; graph and arbitrary properties are not claimed preserved; behavior depends on the installed Open Babel plugin. |

<!-- parsing-basis:end -->
