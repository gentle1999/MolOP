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
| <!-- feature-area:Unknown-extension fallback -->Unknown-extension fallback | Partial | Unknown extensions can be read through OpenBabel-compatible formats and converted to an XYZFile/XYZFileFrame model with detected format ID persistence and first-molecule coordinate preservation. | Only the first parsed molecule is converted; output is coordinate-level and OpenBabel availability controls format reach. |
