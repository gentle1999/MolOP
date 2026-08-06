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
| <!-- feature-area:Reader -->Reader | Supported | Standard multi-frame XYZ files; charge and multiplicity in comments; file-level charge/multiplicity finalized from the first frame. | No graph recovery; malformed atom counts or incomplete frames fail during parse-phase mismatch handling. |
| <!-- feature-area:Writer and conversion -->Writer and conversion | Supported | File/frame XYZ rendering, explicit comment override, and registry conversion from coordinate-bearing parsed files. | Writer uses coordinate semantics; graph-only metadata is not encoded by XYZ. |
| <!-- feature-area:Format mismatch handling -->Format mismatch handling | Supported | Simple-format mismatch detection is deferred to normal parsing to avoid a second IO probe. | No separate file-level fingerprint is required for XYZ. |
