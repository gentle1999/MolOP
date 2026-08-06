# CML Writer

<!-- format-support:cml -->

| Item | Value |
| ---- | ----- |
| Format ID | `cml` |
| Extensions | `.cml` |
| Read | No |
| Write | Yes |
| Registry role | File writer, frame writer |
| Data level | Graph |

Chemical Markup Language writer backed by RDKit or OpenBabel molecule rendering.

Download the [shared ORCA example](../../../assets/examples/water_mp2.out), then run:

```python
from molop import AutoParser

frame = AutoParser("water_mp2.out", n_jobs=1)[0][-1]
cml_text = frame.format_transform("cml")
print(cml_text.splitlines()[0])
```

??? example "Output"

    ```text
    <?xml version="1.0" encoding="utf-8"?>
    ```

No CML reader is currently registered.

| Feature | Support | Scope | Limits |
| ------- | ------- | ----- | ------ |
| <!-- feature-area:File and frame writer -->File and frame writer | Supported | Selected-frame file rendering, negative frame index selection, separate selected-frame blocks, and single-frame rendering through the registry. | No CML reader is registered. |
| <!-- feature-area:Rendering engines and validation -->Rendering engines and validation | Supported | RDKit-backed default rendering, OpenBabel-backed rendering, frame-less file rejection, and unsupported-engine diagnostics. | Requires a recoverable RDKit or OpenBabel molecule; arbitrary XML/CML source preservation is not covered. |
