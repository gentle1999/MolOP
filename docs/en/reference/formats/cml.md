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
| <!-- feature-area:File and frame writer -->File and frame writing | Supported | Selected-frame file rendering, negative frame index selection, separate selected-frame blocks, and single-frame rendering through the registry. | No CML reader is registered. |
| <!-- feature-area:Rendering engines and validation -->Rendering backends and format differences | Supported | RDKit-backed default rendering, OpenBabel-backed rendering, frame-less file rejection, and unsupported-engine diagnostics. | Requires a recoverable RDKit or OpenBabel molecule; arbitrary XML/CML source preservation is not covered. |

## Data sources and writing conventions

<!-- parsing-basis:start -->

The table maps file contents to MolOP data. See [data containers and cross-QM
semantics](../data_model.md) for shared units, indices and missing values.

References: [RDKit file API](https://www.rdkit.org/docs/source/rdkit.Chem.rdmolfiles.html); [Open Babel formats](https://openbabel.org/docs/FileFormats/Overview.html).

| Property | File marker or data source | Reading or writing convention |
| -------- | -------------------------- | ----------------------------- |
| <!-- parsing-basis:File and frame writer -->File and frame writing | `Selected graph-capable frame`; `atom/bond/charge/coordinate fields`; `engine XML serializer` | Writer only: graph/coordinates are existing or reconstructed model fields, not read observations. Selection and concatenation follow MolOP writer parameters; no CML reader is registered. |
| <!-- parsing-basis:Rendering engines and validation -->Rendering backends and format differences | `RDKit MolToMrvBlock`; `Open Babel omol.write("cml")` | Record the actual backend: RDKit currently emits Marvin XML through MolToMrvBlock, whereas Open Babel uses cml. Backend output schemas are not assumed identical; graph presence and engine are validated. |

<!-- parsing-basis:end -->
