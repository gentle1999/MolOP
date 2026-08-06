# API Reference

This section provides the auto-generated API reference for the MolOP library, extracted directly from
the Python source code. Use the task guides for workflows; use these pages for signatures, defaults,
and model members.

!!! note "Reference pages are source-driven"
    A field is public only when it appears in the current generated reference. Format-specific
    availability still depends on the source output; use the [format overview](format_support.md)
    and [scientific field index](model_fields.md) alongside the API pages.

## Key Entry Points

- [AutoParser](api/autoparser.md)
- [FileBatchModelDisk](api/filebatchmodeldisk.md)
- [Registry](api/registry.md)

## Modules

- [molop](api/molop.md)
- [molop.io](api/io.md)
- [molop.cli](api/cli.md)
- [molop.config](api/config.md)
- [molop.structure](api/structure.md)
- [molop.unit](api/unit.md)
- [molop.utils](api/utils.md)

## Minimal public workflow

```python
from molop import AutoParser

batch = AutoParser("water_mp2.out", n_jobs=1)
frame = batch[0][-1]
print(frame.energies.total_energy.m_as("hartree"))
```

??? example "Output"

    ```text
    -74.999374598107
    ```

::: molop
    options:
      members: []
