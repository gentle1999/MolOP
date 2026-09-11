# API Reference

This section provides the auto-generated API reference for the MolOP library, extracted directly from
the Python source code. Use the task guides for workflows; use these pages for signatures, defaults,
and model members.

!!! note "Reference pages are source-driven"
    A field is public only when it appears in the current generated reference. Format-specific
    availability still depends on the source output; use the [format overview](../format_support.md)
    and [scientific field index](../model_fields.md) alongside the API pages.

## Key Entry Points

- [AutoParser](autoparser.md)
- [FileBatchModelDisk](filebatchmodeldisk.md)
- [Registry](../../developer/extensions/registry.md)

## Modules

- [molop](molop.md)
- [molop.io](io.md)
- [molop.cli](cli.md)
- [molop.config](../config.md)
- [molop.structure](structure.md)
- [molop.unit](unit.md)
- [molop.utils](utils.md)

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
