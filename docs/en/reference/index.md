# Reference

This section is for precise lookup rather than sequential learning. Start with
[Quick Start](../getting_started/index.md) or [Workflows](../guides/index.md) when you need an
explanation of a concept.

| Look up | Entry | Contents |
| --- | --- | --- |
| How to choose among related tools | [Tool selection](comparison.md) | Compare the roles and boundaries of MolOP, cclib, and pymatgen |
| Configuration defaults and global policies | [Configuration](config.md) | Logging, parallelism, structure recovery, and drawing |
| CLI parameters and chain rules | [CLI command reference](cli.md) | Global options, parse options, operations, and writer arguments |
| Format capabilities | [Format overview](format_support.md) | Reader/writer status, scientific fields, and limits |
| One format in detail | [Format reference](formats/xyz.md) | Format-specific reading, writing, and boundaries |
| A scientific field | [Scientific field index](model_fields.md) | Find common fields by property |
| Python signatures and members | [Python API](api/index.md) | Core entry points, module APIs, and generated source reference |
| Conversion or serialization boundaries | [Behavior and boundaries](behavior/transforms.md) | Exact transforms, source evidence, and serialization |

## How to use the reference

- Check the [format overview](format_support.md) and individual format pages for capability coverage.
- Use the [Python API](api/index.md) and its generated pages for names, signatures, and defaults.
- Prefer per-call arguments; change the [global configuration](config.md) only for a policy that should affect later calls.
