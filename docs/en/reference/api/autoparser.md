# AutoParser

`AutoParser` is the main entry point for parsing chemical files. It expands paths and glob
patterns, probes extension candidates against file content, and returns a
`FileBatchModelDisk` containing successfully parsed files.

```python
from molop import AutoParser

batch = AutoParser("water_mp2.out", n_jobs=1)
print(len(batch), batch[0].detected_format_id)
```

??? example "Output"

    ```text
    1 orcaout
    ```

Use `parser_detection="orcaout"` when an extension such as `.out` is ambiguous. See
[Parse files](../../guides/parsing.md) for path normalization, frame selection, and failure handling.

::: molop.AutoParser
