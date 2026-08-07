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

Set `return_report=True` when every supplied path needs an explicit outcome:

```python
from molop import AutoParser, ParseOptions

report = AutoParser(
    ["water_mp2.out", "missing.out"],
    parse_options=ParseOptions(only_last_frame=True),
    return_report=True,
)

for outcome in report.outcomes:
    print(outcome.file_path, outcome.status)
```

??? example "Output"

    ```text
    /absolute/path/missing.out missing
    /absolute/path/water_mp2.out ok
    ```

Outcomes use `ok`, `empty`, `missing`, `skipped`, `unsupported`, `mismatch`, or `error`. Reader
warnings and serializable failure details remain attached to each outcome; `report.batch` retains
the existing successful-file `FileBatchModelDisk`.

::: molop.AutoParser
