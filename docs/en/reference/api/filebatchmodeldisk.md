# FileBatchModelDisk

`FileBatchModelDisk` is the collection returned by `AutoParser`. It keeps successfully parsed files
and provides batch-level filtering, summaries, frame selection, conversion, and image previews.

```python
from molop import AutoParser

batch = AutoParser("water_mp2.out", n_jobs=1)
summary = batch.to_summary_df(brief=False, flatten_columns=True)
print(summary.shape)
```

??? example "Complete table result"

    ```text
    (1, 23)
    ```

    [Notebook 02](../../examples/02-batch-summary-filter-select.ipynb) renders the complete HTML
    DataFrame from this call without selecting a temporary subset of columns.

::: molop.io.FileBatchModelDisk
    options:
      show_docstring_parameters: false
      show_docstring_other_parameters: false
      show_docstring_returns: false
      docstring_style: null
