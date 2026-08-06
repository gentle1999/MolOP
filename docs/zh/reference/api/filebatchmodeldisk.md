# FileBatchModelDisk

`FileBatchModelDisk` 是 `AutoParser` 返回的批次容器，保存成功解析的文件，并提供批量过滤、摘要、
frame 选择、格式转换和图像预览操作。

```python
from molop import AutoParser

batch = AutoParser("water_mp2.out", n_jobs=1)
summary = batch.to_summary_df(brief=False, flatten_columns=True)
print(summary.shape)
```

??? example "完整表格结果"

    ```text
    (1, 23)
    ```

    [Notebook 02](../../examples/02-batch-summary-filter-select.ipynb) 会自动渲染这次调用产生的
    完整 HTML DataFrame，不临时截取部分列。

::: molop.io.FileBatchModelDisk
