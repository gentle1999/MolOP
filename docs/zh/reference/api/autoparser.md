# AutoParser

`AutoParser` 是解析化学文件的主要入口。它会展开路径和 glob，根据扩展名候选及文件内容探测
reader，并返回包含成功解析文件的 `FileBatchModelDisk`。

```python
from molop import AutoParser

batch = AutoParser("water_mp2.out", n_jobs=1)
print(len(batch), batch[0].detected_format_id)
```

??? example "输出"

    ```text
    1 orcaout
    ```

扩展名如 `.out` 存在歧义时使用 `parser_detection="orcaout"`。路径规范化、frame 选择和失败文件
处理见[解析文件](../../guides/parsing.md)。

需要为每个输入路径保留明确结果时，使用 `return_report=True`：

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

??? example "输出"

    ```text
    /absolute/path/missing.out missing
    /absolute/path/water_mp2.out ok
    ```

outcome 状态包括 `ok`、`empty`、`missing`、`skipped`、`unsupported`、`mismatch` 和 `error`。
reader warning 与可序列化失败信息保留在各 outcome 中；`report.batch` 仍是只包含成功文件的
`FileBatchModelDisk`。

::: molop.AutoParser
