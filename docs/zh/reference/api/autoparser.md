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

如果不需要构造成功文件集合，可以使用顶层迭代接口逐条消费结果：

```python
from molop import iter_parse_outcomes

for outcome in iter_parse_outcomes(["water_mp2.out", "missing.out"], n_jobs=2):
    print(outcome.input_index, outcome.file_path, outcome.status)
```

??? example "输出"

    ```text
    0 /absolute/path/water_mp2.out ok
    1 /absolute/path/missing.out missing
    ```

迭代结果按完成顺序产生，并保留 `input_index`；重复输入路径会产生重复结果。消费者提前关闭
迭代器会回收 worker 和结果流。设置 `fail_fast=True` 会在首个非 `ok` 结果处抛出
`BatchParseError`，并取消剩余工作。

::: molop.AutoParser
