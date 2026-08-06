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

::: molop.AutoParser
