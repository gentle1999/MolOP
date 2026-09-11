# API 参考

本节提供从 Python 源代码直接提取的 MolOP 自动生成 API 参考。任务流程见使用指南；本节用于
查询签名、默认值和模型成员。

!!! note "参考页以源码为准"
    当前生成参考中出现的字段才属于公开 API。具体格式是否填充某个字段仍取决于源文件内容，
    请同时查看[格式支持概览](../format_support.md)和[科学字段索引](../model_fields.md)。

## 关键入口

- [AutoParser](autoparser.md)
- [FileBatchModelDisk](filebatchmodeldisk.md)
- [Registry](../../developer/extensions/registry.md)

## 模块

- [molop](molop.md)
- [molop.io](io.md)
- [molop.cli](cli.md)
- [molop.config](../config.md)
- [molop.structure](structure.md)
- [molop.unit](unit.md)
- [molop.utils](utils.md)

## 最小公共流程

```python
from molop import AutoParser

batch = AutoParser("water_mp2.out", n_jobs=1)
frame = batch[0][-1]
print(frame.energies.total_energy.m_as("hartree"))
```

??? example "输出"

    ```text
    -74.999374598107
    ```

::: molop
    options:
      members: []
