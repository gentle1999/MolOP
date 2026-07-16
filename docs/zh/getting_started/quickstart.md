# 5 分钟上手

用一份可下载的 ORCA 输出完成“解析、读取最终结构和能量、生成摘要、导出 XYZ”。

## 准备

1. 按[安装指南](installation.md)安装 MolOP。
2. 下载 [water_mp2.out](../../assets/examples/water_mp2.out) 到一个空目录。
3. 在同一目录启动 Python。

样例来自 cclib 的 ORCA 回归数据，许可证和来源见
[样例说明](../../assets/examples/SOURCE.txt)。本页只解析已有输出，不运行 ORCA。

## 解析并读取结果

```python
from molop import AutoParser

batch = AutoParser("water_mp2.out", n_jobs=1)
parsed_file = batch[0]
frame = parsed_file[-1]

print(parsed_file.detected_format_id)
print(len(frame.atoms), frame.coords.shape)
print(frame.energies.total_energy.m_as("hartree"))
```

预期得到：

```text
orcaout
3 (3, 3)
-74.999374598107
```

`AutoParser` 总是返回一个 batch。这里的访问顺序是：

```text
batch -> 第一个输入文件 -> 最后一帧
          batch[0]          [-1]
```

## 生成完整摘要

```python
summary = batch.to_summary_df(
    frame=-1,
    brief=False,
    flatten_columns=True,
)
print(summary[[
    "DiskStorage.FilePath",
    "Status.IsNormal",
    "Energy.total_energy.hartree",
]])
summary.to_csv("summary.csv", index=False)
```

`summary.csv` 会包含一行最后一帧结果。`brief=False` 才会包含能量、热化学和振动等扩展列。

## 导出最终结构

```python
rendered = batch.format_transform("xyz", frame=-1, write_to_disk=False)
print(rendered[parsed_file.file_path])

batch.format_transform(
    "xyz",
    output_dir=".",
    frame=-1,
    write_to_disk=True,
)
```

第一次调用只返回 XYZ 文本；第二次在当前目录写出 `water_mp2.xyz`。

## 换成自己的文件

只需要替换路径：

```python
gaussian = AutoParser("calculation.log")
orca = AutoParser("job.out")
xtb = AutoParser("xtb.out")
many_files = AutoParser("results/*.log")
```

同一扩展名可能有歧义时，显式指定格式：

```python
batch = AutoParser("calculation.out", parser_detection="orcaout")
```

## 下一步

- [Python API 初体验](python-api.md)
- [CLI 初体验](cli.md)
- [读取计算结果](../guides/results.md)
- [格式支持概览](../reference/format_support.md)
