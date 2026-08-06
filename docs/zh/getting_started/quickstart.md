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

??? example "输出"

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
summary.to_csv("summary.csv", index=False)
print(summary.shape)
```

??? example "完整表格结果"

    ```text
    (1, 23)
    ```

    [Notebook 02](../examples/02-batch-summary-filter-select.ipynb) 会把这次调用产生的全部 23 列
    自动渲染为完整 HTML DataFrame，不创建临时列子集。

`summary.csv` 会包含一行最后一帧结果。`brief=False` 才会包含能量、热化学和振动等扩展列。
`DiskStorage.FilePath` 默认保存绝对路径；该列的值随运行目录变化。

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

第一次调用打印：

??? example "输出"

    ```text
    3
    comment charge 0 multiplicity 1
    O               1.7849140000      1.2624220000      0.5119850000
    H               2.6482370000      1.0729290000      0.1316310000
    H               1.1831680000      1.2568160000     -0.2388350000
    ```

第一次调用只返回 XYZ 文本；第二次在当前目录写出 `water_mp2.xyz`。

## 换成自己的文件

将路径替换为你的输入文件或 glob：

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
