# 核心概念

MolOP 读取已有的计算化学文件，并通过统一对象模型暴露其中的内容。最核心的访问关系是：

```text
输入路径或 glob
        |
   AutoParser
        v
      batch  ->  parsed file  ->  frame
                  batch[0]       [-1]
```

## Batch、File 与 Frame

`AutoParser(...)` 返回一个 batch。batch 是成功解析文件的集合，也是汇总、筛选、分组和批量转换
等集合级操作的入口。

```python
from molop import AutoParser

batch = AutoParser("water_mp2.out", n_jobs=1)
parsed_file = batch[0]
frame = parsed_file[-1]

print(len(batch), len(parsed_file), parsed_file.detected_format_id)
print(frame.charge, frame.multiplicity)
print(len(frame.atoms), frame.coords.shape)
```

??? example "输出"

    ```text
    1 1 orcaout
    0 1
    3 (3, 3)
    ```

三个层级的含义不同：

- **Batch** 保存文件集合，提供集合级操作。筛选返回新的 batch。
- **File** 表示一个源文件，可能包含一个 frame、优化轨迹、多个计算 segment，或这些内容的组合。
- **Frame** 表示一次结构和结果快照。`[-1]` 选择保留的最后一帧，但不保证它就是收敛的或科研上
  最适合使用的结构。

MolOP 负责解析文件，不运行 Gaussian、ORCA 或 xTB 计算。

## 格式检测与源文件事实

自动检测会先用扩展名产生候选，再检查文件内容。扩展名并不唯一：`.out` 可能是 ORCA 或 xTB，
`.log` 也可能来自其他程序。扩展名有歧义时显式指定 `parser_detection`：

```python
batch = AutoParser("calculation.out", parser_detection="orcaout")
```

解析模型保留源文件中确实存在的事实，不会为缺失的科学性质填充默认值。字段缺失可能是因为计算
没有请求该性质、程序没有打印该性质，或当前格式没有暴露它。

## 可选结果容器

读取字段前先检查容器：

```python
frame = AutoParser("water_mp2.out", n_jobs=1)[0][-1]

if frame.energies and frame.energies.total_energy is not None:
    print(frame.energies.total_energy.m_as("hartree"))

if frame.vibrations:
    print(frame.vibrations.num_imaginary)
```

??? example "输出"

    ```text
    -74.999374598107
    ```

水分子样例有能量区段，但没有频率区段，因此第二个分支没有输出。`None` 表示源文件没有结构化
值，不是数值零。

## 结构与拓扑

元素和坐标可以独立于分子图存在。访问 `frame.rdmol` 时，如果源文件带有分子图就返回该图；否则
MolOP 可能从坐标恢复拓扑。恢复是惰性的，状态会记录在 frame 上：

```python
mol = frame.rdmol
print(mol.GetNumAtoms(), mol.GetNumBonds())
print(frame.topology_reconstruction_status)
```

??? example "输出"

    ```text
    3 2
    succeeded
    ```

`failed` 表示没有构建出 RDKit 分子。`suspicious_fallback` 表示得到可用候选，但在高可信场景使用
前仍需复核；金属、自由基、离子对和异常几何尤其需要人工检查。见[结构恢复](guides/structure-recovery.md)。

## 汇总与转换是两类操作

用 `to_summary_df(...)` 把解析结果整理为表格，默认 frame 选择器是 `-1`；需要所有 frame 时使用
`frame="all"`。写 CSV 或使用普通列名时可以设置 `flatten_columns=True`，例如
`Energy.total_energy.hartree`。

用 `format_transform(...)` 把选定结构渲染为目标文件格式。Python API 默认返回渲染内容，只有传入
`write_to_disk=True` 才写文件。目标格式决定信息边界：XYZ 无法承载完整能量和振动容器，SDF 转换也
不是所有计算结果的无损序列化。

## 典型工作流

```text
解析一个或多个源文件
        |
选择 frame 或筛选 batch
        |
读取可选科学结果容器
        |
复核状态与拓扑
        |
汇总、渲染或导出
```

从[5 分钟上手](getting_started/quickstart.md)开始，然后选择 [Python API](getting_started/python-api.md)、
[CLI](getting_started/cli.md) 或具体任务教程。
