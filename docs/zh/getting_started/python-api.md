# Python API 初体验

掌握 MolOP 用户代码中最常用的 `batch -> file -> frame` 访问方式。

## 解析单个文件

```python
from molop import AutoParser

batch = AutoParser("calculation.log")
parsed_file = batch[0]
final_frame = parsed_file[-1]
```

`len(batch)` 是成功加入 batch 的文件数，`len(parsed_file)` 是该文件的 frame 数。
对共享 `water_mp2.out`，继续运行：

```python
print(len(batch), len(parsed_file), parsed_file.detected_format_id)
print(len(final_frame.atoms), final_frame.coords.shape)
```

输出：

```text
1 1 orcaout
3 (3, 3)
```

## 解析多个输入

```python
batch = AutoParser([
    "reactants/*.log",
    "products/product.out",
    "structures/*.xyz",
])

for parsed_file in batch:
    frame = parsed_file[-1]
    print(parsed_file.filename, frame.charge, frame.multiplicity)
```

每个成功解析文件输出一行，例如：

```text
water_mp2.out 0 1
```

路径、glob 和路径列表可以混用。结果按绝对路径稳定排序，重复路径只解析一次。

## 检查可选结果

不同计算类型提供的字段不同，读取前检查容器是否存在：

```python
frame = batch[0][-1]

if frame.energies and frame.energies.total_energy is not None:
    print(frame.energies.total_energy.m_as("hartree"))

if frame.vibrations:
    print(frame.vibrations.num_imaginary)

if frame.charge_spin_populations:
    names = frame.charge_spin_populations.population_names
    print(names)
```

共享 MP2 样例的能量分支输出 `-74.999374598107`；它没有频率或布居区段，因此另外两个
分支不输出内容。这正是读取可选容器前需要检查的原因。

MolOP 的数值通常带 Pint 单位。使用 `.m_as("目标单位")` 取得指定单位下的数值。

## 生成表格

```python
df = batch.to_summary_df(
    frame=-1,
    brief=False,
    flatten_columns=True,
)
df.to_csv("summary.csv", index=False)
```

使用 `frame="all"` 汇总每个文件的全部 frame；默认 `frame=-1` 只取最后一帧。

## 下一步

- [解析文件](../guides/parsing.md)
- [读取计算结果](../guides/results.md)
- [批量汇总](../guides/batch.md)
