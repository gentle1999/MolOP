# 格式转换与导出

把选定 frame 渲染为文本，或批量写成结构文件和下一步计算输入。

## 先预览，不写盘

```python
from molop import AutoParser

batch = AutoParser("water_mp2.out", n_jobs=1)
rendered = batch.format_transform("xyz", frame=-1, write_to_disk=False)
print(rendered[batch[0].file_path])
```

输出：

```text
3
comment charge 0 multiplicity 1
O               1.7849140000      1.2624220000      0.5119850000
H               2.6482370000      1.0729290000      0.1316310000
H               1.1831680000      1.2568160000     -0.2388350000
```

返回值是 `{源文件绝对路径: 渲染结果}`。`write_to_disk=False` 时不会创建文件。

## 批量写盘

```python
from pathlib import Path

Path("structures").mkdir(exist_ok=True)
batch.format_transform(
    "xyz",
    output_dir="structures",
    frame=-1,
    write_to_disk=True,
)
```

结果目录：

```text
structures/
  water_mp2.xyz
```

`output_dir` 必须已经存在。输出名称保留文件主名并替换最后一个后缀。

## 选择 frame

```python
final_only = batch.format_transform("xyz", frame=-1)
all_in_one = batch.format_transform("xyz", frame="all", embed_in_one_file=True)
one_per_frame = batch.format_transform("xyz", frame="all", embed_in_one_file=False)
```

`embed_in_one_file=False` 返回字符串列表；写盘时生成带 frame 标识的多个文件，具体命名由
writer 决定。

## 可用目标

| 格式 ID | 典型用途 | 主要信息层级 |
| --- | --- | --- |
| `xyz` | 坐标交换、可视化 | 元素、坐标、电荷/多重度注释 |
| `sdf` | 分子图和属性交换 | 原子、键、坐标 |
| `smi` | SMILES 数据集 | 分子图，不保留 3D 坐标 |
| `gjf` | Gaussian 下一步输入 | 坐标、charge/multiplicity、route/link0 |
| `orcainp` | ORCA 下一步输入 | 坐标、关键词、资源和 `%block` |
| `cml` | XML 化学结构交换 | 分子图和坐标 |
| `fakeg` | Gaussian-like 文本 | 用于兼容性展示，不是 Gaussian 原始输出 |

精确 reader/writer 状态见[格式支持概览](../reference/format_support.md)。

## 生成下一步输入

```python
batch.format_transform(
    "gjf",
    output_dir="gaussian_inputs",
    write_to_disk=True,
    route_section="#p B3LYP/6-31G(d) opt",
    link0_commands="%nprocshared=8",
)
```

```python
batch.format_transform(
    "orcainp",
    output_dir="orca_inputs",
    write_to_disk=True,
    keywords="B3LYP def2-SVP Opt",
    maxcore=2000,
)
```

writer 特定参数会随格式变化；CLI 可通过补全或格式页面查看。

## CLI

```bash
mkdir -p structures
molop -q parse "results/*.out" \
  format-transform --format xyz --output-dir structures --frame -1
```

## 信息边界

XYZ 和 SMILES 无法承载完整能量、热化学、振动或 NMR 容器。SDF 可嵌入部分逐原子 QM 属性，
但普通 `format_transform("sdf")` 不等于无损序列化全部计算结果。需要完整数据交换时使用模型
序列化契约，而不是结构格式。

## 下一步

- [导出结构和下一步输入](../tutorials/export-inputs.md)
- [结构恢复](structure-recovery.md)
- [精确转换行为](../reference/transform_behavior.md)
