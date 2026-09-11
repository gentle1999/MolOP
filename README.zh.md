# MolOP (Molecule OPerator)

[中文](https://github.com/gentle1999/MolOP/blob/main/README.zh.md) |
[English](https://github.com/gentle1999/MolOP/blob/main/README.md)

[![PyPI](https://img.shields.io/pypi/v/molop.svg)](https://pypi.org/project/molop/)
[![Python](https://img.shields.io/pypi/pyversions/molop.svg)](https://pypi.org/project/molop/)
[![Typing: Typed](https://img.shields.io/badge/typing-typed-blue.svg)](https://typing.python.org/en/latest/spec/distributing.html#packaging-type-information)
[![Status](https://img.shields.io/pypi/status/molop.svg)](https://pypi.org/project/molop/)
[![Wheel](https://img.shields.io/pypi/wheel/molop.svg)](https://pypi.org/project/molop/)
[![Downloads](https://img.shields.io/pypi/dm/molop.svg)](https://pypi.org/project/molop/)
[![CI](https://github.com/gentle1999/MolOP/actions/workflows/ci.yaml/badge.svg?branch=main)](https://github.com/gentle1999/MolOP/actions/workflows/ci.yaml)
[![Docs](https://github.com/gentle1999/MolOP/actions/workflows/docs-deploy.yml/badge.svg?branch=main)](https://gentle1999.github.io/MolOP/)
[![License](https://img.shields.io/github/license/gentle1999/MolOP.svg)](https://github.com/gentle1999/MolOP/blob/main/LICENSE)
[![Last commit](https://img.shields.io/github/last-commit/gentle1999/MolOP.svg)](https://github.com/gentle1999/MolOP/commits/main/)
[![Issues](https://img.shields.io/github/issues/gentle1999/MolOP.svg)](https://github.com/gentle1999/MolOP/issues)
[![Stars](https://img.shields.io/github/stars/gentle1999/MolOP.svg)](https://github.com/gentle1999/MolOP/stargazers)
[![Forks](https://img.shields.io/github/forks/gentle1999/MolOP.svg)](https://github.com/gentle1999/MolOP/network/members)

MolOP 是面向计算化学文件的 Python 3.10+ 库和命令行工具。它根据文件内容选择已注册 reader，
把不同量化软件和结构格式统一为 batch、file、frame 及科学结果容器，用于批量筛选、汇总和
格式转换。下游处理不需要为每种软件维护独立的数据提取流程。

## 安装

```bash
pip install molop
```

验证安装：

```bash
python -c "import molop; print(molop.__version__)"
molop --help
```

`molop --help` 应列出 `parse` 命令。

<details>
<summary>验证输出形状</summary>

```text
<version>
Usage: molop [OPTIONS] COMMAND [ARGS]...
...
  parse       Parse files into a FileBatchModelDisk state, then run...
```

</details>

## 统一读取不同量化软件文件

以下流程使用 [Gaussian 16 样例](https://gentle1999.github.io/MolOP/assets/examples/mn_complex_sp.log)
和 [ORCA 6 样例](https://gentle1999.github.io/MolOP/assets/examples/water_mp2.out)作为两种代表性
输入。将文件保存到当前目录后运行：

```python
from molop import AutoParser

batch = AutoParser(["mn_complex_sp.log", "water_mp2.out"], n_jobs=1)
print(type(batch).__name__, len(batch))

for parsed_file in batch:
    frame = parsed_file[-1]
    energy = frame.energies.total_energy.m_as("hartree")
    print(parsed_file.detected_format_id, frame.qm_software, frame.method, energy)
```

<details>
<summary>输出</summary>

```text
FileBatchModelDisk 2
g16log Gaussian DFT -2182.472195
orcaout ORCA MP2 -74.999374598107
```

</details>

两个样例进入同一个 `FileBatchModelDisk`，并通过相同的公共字段访问结果。新增 reader 也遵循
这套容器契约；缺失字段保持 `None`，不会用推测值补齐源文件未提供的科学结果。

## 常用任务

### 批量导出 CSV

```python
from molop import AutoParser

batch = AutoParser("water_mp2.out", n_jobs=1)
summary = batch.to_summary_df(
    frame=-1,
    brief=False,
    flatten_columns=True,
)
summary.to_csv("summary.csv", index=False)
print(summary.shape)
```

对随文档提供的 `water_mp2.out` 样例，完整表格包含一行、23 列。处理批量数据时，将输入替换为
路径或 glob：

<details>
<summary>输出</summary>

```text
(1, 23)
```

</details>

结果是每个成功解析文件一行的 `summary.csv`，列名包含单位，例如
`Energy.total_energy.hartree`。[Notebook 02](https://gentle1999.github.io/MolOP/examples/02-batch-summary-filter-select/)
会自动渲染这次调用产生的完整 DataFrame，不临时截取部分列。

### 筛选并导出结构

```bash
molop -q parse "water_mp2.out" --n-jobs 1 \
  filter-state --state normal \
  format-transform --format xyz --output-dir structures
```

命令会生成：

<details>
<summary>生成文件</summary>

```text
structures/water_mp2.xyz
```

</details>

将 `water_mp2.out` 替换为自己的文件路径或 glob，即可处理一批文件。

### 绘制轨迹并导出过渡态候选

解析后的多 frame 文件可以输出 GIF 或动画 SVG；包含单个虚频的计算文件还可以导出基于几何的
前后体候选：

```python
from molop import AutoParser

trajectory = AutoParser("trajectory.log", n_jobs=1)[0]
ts_frame = next(frame for frame in trajectory if frame.is_TS)
trajectory.draw_animation(file_path="trajectory.gif")
pre_path, post_path = ts_frame.save_pre_post_ts("ts-endpoints", format="sdf")
```

详见[过渡态分析教程](https://gentle1999.github.io/MolOP/tutorials/transition-states/)，其中说明了
frame 筛选、振动动图、批量导出以及“候选结构不是优化结构”的边界。

## 支持范围

- QM 输出：Gaussian log/fchk、ORCA output、xTB output。
- QM 输入：Gaussian input、ORCA input 的读取和规范化写出。
- 结构格式：XYZ、SDF/MOL、SMILES，以及 CML writer。
- 公共结果：结构、能量、热化学、振动、轨道、原子布居、偶极/极化率、
  NMR 和计算状态，具体取决于格式和源文件打印内容。
- 可视化：轨迹和振动 GIF/SVG 动图，以及基于几何的过渡态前后体候选。

精确 reader/writer 状态和字段边界见
[格式支持概览](https://gentle1999.github.io/MolOP/reference/format_support/)。

MolOP 不运行量子化学计算，也不是专用分子查看器或分子动力学引擎。

## 文档

- [5 分钟上手](https://gentle1999.github.io/MolOP/getting_started/quickstart/)
- [读取计算结果](https://gentle1999.github.io/MolOP/guides/results/)
- [批量汇总](https://gentle1999.github.io/MolOP/guides/batch/)
- [过滤与选择](https://gentle1999.github.io/MolOP/guides/filtering/)
- [格式转换与导出](https://gentle1999.github.io/MolOP/guides/conversion/)
- [可选的结构恢复与分子图可视化](https://gentle1999.github.io/MolOP/guides/structure-recovery/)
- [过渡态分析与动图](https://gentle1999.github.io/MolOP/tutorials/transition-states/)
- [贡献指南](https://gentle1999.github.io/MolOP/developer/contributing/)

## 开发

```bash
git clone https://github.com/gentle1999/MolOP.git
cd MolOP
uv sync
make check
```

开发约束和测试门禁见文档站点的“开发者”区域。

## 引用与许可证

如果 MolOP 对研究有帮助，请引用：

> MolOP (Molecule OPerator), <https://github.com/gentle1999/MolOP>

本项目采用 [MIT 许可证](https://github.com/gentle1999/MolOP/blob/main/LICENSE)。
