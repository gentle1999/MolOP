# MolOP (Molecule OPerator)

[中文](README.zh.md) | [English](README.md)

MolOP 是面向计算化学文件的 Python 3.10+ 库和命令行工具。它可以读取
Gaussian、ORCA、xTB 和常见结构格式，提取科学结果，批量筛选并导出结构或
下一步计算输入。

## 安装

MolOP 当前未发布到 PyPI 或 Conda，请从 GitHub 安装：

```bash
python -m pip install git+https://github.com/gentle1999/MolOP.git
```

验证安装：

```bash
python -c "import molop; print(molop.__version__)"
molop --help
```

## 读取一个计算结果

从文档下载
[ORCA 水分子样例](https://gentle1999.github.io/MolOP/assets/examples/water_mp2.out)，
保存为 `water_mp2.out`：

```python
from molop import AutoParser

batch = AutoParser("water_mp2.out", n_jobs=1)
frame = batch[0][-1]

print(batch[0].detected_format_id)
print(len(frame.atoms), frame.coords.shape)
print(frame.energies.total_energy.m_as("hartree"))
```

输出：

```text
orcaout
3 (3, 3)
-74.999374598107
```

## 常用任务

### 批量导出 CSV

```python
batch = AutoParser("results/*.log")
summary = batch.to_summary_df(
    frame=-1,
    brief=False,
    flatten_columns=True,
)
summary.to_csv("summary.csv", index=False)
```

结果是每个成功解析文件一行的 `summary.csv`，列名包含单位，例如
`Energy.total_energy.hartree`。

### 筛选并导出结构

```bash
molop -q parse "results/*.out" \
  filter-state --state normal \
  format-transform --format xyz --output-dir structures
```

结果是在 `structures/` 下按源文件主名生成 `.xyz` 文件。

## 支持范围

- QM 输出：Gaussian log/fchk、ORCA output、xTB output。
- QM 输入：Gaussian input、ORCA input 的读取和规范化写出。
- 结构格式：XYZ、SDF/MOL、SMILES，以及 CML writer。
- 公共结果：结构、能量、热化学、振动、轨道、原子布居、偶极/极化率、
  NMR 和计算状态，具体取决于格式和源文件打印内容。

精确 reader/writer 状态和字段边界见
[格式支持概览](https://gentle1999.github.io/MolOP/reference/format_support/)。

MolOP 不运行量子化学计算，也不是专用分子查看器或分子动力学引擎。

## 文档

- [5 分钟上手](https://gentle1999.github.io/MolOP/getting_started/quickstart/)
- [读取计算结果](https://gentle1999.github.io/MolOP/guides/results/)
- [批量汇总](https://gentle1999.github.io/MolOP/guides/batch/)
- [过滤与选择](https://gentle1999.github.io/MolOP/guides/filtering/)
- [格式转换与导出](https://gentle1999.github.io/MolOP/guides/conversion/)
- [贡献指南](https://gentle1999.github.io/MolOP/contributing/)

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

本项目采用 [MIT 许可证](LICENSE)。
