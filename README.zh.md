# MolOP (Molecule OPerator)

[中文](README.zh.md) | [English](README.md)

MolOP 是一个专为计算化学工作流设计的 Python 3.10+ 库和命令行工具。它旨在连接原始计算输出与结构化的、可供分析的分子数据。

## 核心特性

- **统一解析器**：通过单一接口 (`AutoParser`) 读取 Gaussian log、GJF、XYZ、SDF 等多种格式。
- **结构恢复**：先进的分子图重建算法，能够从坐标中恢复键合信息，对自由基和金属配合物有卓越支持。
- **数据建模**：基于 Pydantic 的模型，提供对能量、振动、轨道和布居分析等数据的类型安全访问。
- **批处理器**：支持数千个文件的并行处理，内置灵活的过滤和格式转换功能。

## 适用场景

- 需要从数百个 Gaussian log 文件中提取热力学数据或分子性质。
- 需要在不同化学文件格式之间转换，同时希望保留或恢复化学键信息。
- 正在构建机器学习流水线，需要从量子化学输出中可靠地提取分子特征。
- 倾向于使用命令行工具快速检查和处理数据，而无需编写 Python 脚本。

## 非目标

- **量子化学求解器**：MolOP 不执行量子化学计算，它只解析和处理计算结果。
- **可视化工具**：虽然集成了 RDKit，但它不是专门的分子查看器。
- **力场引擎**：它并非为运行分子动力学模拟而设计。

## 安装

### 面向终端用户

直接从 GitHub 仓库安装最新版本：

```bash
pip install git+https://github.com/gentle1999/MolOP.git
```

### 面向开发人员

克隆仓库并使用 `uv` 同步环境：

```bash
git clone https://github.com/gentle1999/MolOP.git
cd MolOP
uv sync
```

常用开发命令：

- `make test`: 运行测试套件
- `make format`: 格式化代码
- `make check`: 运行所有质量检查（lint, type-check 等）

欢迎通过 [GitHub Issues](https://github.com/gentle1999/MolOP/issues) 或 Pull Requests 参与贡献。

_注意：MolOP 目前未发布在 PyPI 或 Conda。_

## 快速上手

### Python API

使用 `AutoParser` 批量解析文件并生成摘要：

```python
from molop.io import AutoParser

# 解析指定路径下的所有 Gaussian log 文件
batch = AutoParser("path/to/*.log", n_jobs=-1)

# 获取第一个文件的摘要 DataFrame
df = batch[0].to_summary_df()
print(df)
```

### 命令行界面 (CLI)

MolOP 提供现代化的 Typer CLI：

```bash
# 生成摘要 CSV
molop summary "path/to/*.log" -o summary.csv

# 将分子文件转换为另一种格式
molop transform "path/to/*.log" --to sdf --output-dir ./output --frame -1 --embed
```

## 文档与教程

- **官方文档**: [https://gentle1999.github.io/MolOP/](https://gentle1999.github.io/MolOP/)
- **快速入门**: [中文文档](https://gentle1999.github.io/MolOP/getting_started/quickstart/)
- **命令行界面**: [CLI 文档](https://gentle1999.github.io/MolOP/command_line_interface/)
- **示例 Notebooks**:
  - 01-Gaussian 解析与检查: [文档](https://gentle1999.github.io/MolOP/examples/01-gaussian-parse-and-inspect/) | [源码](https://github.com/gentle1999/MolOP/blob/main/docs/zh/examples/01-gaussian-parse-and-inspect.ipynb)
  - 02-批量摘要、过滤与选择: [文档](https://gentle1999.github.io/MolOP/examples/02-batch-summary-filter-select/) | [源码](https://github.com/gentle1999/MolOP/blob/main/docs/zh/examples/02-batch-summary-filter-select.ipynb)
  - 03-格式转换与导出: [文档](https://gentle1999.github.io/MolOP/examples/03-transform-and-export/) | [源码](https://github.com/gentle1999/MolOP/blob/main/docs/zh/examples/03-transform-and-export.ipynb)

## 引用

如果 MolOP 对您的研究有所帮助，请引用：

> MolOP (Molecule OPerator), <https://github.com/gentle1999/MolOP>

## 许可证

本项目采用 [MIT 许可证](LICENSE)。
