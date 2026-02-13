<!--
 * @Author: TMJ
 * @Date: 2026-02-11 16:48:41
 * @LastEditors: TMJ
 * @LastEditTime: 2026-02-11 20:16:24
 * @Description: 请填写简介
-->
# MolOP (Molecule OPerator)

MolOP 是一个专为计算化学工作流设计的 Python 3.10+ 库和命令行工具。它旨在连接原始计算输出与结构化的、可供分析的分子数据。

## 什么是 MolOP？

- **统一解析器**：通过单一接口 (`AutoParser`) 读取 Gaussian log、GJF、XYZ、SDF 等多种格式。
- **结构恢复**：先进的分子图重建算法，能够从坐标中恢复键合信息，对自由基和金属配合物有卓越支持。
- **数据建模**：基于 Pydantic 的模型，提供对能量、振动、轨道和布居分析等数据的类型安全访问。
- **批处理器**：支持数千个文件的并行处理，内置灵活的过滤和格式转换功能。

## 适用场景

- 需要从数百个 Gaussian log 文件中提取热力学数据或分子性质。
- 需要在不同化学文件格式之间转换，同时希望保留或恢复化学键信息。
- 正在构建机器学习流水线，需要从量子化学输出中可靠地提取分子特征。
- 倾向于使用“链式”命令行工具快速检查和处理数据，而无需编写 Python 脚本。

## 非目标

- **量子化学求解器**：MolOP 不执行量子化学计算，它只解析和处理计算结果。
- **可视化工具**：虽然集成了 RDKit，但它不是专门的分子查看器。
- **力场引擎**：它并非为运行分子动力学模拟而设计。

## 开始使用

- [安装指南](getting_started/installation.md)
- [快速上手](getting_started/quickstart.md)
- [命令行界面](command_line_interface.md)

## 引用

如果 MolOP 对您的研究有所帮助，请引用：

> MolOP (Molecule OPerator), <https://github.com/gentle1999/MolOP>
