# 任务教程

这里的 Markdown 教程覆盖完整操作，不需要先打开 Notebook。

| 科研任务 | 输入 | 输出 |
| --- | --- | --- |
| [汇总 Gaussian、ORCA 与 xTB](qm-summary.md) | 混合计算输出目录 | 统一状态与能量表 |
| [导出能量和热化学 CSV](energy-csv.md) | Gaussian/ORCA 结果 | 可直接分析的 CSV |
| [筛选优化结果与过渡态](select-results.md) | 优化与频率输出 | 选中路径和摘要 |
| [导出结构和下一步输入](export-inputs.md) | 已完成计算输出 | XYZ/SDF/GJF/ORCA input |

## Notebook 补充材料

Notebook 用于较长的交互探索，构建文档时不会执行，也不是完成上述任务的前置条件：

- [Gaussian 解析与检查](../examples/01-gaussian-parse-and-inspect.ipynb)
- [批量汇总、过滤与选择](../examples/02-batch-summary-filter-select.ipynb)
- [转换与导出](../examples/03-transform-and-export.ipynb)

新增 reader/writer 属于开发者任务，见[插件开发](../developer/plugins.md)。
