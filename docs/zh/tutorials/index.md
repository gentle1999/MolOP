# 任务教程

这里的 Markdown 教程覆盖完整操作，不需要先打开 Notebook。

<div class="grid cards" markdown>

-   :material-table-large: **汇总 Gaussian、ORCA 与 xTB**

    混合计算输出目录 → 统一状态与能量表

    [打开教程](qm-summary.md){ .md-button }

-   :material-file-delimited: **导出能量和热化学 CSV**

    Gaussian/ORCA 结果 → 可直接分析的 CSV

    [打开教程](energy-csv.md){ .md-button }

-   :material-filter-check: **筛选优化结果与过渡态**

    优化与频率输出 → 选中路径和摘要

    [打开教程](select-results.md){ .md-button }

-   :material-transit-connection-variant: **分析过渡态**

    含过渡态候选的频率输出 → 虚频与前后体候选

    [打开教程](transition-states.md){ .md-button }

-   :material-file-export: **导出结构和下一步输入**

    已完成计算输出 → XYZ、GJF 或 ORCA input

    [打开教程](export-inputs.md){ .md-button }

</div>

## Notebook 补充材料

??? info "Notebook 的输出由 CI 生成"

    Notebook 是上述任务的可运行交互版本。CI 在构建文档前执行全部代码单元并保存输出；
    `mkdocs-jupyter` 直接渲染这些结果，不手工维护输出单元格。

    - [Gaussian 解析与检查](../examples/01-gaussian-parse-and-inspect.ipynb)
    - [批量汇总、过滤与选择](../examples/02-batch-summary-filter-select.ipynb)
    - [转换与导出](../examples/03-transform-and-export.ipynb)

新增 reader/writer 属于开发者任务，见[插件开发](../developer/plugins.md)。
