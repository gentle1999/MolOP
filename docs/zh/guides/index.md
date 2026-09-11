# 工作流

本区按“要完成什么任务”组织。每页先给出最短可运行示例，再说明常用变体、输出形状和信息边界。

## 日常操作

| 任务 | 页面 | 适合场景 |
| --- | --- | --- |
| 解析文件 | [解析文件](parsing.md) | 单文件、glob、混合输入和解析选项 |
| 读取结果 | [读取计算结果](results.md) | 能量、热化学、振动、轨道、布居和 NMR |
| 批量汇总 | [批量汇总](batch.md) | DataFrame、CSV/JSON 和并行处理 |
| 过滤与选择 | [过滤与选择](filtering.md) | 正常结束、优化结果、过渡态和自定义条件 |
| 转换与导出 | [格式转换与导出](conversion.md) | XYZ、SDF、Gaussian/ORCA 输入和写盘 |
| 结构恢复 | [结构恢复](structure-recovery.md) | 从坐标恢复拓扑、配位键和过渡态候选 |
| CLI 配方 | [CLI 常用任务](cli-recipes.md) | 用 shell 命令复现常见操作链 |
| 排查问题 | [常见问题](troubleshooting.md) | 路径、格式、字段、转换和原生库诊断 |

## 完整任务教程

当任务需要组合多个操作时，进入[完整任务教程](../tutorials/index.md)。教程从输入文件开始，
一路完成筛选、汇总或导出，适合复制后改成自己的工作流。

## 选择阅读方式

- 只想读取一个文件：从[解析文件](parsing.md)开始。
- 已经拿到 batch：先看[读取计算结果](results.md)、[批量汇总](batch.md)或[过滤与选择](filtering.md)。
- 需要写出新文件：看[格式转换与导出](conversion.md)，再按目标格式查阅[格式参考](../reference/index.md)。
- 需要调整日志、并行度或拓扑恢复策略：看[配置参考](../reference/config.md)。
