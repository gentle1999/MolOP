# 参考

本区用于精确查询，不按学习顺序展开。需要了解一个概念时先看[快速开始](../getting_started/index.md)或[工作流](../guides/index.md)。

| 需要查询 | 入口 | 内容 |
| --- | --- | --- |
| 如何选择相近工具 | [工具选型](comparison.md) | 比较 MolOP、cclib 和 pymatgen 的定位与使用边界 |
| 配置默认值和全局策略 | [配置](config.md) | 日志、并行度、结构恢复和绘图配置 |
| CLI 参数和命令链规则 | [CLI 命令参考](cli.md) | 全局选项、解析选项、操作和 writer 参数 |
| 某种格式支持什么 | [格式支持概览](format_support.md) | reader/writer 状态、科学字段和限制 |
| 某个格式的细节 | [格式参考](formats/xyz.md) | 格式专用读取、写出和边界说明 |
| 某个科学字段在哪里 | [科学字段索引](model_fields.md) | 按性质查找公共字段 |
| Python 签名和成员 | [Python API](api/index.md) | 核心入口、模块 API 和源码生成参考 |
| 转换或序列化边界 | [行为与边界](behavior/transforms.md) | 精确转换、source evidence 和序列化行为 |

## 参考页的使用方式

- 能力是否存在，以[格式支持概览](format_support.md)和格式页为准。
- API 名称、签名和默认值，以[Python API](api/index.md)及其源码生成页面为准。
- 单次调用优先使用显式参数；需要影响后续调用时才修改[全局配置](config.md)。
