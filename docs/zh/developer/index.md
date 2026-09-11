# 开发者文档

本区面向修改 MolOP 源码、扩展 reader/writer 或维护数据模型的贡献者。

## 责任边界

```text
源文件 -> codec/reader -> File/Frame 公共模型 -> summary/transform/DTO
                         ↘ 格式专用模型
```

- reader 负责从原文提取证据并填充公共或格式专用模型。
- 公共模型定义跨格式可消费的科学结果与单位。
- writer 只承诺目标格式能表达的信息。
- 数据库标识、准入策略和反应路径构建不应塞进 parser。

## 开发入口

| 任务 | 文档 |
| --- | --- |
| 新增 reader/writer | [插件开发](extensions/plugins.md) |
| 查看完整 parser/模型约束 | [Parser 契约](contracts/parser.md) |
| 修改文档 | [文档贡献](contributing/documentation.md) |
| 运行测试和质量门禁 | [开发环境与质量门禁](contributing/quality.md) |
| 理解 QM 公共数据模型 | [QM 公共数据模型](architecture/qm_data_model.md) |
| 理解 ORCA 输入模型 | [ORCA 输入文件模型](architecture/orca_input_model.md) |
| 理解 source lifecycle | [API 契约](contracts/api.md) |

用户任务说明不应依赖本区内容；公共 API 的变化需要同步更新用户指南和可执行示例。
