# Gaussian-Like Renderer

<!-- format-support:fakeg -->

| 项目 | 值 |
| ---- | -- |
| 格式 ID | `fakeg` |
| 扩展名 | `.fakeg` |
| 读取 | 否 |
| 写入 | 是 |
| Registry 角色 | 文件 writer |
| 数据层级 | 已解析 Gaussian 输出数据 |

`fakeg` 将已解析的 Gaussian 输出数据渲染为 Gaussian-like 文本，适合人工检查、
兼容性测试和下游工具衔接。它不是 Gaussian log 的逐字复现器。

`format_transform("fakeg")` 遵循通用转换默认值 `frameID=-1`；需要全文件
Gaussian-like 渲染时传入 `frameID="all"`。

| 能力 | 支持程度 | 可渲染内容 | 边界 |
| ---- | -------- | ---------- | ---- |
| <!-- feature-area:File-level Gaussian-like writer -->文件级 Gaussian-like writer | 部分支持 | 从已解析 Gaussian 输出数据写出 `.fakeg` 文件。 | 只支持文件级渲染；不提供帧级 fakeG writer；选中帧遵循 `frameID`；不保证逐字复现 Gaussian log。 |
| <!-- feature-area:Structure and SCF energy rendering -->结构与 SCF 能量渲染 | 部分支持 | 从坐标和能量字段渲染归一化 orientation 与 SCF-cycle 文本。 | 渲染结果是语义化 Gaussian-like 文本，不是原始 log 文本。 |
| <!-- feature-area:Vibrational frequency rendering -->振动频率渲染 | 部分支持 | frequency、reduced mass、force constant、IR intensity 和逐模式位移区段。 | 当前只声明 frequency 与 IR 相关字段；不声明 Raman/VCD 区段。 |
| <!-- feature-area:Thermochemistry rendering -->热力学渲染 | 部分支持 | temperature、correction、energy、entropy、heat capacity、mass、inertia 和转动/振动 metadata。 | 这是归一化的热力学摘要，不是完整 Gaussian thermochemistry pretty-printer。 |
| <!-- feature-area:Reparseable frequency and thermochemistry output -->可再次解析的频率与热力学输出 | 已支持 | 渲染出的频率和热力学内容可以再次解析回帧字段。 | Round trip 证明的是支持的语义字段，不表示与 Gaussian 原始输出逐字等价。 |
