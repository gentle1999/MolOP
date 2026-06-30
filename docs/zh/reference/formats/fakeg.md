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

从已解析 Gaussian 输出数据渲染 Gaussian-like 文本。

| 特性 | 支持程度 | 支持范围 | 明确边界 | 测试证据 |
| ---- | -------- | -------- | -------- | -------- |
| <!-- feature-area:Rendered Gaussian-like output sections -->Rendered Gaussian-like output sections | 部分支持 | 从已解析 Gaussian 输出数据渲染 Gaussian-like 结构、电荷/自旋多重度、频率和热力学区段。 | 只是兼容渲染器，不逐字复现 Gaussian log。 | `tests/test_g16log_parser_v3_render.py::test_g16log_v3_render_contains_expected_gaussian_sections`<br>`tests/test_g16log_parser_v3_render.py::test_g16log_file_model_can_render_fakeg_for_all_frames` |
| <!-- feature-area:Frequency and thermochemistry reparse -->Frequency and thermochemistry reparse | 已支持 | 渲染出的频率和热力学内容可以再次解析回帧字段。 | 帧级 fakeg 渲染有意不支持。 | `tests/test_g16log_parser_v3_render.py::test_g16log_file_format_transform_registers_fakeg_file_renderer_only`<br>`tests/test_g16log_parser_v3_render.py::test_g16log_rendered_fakeg_can_be_reparsed_with_frequency_and_thermochemistry` |
