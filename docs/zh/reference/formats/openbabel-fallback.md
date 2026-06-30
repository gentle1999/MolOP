# OpenBabel Fallback

<!-- format-support:openbabel-fallback -->

| 项目 | 值 |
| ---- | -- |
| 格式 ID | 支持矩阵中为 `openbabel-fallback`；运行时 reader 报告 `openbabel` |
| 扩展名 | 未知扩展名路径，然后尝试 OpenBabel 支持的候选格式 |
| 读取 | 是 |
| 写入 | 否 |
| Registry 角色 | 特殊 fallback reader |
| 数据层级 | 坐标 |

当 OpenBabel 能解析源文件时，用于未知扩展名的 fallback reader。

| 特性 | 支持程度 | 支持范围 | 明确边界 | 测试证据 |
| ---- | -------- | -------- | -------- | -------- |
| <!-- feature-area:Unknown-extension fallback -->Unknown-extension fallback | 部分支持 | 未知扩展名可以通过 OpenBabel 兼容格式读取，并转换为带 detected format ID、且保留首个分子坐标的 `XYZFile`/`XYZFileFrame` 模型。 | 只转换第一个成功解析的分子；输出是坐标层级；实际格式范围取决于 OpenBabel 安装。 | `tests/test_structure_format_features.py::test_openbabel_unknown_extension_fallback_keeps_first_molecule_coordinates`<br>`tests/test_io_convert_registry_smoke.py::test_unknown_extension_fallback_smoke`<br>`tests/test_filebatchmodeldisk_codec_filter.py::test_filter_by_codec_id_uses_detected_format_id` |
