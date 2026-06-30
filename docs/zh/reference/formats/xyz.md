# XYZ

<!-- format-support:xyz -->

| 项目 | 值 |
| ---- | -- |
| 格式 ID | `xyz` |
| 扩展名 | `.xyz` |
| 读取 | 是 |
| 写入 | 是 |
| Registry 角色 | Reader、文件 writer、帧 writer |
| 数据层级 | 坐标 |

标准 XYZ 坐标读写，支持从 comment 行读取电荷和自旋多重度。

| 特性 | 支持程度 | 支持范围 | 明确边界 | 测试证据 |
| ---- | -------- | -------- | -------- | -------- |
| <!-- feature-area:Reader -->Reader | 已支持 | 标准多帧 XYZ；comment 中的电荷和自旋多重度；文件级电荷/自旋多重度从第一帧 finalize。 | 不恢复分子图；非法原子数或不完整帧在解析阶段按格式不匹配处理。 | `tests/test_parsers_smoke.py::test_autoparser_xyz_smoke`<br>`tests/test_io_registry_and_batch_more.py::test_base_file_parser_finalizes_file_charge_and_multiplicity_from_first_frame` |
| <!-- feature-area:Writer and conversion -->Writer and conversion | 已支持 | 文件/帧 XYZ 渲染、显式 comment 覆盖，以及从带坐标的已解析文件通过 registry 转换为 XYZ。 | Writer 使用坐标语义；XYZ 不编码 graph-only 元数据。 | `tests/test_structure_format_features.py::test_xyz_frame_writer_uses_explicit_comment_override`<br>`tests/test_io_convert_registry_smoke.py::test_registry_convert_g16log_smoke`<br>`tests/test_format_transform_output_dir.py::test_format_transform_batch_output_dir_writes_under_requested_dir` |
| <!-- feature-area:Format mismatch handling -->Format mismatch handling | 已支持 | 简单格式的不匹配检测延迟到正常解析阶段，避免额外 IO probe。 | XYZ 不要求独立文件级指纹。 | `tests/test_io_registry_and_batch_more.py::test_simple_coordinate_readers_defer_mismatch_to_parse_phase` |
