# SDF/MOL

<!-- format-support:sdf -->

| 项目 | 值 |
| ---- | -- |
| 格式 ID | `sdf` |
| 扩展名 | `.sdf`, `.sd`, `.mol` |
| 读取 | 是 |
| 写入 | 是 |
| Registry 角色 | Reader、文件 writer、帧 writer |
| 数据层级 | 读取为坐标；写出要求 graph |

基于 RDKit graph 提取和 strict graph writer 语义的 SDF/MOL 结构读写。

| 特性 | 支持程度 | 支持范围 | 明确边界 | 测试证据 |
| ---- | -------- | -------- | -------- | -------- |
| <!-- feature-area:Reader -->Reader | 已支持 | SDF/MOL block 解析为原子、坐标、键、形式电荷、自由基、总电荷和自旋多重度字段。 | 测试契约不声明任意 SD data field 的 round-trip。 | `tests/test_autoparser_sdf_tmpfile.py::test_autoparser_sdf_tmpfile` |
| <!-- feature-area:Writer graph policy -->Writer graph policy | 已支持 | SDF 文件/帧 writer 要求 graph-capable 输入，默认 strict graph 语义，并支持 RDKit/OpenBabel 渲染后端。 | 当不存在 coords-only SDF writer 时，coords-only override 会被拒绝。 | `tests/test_structure_format_features.py::test_sdf_frame_writer_supports_rdkit_openbabel_engines_and_rejects_unknown_engine`<br>`tests/test_io_registry_and_batch_more.py::test_registry_builtin_default_graph_policies_match_format_semantics`<br>`tests/test_io_registry_and_batch_more.py::test_registry_raises_when_coords_override_has_no_coords_writer` |
| <!-- feature-area:Registry conversion -->Registry conversion | 已支持 | 已解析结构可以通过 codec registry 转换为 SDF。 | 转换质量取决于已恢复或已转换的 graph 数据。 | `tests/test_io_convert_registry_smoke.py::test_registry_convert_xyz_smoke` |
