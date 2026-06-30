# Gaussian 输出

<!-- format-support:g16log -->

| 项目 | 值 |
| ---- | -- |
| 格式 ID | `g16log` |
| 扩展名 | `.log`, `.out`, `.g16`, `.gal`, `.irc`, `.gau` |
| 读取 | 是 |
| 写入 | 否 |
| Registry 角色 | Reader |
| 数据层级 | 坐标和 QM 结果 |

Gaussian 输出解析以用户可获得的 QM 结果属性为口径，例如结构、能量、热力学、振动、
轨道、布居、力和响应性质。

| 特性 | 支持程度 | 支持范围 | 明确边界 | 测试证据 |
| ---- | -------- | -------- | -------- | -------- |
| <!-- feature-area:Metadata, geometry, charge, multiplicity, and status -->Metadata, geometry, charge, multiplicity, and status | fixture 覆盖 | 代表性 Gaussian 输出文件暴露原子、坐标或标准取向坐标、电荷、自旋多重度、title/options/keywords、运行时间和终止状态。 | 只读格式；覆盖 fixture 之外的输出区段可能仍是非结构化。 | `tests/test_parsers_smoke.py::test_autoparser_g16log_smoke`<br>`tests/test_g16log_parser_v3_node_payloads.py::test_g16log_v3_first_frame_parses_link1_header_payloads`<br>`tests/test_g16log_numeric_regression.py::test_g16log_file_running_time_is_accumulated_into_model_field`<br>`tests/test_io_registry_and_batch_more.py::test_g16_file_parser_finalizes_status_from_last_frame` |
| <!-- feature-area:Model chemistry and route-derived requests -->Model chemistry and route-derived requests | 部分支持 | 已解析输出帧暴露共享 Gaussian route 语义，包括模型化学、任务类型、色散、溶剂化、布居请求和 HF 处理。 | Route 派生语义受共享 Gaussian route parser 覆盖范围限制；测试外的方法专有关键字可能仍为 raw 或非结构化。 | `tests/test_gaussian_route_semantics.py::test_g16log_frame_exposes_shared_semantic_route` |
| <!-- feature-area:Energies -->Energies | 部分支持 | 代表性输出帧填充 reference/electronic 与方法特异能量；live reference energy 优先于 archive 值；terminal archive-only frame 不会凭空产生能量。 | 只声明测试覆盖的能量字段；不保证每个 Gaussian post-HF 或 correction energy 表都结构化。 | `tests/test_g16log_parser_v3_regression.py::test_g16log_v3_matches_v1_on_representative_fixtures`<br>`tests/test_g16log_parser_v3_regression.py::test_g16log_v3_preserves_live_reference_energy_over_archive_value`<br>`tests/test_g16log_parser_v3_regression.py::test_g16log_v3_does_not_invent_archive_energies_on_terminal_frame` |
| <!-- feature-area:Thermochemistry -->Thermochemistry | 部分支持 | Thermal information 包含 temperature、molecular mass、moments of inertia、rotational symmetry number、rotational/vibrational temperatures、rotational constants、ZPVE、thermal energy、enthalpy、Gibbs free energy、entropy 和 heat capacity。 | 热力学覆盖基于 fixture，主要对应 frequency 类型 Gaussian 输出。 | `tests/test_g16log_parser_v3_node_payloads.py::test_g16log_v3_l716_children_decompose_parent_payloads`<br>`tests/test_g16log_parser_v3_render.py::test_g16log_rendered_fakeg_can_be_reparsed_with_frequency_and_thermochemistry` |
| <!-- feature-area:Vibrations and IR data -->Vibrations and IR data | 部分支持 | 暴露振动频率、约化质量、力常数、IR 强度、逐模式位移向量和虚频标记。 | Raman/VCD 等其他谱学变体不在已声明范围内，除非后续有明确测试。 | `tests/test_g16log_numeric_regression.py::test_default_g16log_parser_handles_concatenated_frequency_values`<br>`tests/test_g16log_parser_v3_node_payloads.py::test_g16log_v3_l716_vibration_mode_children_expose_per_mode_semantics` |
| <!-- feature-area:Molecular orbitals and charge/spin populations -->Molecular orbitals and charge/spin populations | 部分支持 | 在代表性闭壳层和开壳层输出中暴露分子轨道能级/占据数和电荷/自旋布居数据。 | 不声明轨道系数矩阵或每种布居分析形式都已覆盖。 | `tests/test_g16log_mo_regression.py::test_default_g16log_parser_preserves_open_shell_mo_counts`<br>`tests/test_g16log_mo_regression.py::test_default_g16log_parser_handles_concatenated_orbital_energies`<br>`tests/test_g16log_parser_v3_regression.py::test_g16log_v3_matches_v1_on_representative_fixtures`<br>`tests/test_g16log_parser_v3_node_payloads.py::test_g16log_v3_major_components_parse_payloads` |
| <!-- feature-area:Forces, Hessian, optimization status, and polarizability -->Forces, Hessian, optimization status, and polarizability | 部分支持 | 代表性输出保留 forces、Hessian 是否存在、几何优化状态和 polarizability。 | 只声明 fixture 覆盖字段；dipole detail 和扩展 polarizability 子区段不在当前声明范围内。 | `tests/test_g16log_parser_v3_regression.py::test_g16log_v3_matches_v1_on_representative_fixtures`<br>`tests/test_g16log_numeric_regression.py::test_extract_labeled_float_tokens_handles_concatenated_polarizability_values` |
| <!-- feature-area:Link1 metadata propagation -->Link1 metadata propagation | 部分支持 | Link1 section 元数据会传播到后续帧，使多步 Gaussian 任务保留期望的逐帧上下文。 | 这一项仅限用户可见的逐帧元数据。 | `tests/test_g16log_link1_metadata.py::test_g16log_link1_sections_propagate_section_metadata_to_later_frames` |
| <!-- feature-area:Registry conversion -->Registry conversion | 已支持 | 已解析 Gaussian 输出可以通过 registry 转换为坐标格式和 graph 格式。 | 转换取决于结构和 graph 的成功恢复。 | `tests/test_io_convert_registry_smoke.py::test_registry_convert_g16log_smoke` |
