# Gaussian 输入

<!-- format-support:gjf -->

| 项目 | 值 |
| ---- | -- |
| 格式 ID | `gjf` |
| 扩展名 | `.gjf`, `.gif`, `.com`, `.gau`, `.gjc` |
| 读取 | 是 |
| 写入 | 是 |
| Registry 角色 | Reader、文件 writer、帧 writer |
| 数据层级 | 坐标和 QM 输入语义 |

Gaussian 输入解析/渲染，覆盖 route、分子规格、Link1 和附加区段。

| 特性 | 支持程度 | 支持范围 | 明确边界 | 测试证据 |
| ---- | -------- | -------- | -------- | -------- |
| <!-- feature-area:Fixture parse/render inventory -->Fixture parse/render inventory | fixture 覆盖 | 当前维护的 Gaussian 输入 fixtures 可解析、可渲染。 | 覆盖是 fixture-anchored，不等价于完整 Gaussian grammar 证明。 | `tests/test_g16gjf_fixtures_coverage.py::test_g16gjf_all_fixtures_are_parseable_and_renderable` |
| <!-- feature-area:Route semantics -->Route semantics | 部分支持 | 模型化学、任务类型、色散、溶剂化、布居请求和 HF 处理会解析到共享语义容器。 | `functional` 只为 DFT/hybrid/double-hybrid 语义填充；测试外的 Gaussian 方法关键字可能保留为 raw 或非结构化。 | `tests/test_gaussian_route_semantics.py::test_shared_gaussian_route_parser_extracts_model_chemistry_and_capabilities`<br>`tests/test_gaussian_route_semantics.py::test_hf_route_does_not_populate_functional`<br>`tests/test_gaussian_route_semantics.py::test_shared_semantic_route_exposes_structured_dispersion_and_solvation`<br>`tests/test_gaussian_route_semantics.py::test_shared_route_parser_supports_more_gaussian_job_types_and_params` |
| <!-- feature-area:Molecule specification -->Molecule specification | 部分支持 | Cartesian、整数 Cartesian、Z-matrix 变量、free delimiters、fragments 和非法 domain 诊断。 | 支持范围是测试覆盖的显式 molecule-spec grammar 变体。 | `tests/test_gjf_spec_cases.py::test_zmat_labeled_variable_blocks_are_applied_to_atom_lines`<br>`tests/test_gjf_spec_cases.py::test_zmat_unlabeled_variable_blocks_are_applied_by_position`<br>`tests/test_gjf_spec_cases.py::test_free_format_delimiters_are_parsed_and_rendered_canonically`<br>`tests/test_gjf_spec_cases.py::test_integer_cartesian_coordinates_are_not_misparsed_as_zmatrix`<br>`tests/test_gjf_spec_cases.py::test_multifragment_assignments_match_declared_pairs`<br>`tests/test_gjf_spec_cases.py::test_invalid_zmatrix_domains_raise_clear_errors` |
| <!-- feature-area:Additional sections -->Additional sections | 部分支持 | ModRedundant、GIC、NBO、结构化渲染和混合区段诊断。 | 混合 GIC/ModRedundant 区段按设计回退到 unknown，并保留诊断。 | `tests/test_gjf_spec_cases.py::test_modredundant_additional_section_detected_as_typed_section`<br>`tests/test_gjf_spec_cases.py::test_gic_additional_section_detected_as_typed_section`<br>`tests/test_gjf_spec_cases.py::test_nbo_additional_section_detected_as_typed_section`<br>`tests/test_gjf_spec_cases.py::test_mixed_gic_and_modredundant_section_falls_back_to_unknown_with_diagnostic`<br>`tests/test_gjf_spec_cases.py::test_render_prefers_structured_additional_sections` |
| <!-- feature-area:Link1, includes, and writer options -->Link1, includes, and writer options | 部分支持 | Link1 分割和空行校验；Geom=AllCheck；磁盘 `@include` 展开；转换时 checkpoint 传播。 | `@include` 展开需要源路径上下文；并非所有 Gaussian writer 选项都已覆盖。 | `tests/test_gjf_spec_strictness.py::test_link1_requires_blank_line_before_separator`<br>`tests/test_gjf_spec_strictness.py::test_geom_allcheck_allows_missing_title_and_molecule_sections`<br>`tests/test_gjf_include_handling.py::test_gjf_disk_parser_expands_at_include`<br>`tests/test_format_transform_output_dir.py::test_format_transform_gjf_chk_propagation_single_file` |
