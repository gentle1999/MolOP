# ORCA 输入

<!-- format-support:orcainp -->

| 项目 | 值 |
| ---- | -- |
| 格式 ID | `orcainp` |
| 扩展名 | `.inp` |
| 读取 | 是 |
| 写入 | 否 |
| Registry 角色 | Reader |
| 数据层级 | 坐标和 QM 输入语义 |

ORCA 输入解析，覆盖坐标、多任务分割、模型化学、任务族 fixture 和结构化请求容器。

| 特性 | 支持程度 | 支持范围 | 明确边界 | 测试证据 |
| ---- | -------- | -------- | -------- | -------- |
| <!-- feature-area:Geometry and job splitting -->Geometry and job splitting | 部分支持 | 直接坐标、点电荷、外部 `xyzfile`/`pdbfile` 引用，以及 `$new_job` 分割。 | 外部几何作为引用解析；解析器不要求读取外部文件。 | `tests/test_autoparser_orcainp_tmpfile.py::test_autoparser_orcainp_minimal_xyz_parses_atoms_and_coords`<br>`tests/test_autoparser_orcainp_tmpfile.py::test_autoparser_orcainp_point_charges_parse_into_frame_geometry`<br>`tests/test_autoparser_orcainp_tmpfile.py::test_autoparser_orcainp_xyzfile_does_not_require_external_file`<br>`tests/test_autoparser_orcainp_tmpfile.py::test_autoparser_orcainp_new_job_splits_frames` |
| <!-- feature-area:Model chemistry and options -->Model chemistry and options | 部分支持 | Method、functional、basis、辅助基组、色散作为 functional 后缀、混合基组、print 设置、PARAS 变量和扫描坐标。 | 未支持的 ORCA block 选项可能保留在 raw/resource 字段，而不是专门语义容器中。 | `tests/test_autoparser_orcainp_tmpfile.py::test_autoparser_orcainp_metadata_population_recognizes_d4`<br>`tests/test_autoparser_orcainp_tmpfile.py::test_autoparser_orcainp_mixed_basis_and_output_print_settings`<br>`tests/test_autoparser_orcainp_tmpfile.py::test_autoparser_orcainp_paras_structures_scan_and_resolves_cartesian_variables`<br>`tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_orca6_manual_dispersion_is_functional_suffix` |
| <!-- feature-area:ORCA 6 manual fixture families -->ORCA 6 manual fixture families | fixture 覆盖 | ORCA 6 手册中的单点、SCF stability、优化、频率、激发态、MRCI 和 Solvator 输入示例。 | 缺少完整显式几何的手册片段会被排除。 | `tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_single_point_inputs_include_orca6_manual_fixtures`<br>`tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_scf_stability_inputs_include_orca6_manual_fixtures`<br>`tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_optimization_inputs_include_orca6_manual_fixtures`<br>`tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_frequency_inputs_include_orca6_manual_fixtures`<br>`tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_excited_state_inputs_include_orca6_manual_fixtures`<br>`tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_mrci_inputs_include_orca6_manual_fixtures`<br>`tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_solvator_inputs_include_orca6_manual_fixtures` |
| <!-- feature-area:Excited-state and multireference requests -->Excited-state and multireference requests | 部分支持 | 结构化激发态和多参考任务语义，包括 MRCI multi-job fixture。 | 覆盖示例之外的 ORCA 激发态和多参考关键字族可能仍为 raw。 | `tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_excited_state_manual_fixtures_are_structured`<br>`tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_mrci_manual_fixtures_are_structured`<br>`tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_mrci_multijob_manual_fixture_splits_frames_and_structures_last_job` |
| <!-- feature-area:Optimization, coordinates, frequency, and Solvator structures -->Optimization, coordinates, frequency, and Solvator structures | 部分支持 | 优化约束、内坐标、fragments、NEB、频率 restart 和 Solvator 示例。 | 覆盖范围来自显式 ORCA 手册 fixtures 和针对性回归示例。 | `tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_orca6_optimization_constraints_block_is_not_truncated`<br>`tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_orca6_manual_internal_coords_fill_frame_atoms_and_coords`<br>`tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_orca6_optimization_fragment_mixed_basis_is_structured`<br>`tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_orca6_optimization_neb_block_and_xtb_geometry_are_structured`<br>`tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_orca6_frequency_restart_example_is_structured`<br>`tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_orca6_solvator_fixture_is_structured` |
| <!-- feature-area:Writer availability -->Writer availability | 有意不支持 | ORCA 输入 writer 未注册。 | 在结构化 renderer 能保留 ORCA 输入语义前禁用渲染。 | `tests/test_autoparser_orcainp_tmpfile.py::test_orcainp_writer_is_not_registered` |
