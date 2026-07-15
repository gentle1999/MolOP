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

Gaussian 输入解析和规范化渲染，覆盖 Link0、route、标题、分子规格、Link1 和附加区段。
writer 可以从已有 GJF frame 渲染，也可以从带坐标的其他格式构建 Gaussian 输入：

```python
rendered = frame.format_transform(
    "gjf",
    link0_commands={"nprocshared": "16", "mem": "32GB"},
    route_section="#p wb97xd/def2tzvp opt freq",
    title_card="optimization and frequency",
    coords_type="cartesian",
    chk=True,
)
```

跨格式转换时应显式提供 Gaussian `route_section`。writer 不会验证计算方法、基组、
任务与 Gaussian 版本之间的兼容性。

| 特性 | 支持程度 | 支持范围 | 明确边界 | 测试证据 |
| ---- | -------- | -------- | -------- | -------- |
| <!-- feature-area:Fixture parse/render inventory -->Fixture parse/render inventory | fixture 覆盖 | 当前维护的 Gaussian 输入 fixtures 可解析、可渲染。 | 覆盖是 fixture-anchored，不等价于完整 Gaussian grammar 证明。 | `tests/test_g16gjf_fixtures_coverage.py::test_g16gjf_all_fixtures_are_parseable_and_renderable` |
| <!-- feature-area:Route semantics -->Route semantics | 部分支持 | 模型化学、任务类型、色散、溶剂化、布居请求和 HF 处理会解析到共享语义容器。 | `functional` 只为 DFT/hybrid/double-hybrid 语义填充；测试外的 Gaussian 方法关键字可能保留为 raw 或非结构化。 | `tests/test_gaussian_route_semantics.py::test_shared_gaussian_route_parser_extracts_model_chemistry_and_capabilities`<br>`tests/test_gaussian_route_semantics.py::test_hf_route_does_not_populate_functional`<br>`tests/test_gaussian_route_semantics.py::test_shared_semantic_route_exposes_structured_dispersion_and_solvation`<br>`tests/test_gaussian_route_semantics.py::test_shared_route_parser_supports_more_gaussian_job_types_and_params` |
| <!-- feature-area:Molecule specification -->Molecule specification | 部分支持 | Cartesian、整数 Cartesian、Z-matrix 变量、free delimiters、fragments 和非法 domain 诊断。 | 支持范围是测试覆盖的显式 molecule-spec grammar 变体。 | `tests/test_gjf_spec_cases.py::test_zmat_labeled_variable_blocks_are_applied_to_atom_lines`<br>`tests/test_gjf_spec_cases.py::test_zmat_unlabeled_variable_blocks_are_applied_by_position`<br>`tests/test_gjf_spec_cases.py::test_free_format_delimiters_are_parsed_and_rendered_canonically`<br>`tests/test_gjf_spec_cases.py::test_integer_cartesian_coordinates_are_not_misparsed_as_zmatrix`<br>`tests/test_gjf_spec_cases.py::test_multifragment_assignments_match_declared_pairs`<br>`tests/test_gjf_spec_cases.py::test_invalid_zmatrix_domains_raise_clear_errors` |
| <!-- feature-area:Additional sections -->Additional sections | 部分支持 | ModRedundant、GIC、NBO、结构化渲染和混合区段诊断。 | 混合 GIC/ModRedundant 区段按设计回退到 unknown，并保留诊断。 | `tests/test_gjf_spec_cases.py::test_modredundant_additional_section_detected_as_typed_section`<br>`tests/test_gjf_spec_cases.py::test_gic_additional_section_detected_as_typed_section`<br>`tests/test_gjf_spec_cases.py::test_nbo_additional_section_detected_as_typed_section`<br>`tests/test_gjf_spec_cases.py::test_mixed_gic_and_modredundant_section_falls_back_to_unknown_with_diagnostic`<br>`tests/test_gjf_spec_cases.py::test_render_prefers_structured_additional_sections` |
| <!-- feature-area:Link1, includes, and writer options -->Link1、include 与 writer 选项 | 部分支持 | Link1 使用原文精确 span，并校验空行；支持 Geom=AllCheck 和转换时 checkpoint 传播。 | `@include` 涉及多个 source artifact，当前单 artifact `SourceSpan` 无法无歧义表达，因此 parser 以 `MOL.PARSE.GJF_INCLUDE_PROVENANCE_UNSUPPORTED` 明确拒绝；并非所有 Gaussian writer 选项都已覆盖。 | `tests/test_gjf_spec_strictness.py::test_link1_requires_blank_line_before_separator`<br>`tests/test_gjf_spec_strictness.py::test_geom_allcheck_allows_missing_title_and_molecule_sections`<br>`tests/test_gjf_include_handling.py::test_gjf_disk_parser_rejects_include_without_multi_artifact_spans`<br>`tests/test_format_transform_output_dir.py::test_format_transform_gjf_chk_propagation_single_file` |

## Gaussian 输入能力覆盖矩阵

下表中的“原文保真”表示能否在解析后保持原始顺序、空白、换行和 token 表示重新输出，
不等同于结构化字段是否保存了同一计算含义。

| Gaussian 输入能力 | 解析 | 结构化语义 | 规范化渲染 | 原文保真 | 当前边界 |
| ----------------- | ---- | ---------- | ---------- | -------- | -------- |
| Link0 `%key=value` | 支持 | 部分 | 支持 | 不支持 | 保存为有序 `GJFLink0Command`；CPU 和 memory 可投影为资源请求，其他 key 不做专用校验。 |
| `%chk` / `%oldchk` 渲染覆盖 | 不适用 | 支持 | 支持 | 不适用 | `chk`、`old_chk` 在临时 Link0 副本末尾追加，不替换或去重已有同名指令。 |
| 单行和多行 route section | 支持 | 部分 | 支持 | 不支持 | route 会归一化为以 `#` 开头的文本；原始折行、空白和 token 大小写不保证保留。 |
| method / functional / basis | 部分 | 部分 | 支持 | 不适用 | 识别 HF、DFT/hybrid/double-hybrid、部分 post-HF 和常见基组；词表外内容仍保留在 route，但可能没有结构化字段。 |
| OPT / FREQ / SP / IRC / scan 等任务 | 部分 | 部分 | 支持 | 不适用 | 常见 job type 和已覆盖参数进入 `task_requests`；Gaussian 全部任务组合未穷举。 |
| `Opt(...)` 选项 | 部分 | 部分 | 支持 | 不适用 | 覆盖 TS、CalcFC、ReadFC、ModRedundant、MaxCycles、coordinate system 等已建模选项。 |
| `Freq(...)` 选项 | 部分 | 部分 | 支持 | 不适用 | 覆盖 anharmonic、projected、hindered rotor、VCD/ROA/Raman、温压和 atom selector 等已识别选项；不执行任务兼容性检查。 |
| TD、SCRF、Pop、Geom 选项 | 部分 | 部分 | 支持 | 不适用 | 已识别选项进入对应语义容器；未识别参数仍在 raw route 中。 |
| title card | 支持 | 支持 | 支持 | 不支持 | 最多五行；模型会清理规定的非法字符，空标题渲染时使用文件名或 `title`。 |
| 电荷与自旋多重度 | 支持 | 支持 | 支持 | 不支持 | 单 fragment 和多 fragment charge/multiplicity 均有结构化模型，输出统一使用规范空白。 |
| Cartesian 坐标 | 支持 | 支持 | 支持 | 不支持 | 整数和浮点坐标均可解析；writer 固定为六位小数，不保留原有效数字和分隔符。 |
| atom label、atom type、原子电荷和参数 | 部分 | 部分 | 部分 | 不支持 | 保存 Gaussian atom specification 的 label、type、charge、`(...)` 参数和 frozen tag；未覆盖所有 ONIOM/MM 扩展语法。 |
| dummy atom 和 ghost atom | 部分 | 部分 | 部分 | 不支持 | `X` 和 `Bq` 形式可投影；包含 dummy/ghost 时拒绝 Cartesian 与 internal 的隐式互转。 |
| 数值 Z-matrix | 支持 | 支持 | 支持 | 不支持 | 支持普通 dihedral 和尾部 `0/1` alternate-angle 标记，并校验长度、角度和引用 domain。 |
| Z-matrix 变量 | 部分 | 部分 | 部分 | 不支持 | 支持已覆盖的 labeled/unlabeled 变量区段；解析时变量被数值替换，writer 不恢复原变量名。 |
| free-format 坐标分隔符 | 部分 | 部分 | 支持 | 不支持 | 已覆盖逗号、tab 和 `/` 等形式；输出统一为空格分隔。 |
| Cartesian/internal 坐标互转 | 不适用 | 支持 | 部分 | 不适用 | `coords_type` 支持普通真实原子 fragment 双向转换；dummy/ghost 或无法安全对齐时抛出 `ValueError`。 |
| fragment molecule specification | 部分 | 支持 | 支持 | 不支持 | 要求每个原子声明连续的 `Fragment=n`，并校验声明的 fragment charge/multiplicity 数量。 |
| `Geom=AllCheck` / `Geom=Checkpoint` | 支持 | 部分 | 部分 | 不支持 | parser 允许缺少 title 或 molecule block；canonical writer 仍按模型字段和默认标题重新组装。 |
| connectivity section | 部分 | 不支持 | 部分 | 不支持 | 原始 connectivity 没有专用字段；`add_gjf_connectivity=True` 会从恢复的分子图重新生成连接表。 |
| ModRedundant | 部分 | 部分 | 支持 | 不支持 | 覆盖已测试的行结构；复杂或未知命令可能作为 unknown/raw additional section。 |
| GIC | 部分 | 部分 | 支持 | 不支持 | 解析 label、options、function、arguments 和 standalone action；不声明完整 Gaussian GIC grammar。 |
| `$NBO ... $END` | 部分 | 部分 | 支持 | 不支持 | 保存 header、commands 和 footer；NBO 命令内容不做进一步语义验证。 |
| Gen/GenECP、自定义基组和其他附加输入 | 部分 | 不支持 | Raw | 部分 | route 可识别 `Gen` / `GenECP`；具体基组/ECP 文本通常作为 unknown additional section 保存和输出。 |
| 混合或未知 additional section | 支持 | 不支持 | Raw | 部分 | 混合 GIC/ModRedundant 会降级为 unknown 并产生 diagnostic；raw 内容可输出但不解释。 |
| `--Link1--` 多任务 | 支持 | 支持 | 支持 | 部分 | 使用精确 source span 分帧并要求 separator 前有空行；文件 writer 使用规范 separator 合并 frame。 |
| Link1 checkpoint 传播 | 不适用 | 部分 | 支持 | 不适用 | 批量/多帧写出可生成 `%chk` / `%oldchk` 链；文件命名和输出目录遵循通用转换参数。 |
| `@include` | 不支持 | 不支持 | 不支持 | 不支持 | 因缺少多 artifact provenance 模型，disk 和 memory parser 都明确拒绝，不递归读取 include。 |
| 从通用坐标 frame 构建 GJF | 不适用 | 部分 | 支持 | 不适用 | 从 `atoms`、`coords`、charge、multiplicity 构建单 fragment Cartesian molecule block；计算 route 应由调用方提供。 |
| Gaussian 版本和关键词兼容性 | 不支持 | 不支持 | 不支持 | 不适用 | `qm_software_version` 记为 `Any`，不判断 route 是否被特定 Gaussian 版本接受。 |

当前 GJF writer 是 **canonical writer**：适合构建和规范化 Gaussian 输入，不是
**source-preserving editor**。需要只修改一个字段且保持其余输入逐字不变时，不应依赖
`parse -> render` 作为无损编辑流程。
