# ORCA 输入

<!-- format-support:orcainp -->

| 项目 | 值 |
| ---- | -- |
| 格式 ID | `orcainp` |
| 扩展名 | `.inp` |
| 读取 | 是 |
| 写入 | 是 |
| Registry 角色 | Reader、file writer、frame writer |
| 数据层级 | 坐标和 QM 输入语义 |

ORCA 输入解析，覆盖坐标、多任务分割、模型化学、任务族 fixture 和结构化请求容器。

带坐标的 frame 可以通过显式 ORCA 计算设置生成输入文件：

```python
rendered = frame.format_transform(
    "orcainp",
    keywords="wB97M-V def2-TZVPP RIJCOSX def2/J TightSCF DefGrid3 NoAutoStart SP",
    nprocs=16,
    maxcore=4000,
    blocks={"scf": {"MaxIter": 300, "STABPerform": True}},
)
```

`rendered` 是包含 ORCA simple input line、`%pal/%maxcore/%scf` 和 `* xyz` 坐标块的字符串。
例如第一行以 `! wB97M-V def2-TZVPP` 开始。

跨格式转换必须提供 `keywords`；MolOP 不会复用其他 QM 软件的关键字语法。
`blocks` 也接受原始 `%block` 文本。

| 特性 | 支持程度 | 支持范围 | 明确边界 | 测试证据 |
| ---- | -------- | -------- | -------- | -------- |
| <!-- feature-area:Geometry and job splitting -->Geometry and job splitting | 部分支持 | 直接坐标、点电荷、外部 `xyzfile`/`pdbfile` 引用，以及 `$new_job` 分割。 | 外部几何作为引用解析；解析器不要求读取外部文件。 | `tests/test_autoparser_orcainp_tmpfile.py::test_autoparser_orcainp_minimal_xyz_parses_atoms_and_coords`<br>`tests/test_autoparser_orcainp_tmpfile.py::test_autoparser_orcainp_point_charges_parse_into_frame_geometry`<br>`tests/test_autoparser_orcainp_tmpfile.py::test_autoparser_orcainp_xyzfile_does_not_require_external_file`<br>`tests/test_autoparser_orcainp_tmpfile.py::test_autoparser_orcainp_new_job_splits_frames` |
| <!-- feature-area:Model chemistry and options -->Model chemistry and options | 部分支持 | Method、functional、basis、辅助基组、色散作为 functional 后缀、混合基组、print 设置、PARAS 变量和扫描坐标。 | 未支持的 ORCA block 选项可能保留在 raw/resource 字段，而不是专门语义容器中。 | `tests/test_autoparser_orcainp_tmpfile.py::test_autoparser_orcainp_metadata_population_recognizes_d4`<br>`tests/test_autoparser_orcainp_tmpfile.py::test_autoparser_orcainp_mixed_basis_and_output_print_settings`<br>`tests/test_autoparser_orcainp_tmpfile.py::test_autoparser_orcainp_paras_structures_scan_and_resolves_cartesian_variables`<br>`tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_orca6_manual_dispersion_is_functional_suffix` |
| <!-- feature-area:ORCA 6 manual fixture families -->ORCA 6 manual fixture families | fixture 覆盖 | ORCA 6 手册中的单点、SCF stability、优化、频率、激发态、MRCI 和 Solvator 输入示例。 | 缺少完整显式几何的手册片段会被排除。 | `tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_single_point_inputs_include_orca6_manual_fixtures`<br>`tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_scf_stability_inputs_include_orca6_manual_fixtures`<br>`tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_optimization_inputs_include_orca6_manual_fixtures`<br>`tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_frequency_inputs_include_orca6_manual_fixtures`<br>`tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_excited_state_inputs_include_orca6_manual_fixtures`<br>`tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_mrci_inputs_include_orca6_manual_fixtures`<br>`tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_solvator_inputs_include_orca6_manual_fixtures` |
| <!-- feature-area:Excited-state and multireference requests -->Excited-state and multireference requests | 部分支持 | 结构化激发态和多参考任务语义，包括 MRCI multi-job fixture。 | 覆盖示例之外的 ORCA 激发态和多参考关键字族可能仍为 raw。 | `tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_excited_state_manual_fixtures_are_structured`<br>`tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_mrci_manual_fixtures_are_structured`<br>`tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_mrci_multijob_manual_fixture_splits_frames_and_structures_last_job` |
| <!-- feature-area:Optimization, coordinates, frequency, and Solvator structures -->Optimization, coordinates, frequency, and Solvator structures | 部分支持 | 优化约束、内坐标、fragments、NEB、频率 restart 和 Solvator 示例。 | 覆盖范围来自显式 ORCA 手册 fixtures 和针对性回归示例。 | `tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_orca6_optimization_constraints_block_is_not_truncated`<br>`tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_orca6_manual_internal_coords_fill_frame_atoms_and_coords`<br>`tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_orca6_optimization_fragment_mixed_basis_is_structured`<br>`tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_orca6_optimization_neb_block_and_xtb_geometry_are_structured`<br>`tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_orca6_frequency_restart_example_is_structured`<br>`tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_orca6_solvator_fixture_is_structured` |
| <!-- feature-area:Writer availability -->Writer availability | 部分支持 | 从带坐标的模型规范化渲染 file/frame，支持显式 keywords、PAL、maxcore、任意 `%block`、笛卡尔坐标和 `.inp` 输出路径。 | 尚不规范化非笛卡尔内联几何和任意原文语句顺序；跨格式转换必须显式提供 ORCA keywords。 | `tests/test_autoparser_orcainp_tmpfile.py::test_orcainp_writer_builds_common_sp_input`<br>`tests/test_autoparser_orcainp_tmpfile.py::test_orcainp_writer_uses_inp_output_extension` |

## ORCA 6.1 能力覆盖矩阵

下表区分四种能力，避免把“能够保留 raw 文本”误认为“已经结构化支持”：

- **支持**：对所列语法有明确模型和测试，可以读取或规范化写出。
- **部分**：只覆盖已知形式或 fixture，不能承诺覆盖 ORCA 6.1 的全部变体。
- **Raw**：可以通过原始 `%block` 文本传递，但 MolOP 不解释或校验其中语义。
- **不支持**：当前会拒绝、丢失信息，或者不能从模型重新生成。

| ORCA 输入能力 | 解析 | 结构化语义 | 规范化渲染 | 原文保真 | 当前边界 |
| ------------- | ---- | ---------- | ---------- | -------- | -------- |
| `!` 简单关键词行 | 支持 | 部分 | 支持 | 不支持 | 可保留关键词文本；method、functional 和 basis 仅识别内置词表，输出时重新排版。 |
| 顶层 `#` 注释 | 部分 | 不适用 | 支持 | 不支持 | 保存独立注释行；不保留全局位置，也不支持 `# comment #` 后继续解析同一行。 |
| `%maxcore`、`%pal` | 支持 | 支持 | 支持 | 部分 | 支持 `maxcore`、`nprocs` 结构化覆盖；输出采用规范格式。 |
| 普通非嵌套 `%block ... end` | 支持 | 部分 | Raw | 部分 | 未建模的 block 可保存在 `ORCABlock.raw_text`，但不保证其与其他语句的原始相对位置。 |
| 单行 `%` 指令，如 `%moinp`、`%base` | 部分 | 不支持 | Raw | 部分 | 通用 block 扫描可以捕获；除资源字段外通常没有专用模型或参数校验。 |
| 嵌套 block，如 `%scf/SOSCF`、`%basis/NewGTO` | 部分 | 部分 | 部分 | 不支持 | 嵌套识别目前依赖有限白名单；未识别的内层 `end` 可能截断外层 block。 |
| 重复 block 和输入优先级 | 部分 | 不支持 | 部分 | 不支持 | 解析保留 block 列表顺序；按名称覆盖 block 时会移除旧实例并把新实例追加到末尾。 |
| 未分类顶层语句 | 部分 | 不支持 | 部分 | 不支持 | 存入 `trailing_lines`，渲染时统一移动到几何之后。 |
| `* xyz`、`* cart`、`* cartesian` | 支持 | 支持 | 支持 | 不支持 | 坐标规范化为 `* xyz`；空白、有效数字和原始坐标类型不保留。 |
| `%coords` Cartesian 几何 | 支持 | 支持 | 支持 | 不支持 | 解析后规范化为 `* xyz`，不重新生成 `%coords`。 |
| 点电荷 `Q q x y z` | 支持 | 支持 | 支持 | 不支持 | 输出电荷使用固定六位小数，坐标使用 writer 的精度设置。 |
| ghost、dummy、fragment、frozen、原子级基组 | 部分 | 部分 | 部分 | 不支持 | 覆盖现有 fixture 形式；原始 token、引号和排列不保证保留。 |
| isotope `M=`、nuclear charge `Z=` | 支持 | 支持 | 不支持 | 不支持 | parser 写入原子模型，但 Cartesian writer 尚未输出这两个属性。 |
| `* int`、`* internal`、`* gzmt` | 部分 | 部分 | 不支持 | 不支持 | 可解析已覆盖的数值内坐标并投影 Cartesian 坐标；writer 会抛出 `NotImplementedError`。 |
| `* xyzfile`、`* gzmtfile`、`* pdbfile` | 支持 | 支持 | 支持 | 部分 | 只保存并重新输出路径引用，不读取或验证外部文件。 |
| `%paras` 和坐标参数表达式 | 部分 | 部分 | 部分 | 不支持 | 支持简单参数、范围和 `name +/- offset`；几何渲染时表达式会被解析后的数值替代。 |
| `$new_job` 多任务 | 支持 | 支持 | 支持 | 不支持 | 每个 job 映射为一个 frame；输出统一使用 `$new_job` 分隔符。 |
| `%Compound` 工作流 | 部分 | 不支持 | 不支持 | 不支持 | 无顶层几何的合法 Compound-only 输入会被格式检测拒绝，writer 也要求关键词和几何。 |
| 模型化学与常见任务 | 部分 | 部分 | 部分 | 不适用 | 覆盖常见 SP、OPT、TS/NEB、FREQ、IRC、gradient、DFT、部分波函数方法和色散关键词。 |
| 激发态、多参考、Solvator | 部分 | 部分 | Raw | 部分 | parser 提供结构化请求；writer 依赖已解析的 raw block 或调用方显式传入 block，不从请求模型完整反向生成。 |
| 完整 `%basis`、ECP、AuxJ/AuxJK/AuxC/CABS | 部分 | 不支持 | Raw | 部分 | 简单关键词和原子级覆盖可识别；完整 block 内容不做字段级建模与校验。 |
| QMMM、MD、GOAT、DOCKER、RESP、Raman 等 ORCA 6.1 功能 | 部分 | 不支持 | Raw | 部分 | 非嵌套 raw block 可以传递；不承诺嵌套语法、任务兼容性或字段级构建。 |
| ORCA 版本和兼容性校验 | 不支持 | 不支持 | 不支持 | 不适用 | 当前版本记为 `Any`，不检查关键词、方法、基组和 ORCA 版本之间的兼容性。 |

当前 writer 是 **canonical writer**，适合从坐标和显式计算设置构建新输入；它不是
**source-preserving writer**。需要读取现有复杂输入、只修改一个字段并保持其他文本不变时，
不应依赖当前 round-trip 行为。

参考：[ORCA 6.1 输入结构](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/input.html)、
[坐标](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/coordinates.html)、
[基组](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/basisset.html)和
[Compound](https://www.faccts.de/docs/orca/6.1/manual/contents/workflowsautomatization/compound.html)。
