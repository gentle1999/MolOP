# 模型字段

本页从公开 API 角度说明常用计算化学数据在 MolOP 文件对象和帧对象中的位置。这里列出的字段
适合在分析脚本、notebook 和下游转换中依赖。

## Gaussian 输出契约

`AutoParser(..., parser_detection="g16log")` 返回一个 batch。batch 中每个元素是一个
Gaussian 输出文件对象，文件对象中包含多个 frame。解析器会把 Gaussian 输出归一化成模型
字段；不承诺逐字节保留原始 log。

### 文件级字段

文件级字段适合读取整次任务的元数据，以及从最后一帧 finalize 出来的文件级状态。

| 数据 | 字段 | 说明 |
| ---- | ---- | ---- |
| 软件信息 | `file.qm_software`, `file.qm_software_version` | 通常为 `"Gaussian"` 和 Gaussian revision 字符串。 |
| Route 文本 | `file.keywords` | 归一化后的 route line 文本。 |
| 结构化 route | `file.semantic_route`, `file.model_chemistry`, `file.task_requests` | 推荐用于读取 method family、basis、任务类型、溶剂化、色散等语义。 |
| 兼容 route 投影 | `file.method`, `file.basis_set`, `file.functional` | 为兼容旧代码保留；能推断时从结构化 route 数据投影得到。 |
| 电荷与多重度 | `file.charge`, `file.multiplicity` | 从已解析 frame finalize 得到。 |
| 运行时间 | `file.running_time` | 当 Gaussian 打印 timing 记录时累计为 `pint` quantity。 |
| 终止状态 | `file.status` | 从最后一帧 finalize 得到的文件级状态。 |

### 帧级字段

帧级字段适合读取结构、优化步骤和逐帧 QM 结果。

| 数据 | 字段 | 说明 |
| ---- | ---- | ---- |
| 帧标识 | `frame.frame_id` | 文件内从 0 开始的 frame index。 |
| 结构 | `frame.atoms`, `frame.atom_symbols`, `frame.coords` | 原子序数、元素符号和 input-orientation 坐标。 |
| 标准取向 | `frame.standard_coords`, `frame.standard_orientation_transformation_matrix` | Gaussian 打印 standard orientation 或可重建时存在。 |
| 能量 | `frame.energies` | `Energies` 对象，包含 `reference_energy`、`electronic_energy`、post-HF 能量和计算字段 `total_energy`。 |
| 热力学 | `frame.thermal_informations` | `ThermalInformations` 对象，包含 ZPVE、热校正、热力学能量、熵、热容、质量以及转动/振动 metadata。 |
| 振动 | `frame.vibrations` | 频率、约化质量、力常数、IR 强度、振动模式向量和虚频数量。 |
| 分子轨道 | `frame.molecular_orbitals` | 轨道能级、占据数、对称性和派生前线轨道量。 |
| 布居 | `frame.charge_spin_populations` | Mulliken、Lowdin、Hirshfeld、CM5、NPA 和自旋布居字段，取决于 Gaussian 输出。 |
| 响应性质 | `frame.polarizability` | dipole、polarizability、electronic spatial extent 和多极矩。 |
| 力与 Hessian | `frame.forces`, `frame.hessian` | 归一化单位后的 Cartesian 数组。 |
| 几何优化 | `frame.geometry_optimization_status` | Berny 收敛数值、原始阈值和派生优化结果。 |
| 状态 | `frame.status`, `frame.is_error`, `frame.is_normal`, `frame.is_TS`, `frame.is_optimized` | 用于常见过滤和判断的公开状态字段。`is_optimized` 接受已优化的极小值和过渡态，但会拒绝超过一个虚频的 frame。 |
| 运行时间 | `frame.running_time` | 该帧 timing 信息，存在时为 `pint` quantity。 |

### 能量选择

`frame.energies.total_energy` 是计算字段，不是输入字段。它按以下优先级选择最具体的可用能量：

`ccsd_energy -> mp5_energy -> mp4_energy -> mp3_energy -> mp2_energy -> electronic_energy -> reference_energy`

如果需要明确能量来源，应读取具体方法字段；如果只需要“当前最佳可用标量能量”，使用
`total_energy`。

### 优化状态

当解析输出没有提供更强的结果时，`frame.geometry_optimization_status.geometry_optimized`
会从已有收敛指标推导。每个指标按
`convergence_multiplier * threshold` 判定，默认倍率为 `2.0`，比较时使用指标绝对值。

对于带频率结果的 frame，`frame.is_optimized` 接受极小值的零个虚频，也接受过渡态的
恰好一个虚频；超过一个虚频会被视为未优化。

### Summary 表

`file.to_summary_df()` 与 `batch.to_summary_df()` 默认汇总最后一帧（`frame=-1`）。
使用 `frame="all"` 时返回每个 frame 一行。默认 `brief=True` 时，稳定列覆盖存储信息、
电荷/多重度、结构摘要、route 元数据、环境和状态：

| 列组 | 示例 |
| ---- | ---- |
| `DiskStorage` | `FilePath`, `FileFormat` |
| `General` | `Charge`, `Multiplicity`, `CanonicalSMILES`, `NumAtoms`, `FrameID` |
| `Calc Parameter` | `Software`, `Version`, `Method`, `BasisSet`, `Functional`, `Keywords` |
| `Environment` | `SolventModel`, `Solvent`, `Temperature`, `Pressure` |
| `Status` | `IsError`, `IsNormal`, `IsTS`, `IsOptimized` |

Summary DataFrame 的列使用三层 `MultiIndex`：`(group, field, unit)`。无单位字段的
第三层为空；带单位字段把归一化单位放在第三层，例如
`("Energy", "total_energy", "hartree")`。

使用 `brief=False` 会额外加入结果密集字段，例如 `Energy`、`Thermal`、
`GeometryOptimizationStatus` 和 `Vibration`。优化 summary 列包含原始指标数值和原始阈值，
而不是只输出经过倍率处理后的收敛布尔值。

Batch 入口可使用 `batch.to_summary_df(frame="all", flatten_columns=True)` 获取所有选中
frame，并将列名展开为 `General.FrameID`、`Status.IsError`、
`Energy.total_energy.hartree` 这样的点分形式。

### fakeG 渲染

`file.format_transform("fakeg")` 遵循通用转换默认值 `frameID=-1`，因此默认渲染最后一帧。
如需全文件 Gaussian-like 输出，使用 `file.format_transform("fakeg", frameID="all")` 或
`file.render_fakeg()`。

`fakeg` 输出是语义化、归一化文本，适合人工检查、兼容性测试和重新解析检查；它不是
Gaussian log 的逐字复现。

### 兼容字段与内部数据

公开 API 保留了一些扁平兼容字段：

| 兼容字段 | 推荐结构化字段 |
| -------- | -------------- |
| `method`, `basis_set`, `functional` | `model_chemistry` |
| `keywords` | `semantic_route` / `model_chemistry` / task request 容器 |
| `is_error`, `is_normal`, `is_TS`, `is_optimized` | `status`, `vibrations`, `geometry_optimization_status` |

用于重建 Gaussian-like 文本的派生数据属于内部实现，不会进入 `model_dump()`。用户代码应依赖
文件/帧字段，而不是渲染实现细节。

## 通用字段映射

| 数据 | 位置 | 类型 |
| ---- | ---- | ---- |
| 结构 | `frame.atoms`, `frame.coords` | `list[int]`, `NumpyQuantity` |
| 键 | `frame.bonds` | `list` |
| SMILES | `frame.to_SMILES()` | `str` |
| 总能量 | `frame.energies.total_energy` | `PlainQuantity | None` |
| 参考能量 | `frame.energies.reference_energy` | `PlainQuantity | None` |
| 电子能量 | `frame.energies.electronic_energy` | `PlainQuantity | None` |
| 频率 | `frame.vibrations.frequencies` | `NumpyQuantity` |
| 虚频数量 | `frame.vibrations.num_imaginary` | `int` |
| 轨道 | `frame.molecular_orbitals` | `MolecularOrbitals | None` |
| 电荷/自旋布居 | `frame.charge_spin_populations` | `ChargeSpinPopulations | None` |
| 优化状态 | `frame.geometry_optimization_status` | `GeometryOptimizationStatus | None` |
| 计算状态 | `file.status` 或 `frame.status` | `Status | None` |
| QM 软件 | `file.qm_software`, `frame.qm_software` | `str` |
