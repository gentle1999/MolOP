# QM 公共数据模型设计

本文定义 MolOP 中 QM 输入与结果数据容器的职责边界。目标是让 Gaussian、ORCA 和后续程序的解析结果能进入同一套公共模型，同时避免同一语义在多个字段之间互相覆盖。

## 分层原则

公共 QM 数据分三层：

1. 兼容字段：`method`、`functional`、`basis_set`、`keywords`、`request_num_cpu`、`request_memory`、`resources_raw`。这些字段保持现有 API 可用，但不是新解析逻辑的权威来源。
2. 输入请求容器：`model_chemistry`、`task_requests`、`resource_request`、`excited_state_requests`、`multireference_requests`。解析器应优先写入这些结构化字段。
3. 结果容器：`energies`、`molecular_orbitals`、`vibrations`、`single_point_properties`、`electronic_states`、`multireference_result`。这些字段描述计算输出，不反向推断输入请求。

## 权威字段

结构化容器是权威字段。解析器完成语义识别后，应调用 `project_common_qm_fields()`，把结构化容器投影到兼容字段：

- `model_chemistry.raw_keywords` -> `keywords`
- `model_chemistry.method_family/method` -> `method`
- `model_chemistry.functional` -> `functional`
- `model_chemistry.basis_set` -> `basis_set`
- `resource_request.num_cpu/memory/raw` -> `request_num_cpu/request_memory/resources_raw`

旧的 `refresh_common_qm_containers()` 仅用于兼容旧调用方：它先在结构化字段为空时从兼容字段补全，
再把结构化字段投影回兼容字段；已经解析出的结构化语义不会被平铺字段覆盖。需要修改规范化语义时，
应直接修改结构化容器，再调用 `project_common_qm_fields()`。

## 输入请求容器

`QMModelChemistry` 描述模型化学级别：

- `method` 是用于匹配和计算协议 identity 的稳定、具体电子结构方法名。Gaussian
  的 `R`、`U`、`RO` 前缀不写入该值，而是单独保存在 `spin_treatment`；原始 route
  token 仍保留在 Gaussian semantic route evidence 中。
- 泛函与色散属于同一模型化学表达。Gaussian 和 ORCA 都应把色散后缀写进 `functional`，例如 `B3LYP-GD3BJ`。
- `dispersion_correction` 同时保留规范化的色散标识，例如 `GD3BJ` 或 `D4`。
- `basis_sets` 用于混合基组、辅助基组、ECP 和原子级覆盖；`basis_set` 只表示主全局轨道基组。
- `solvation_model` 和 `solvent` 表示输入请求，不等同于输出层的 `BaseCalcFrame.solvent`。

`QMTaskRequest` 只描述通用任务类型与通用选项，例如 `sp`、`opt`、`freq`、`excited_state`、`multi_reference`。激发态和多参考的专有参数放在专门容器中。

`ExcitedStateRequest` 描述输入侧的根、态数、自旋通道、root following 和请求性质。

`MultireferenceRequest` 描述输入侧的多参考方法、参考方法、active space、state blocks、阈值和修正。

## 结果容器

`ElectronicStates` 表示态/根分辨的输出结果，适合激发态、MRCI、CASSCF 等结果。它不替代 `energies`，后者仍用于总体或基态能量。

`MultireferenceResult` 表示多参考输出结果，包括多参考方法、active space、态分辨结果、修正和诊断。它不替代输入侧的 `MultireferenceRequest`。

## 跨格式结果语义

以下约定统一用于 Gaussian 日志、FCHK、ORCA 日志和 xTB 标准输出。
回归测试位于 `tests/test_qm_result_semantics.py` 及各格式的数值测试中。

- **状态身份**：`state_index` 是内部唯一的非负编号，`root` 保留程序根编号。
  `is_ground_state=True` 明确表示基态，`False` 表示激发态，`None` 表示未知。
  不按列表位置或根编号零猜测基态。只有激发态的集合仍导出全部已知跃迁；
  未存储基态记录时，跃迁的内部初态编号为 `None`，打印的初态标签及根编号保留在跃迁属性中。
- **缺失观测值**：未知轨道能量、键级及耦合矩阵项在内部使用数值型 `np.nan`；
  未知占据数和对称性标签使用 `None`。仅明确允许缺失的字段将 NaN 导出为 JSON `null`，
  二进制 ndarray 导出保留 NaN；无穷值仍不能导出为 JSON。实际打印的零保持零。
  `BondOrders` 的每条 `BondOrderMatrix` 保存打印阈值和缺失原因，RDKit 嵌入跳过未观测原子对。
- **自旋通道**：公共轨道校验器补齐较短或缺席的通道，但不据此推断限制性。
  解析器须明确识别共享的限制性轨道；`orbital_basis`、`occupancy_origin` 和
  `total_occupancies` 区分打印的自旋占据数与从共享空间轨道推导的占据数。
  闭壳层总占据数等分；有依据的整数高自旋总占据数可拆分到两侧；其他开壳层情况保留未知自旋占据数。
  xTB 未打印的轨道区间不视为未占据。打印的 HOMO/LUMO 标记写入
  `reported_homo_index`、`reported_lumo_index`；仅凭分数占据数不推断前线轨道或 SOMO。
- **振动模式**：解析后的 `mode_indices` 从零开始，作为稳定内部编号；
  `source_mode_indices` 单独保存打印编号，只有数组而没有打印编号时为 `None`。
  筛选保留两种编号。单个 `Vibration` 及 IR/Raman 谱带继续携带来源程序、方向、归一化和质量加权信息。
  不通过缩放源向量来强行统一不同程序的归一化结果。
- **热化学**：修正量与总量保持既有定义，例如 `U_T = Eel + TCE`、
  `H_T = Eel + TCH`、`G_T = Eel + TCG`。容器在来源可确认时同时保存温度、压力、
  `entropy_model`、频率缩放、排除模式的截断、谐振子/转子插值参考频率及虚频截断，
  property bundle 导出继续保留这些条件。转子插值参考频率不等于排除低频模式的阈值；
  仅有 FCHK 热化学总量时，不猜测熵模型。

Gaussian 拉曼活性与 P/U 退偏振比独立存储：
`raman_depolarization_ratios_polarized` 为平面偏振，`raman_depolarization_ratios_unpolarized` 为非偏振；
无后缀字段仅用于来源未明确偏振约定的比值。xTB PTB 频率和活性分别保留
`frequency_method`、`raman_activity_method`，FCHK 热校正的推导式保存在 `derivations`。
公共容器概览和用户侧比较规则见[数据容器与跨 QM 语义](../../reference/data_model.md)。

迁移现有调用时，手动构造电子态需要显式设置 `is_ground_state`；需要原始振动编号时，
改用 `source_mode_indices` 或单模式的 `Vibration.source_mode_index`。跨程序比较时，不将未知矩阵项替换为零。

## 可扩展键级记录

键级与 population 使用相同的容器设计：`BondOrders.bond_orders` 是稳定标识符到
`BondOrderMatrix` 的字典。`scheme` 是可扩展字符串；自旋通道、来源标签、打印阈值、
缺失原因和程序元数据保存在对应记录中，不另设按名称关联的平行字典。

```python
import numpy as np
from molop.io.base_models.DataClasses import BondOrderMatrix, BondOrders

orders = BondOrders(bond_orders={
    "mayer_bond_order": BondOrderMatrix(
        scheme="mayer",
        values=np.array([[np.nan, 0.9], [0.9, np.nan]]),
        source_label="Mayer bond orders larger than 0.1",
        print_threshold=0.1,
        missing_reason="未打印的原子对",
    ),
})
matrix = orders["mayer_bond_order"].values
```

新增键级方案直接增加字典记录，不修改 serializer 或导出字段清单。
各矩阵须为同尺寸、非空方阵；帧有结构时，尺寸还须与源原子数一致。
`bond_order_items()`、`bond_order_names` 和 `get_bond_order()` 提供通用访问；
property bundle 和 RDKit 嵌入均遍历这些记录。RDKit 跳过未知项，对已有键保留已观测的零值。

旧调用 `BondOrders(mayer_bond_order=matrix)` 迁移为上述字典构造；
`frame.bond_orders.mayer_bond_order` 改为 `frame.bond_orders["mayer_bond_order"].values`。
旧平铺字段不再接受，避免记录被忽略。序列化的矩阵位于
`bond_orders[记录名]["values"]`，其 NaN 仍导出为 JSON `null`。

## 解析器约束

新增解析逻辑应遵守以下流程：

1. file parser 声明规范 `format_id`，实现 `_quick_check_file_format()` 和 `_locate_segments()`，直接定位原文区间。
2. 仅在需要时实现 `_parse_artifact_metadata()` 或 `_parse_segment_metadata()`。
3. 解析程序特有语法，得到程序特有语义模型并构造公共结构化容器。
4. 如仍有旧字段上解析出的资源或关键词，只用 `backfill_common_qm_containers_from_legacy()` 补空字段。
5. 调用 `project_common_qm_fields()` 投影兼容字段，并同时测试结构化字段与投影一致。

不得重新引入 `_split_file()` 或 `_parse_metadata()`。locator 是 segment/frame 边界和送入
frame parser 的原文内容的唯一事实源。

不要把新语义只塞进 `keywords`、`resources_raw` 或 `options`。`options` 只保存尚未抽象成公共字段的程序特有细节。
