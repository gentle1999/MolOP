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

旧的 `refresh_common_qm_containers()` 仅用于兼容旧调用方：它只在结构化字段为空时从兼容字段补全，不允许覆盖已经解析出的结构化语义。

## 输入请求容器

`QMModelChemistry` 描述模型化学级别：

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

## 解析器约束

新增解析逻辑应遵守以下流程：

1. 解析程序特有语法，得到程序特有语义模型。
2. 构造公共结构化容器。
3. 如仍有旧字段上解析出的资源或关键词，只用 `backfill_common_qm_containers_from_legacy()` 补空字段。
4. 调用 `project_common_qm_fields()` 投影兼容字段。
5. 测试同时断言结构化字段与兼容字段一致。

不要把新语义只塞进 `keywords`、`resources_raw` 或 `options`。`options` 只保存尚未抽象成公共字段的程序特有细节。
