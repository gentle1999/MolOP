# 格式转换与信息保留

`format_transform` 方法可在 `FileBatchModelDisk` 和单个文件对象上使用，用于在不同化学文件格式之间转换，并尽量保留目标格式能够表达的结构信息。

参数默认值、返回值与错误策略见 [API 契约](api_contracts.md)。

## 核心行为

- **帧选择**：默认情况下，仅转换最后一帧 (`frameID=-1`)。可以指定 `frameID="all"` 或一组帧 ID 来转换更多帧。
- **合并输出**：如果 `embed_in_one_file=True`（默认），多个帧将合并到一个输出文件中（如果格式支持，如 SDF 或多帧 XYZ）。
- **文件输出**：Python API 需要显式传入 `write_to_disk=True` 才会写盘。写盘时，
  `file_path` 或 batch 的 `output_dir` 用于选择输出位置；如果没有提供路径，则写到
  源文件所在目录。`write_to_disk=False` 时，`file_path` 和 `output_dir` 会被忽略，
  转换只返回渲染字符串或字符串列表。
- **结构层级**：
  - **COORDS (坐标级)**：`xyz` 和 `gjf` 等格式主要保留原子坐标和元素信息。`orcainp` 当前只提供结构化 reader。
  - **GRAPH (图级)**：`sdf`、`smi` 和 `cml` 等格式保留成键信息（分子图）。如果源文件仅包含坐标（例如 `.log` 文件），MolOP 将自动尝试使用其内置算法重建分子图。
- **元数据保留**：
  - `gjf` writer 会保留结构化 Gaussian 指令和关键字。`orcainp` reader 会把 ORCA 关键词、块和几何段解析到 frame 字段，但当前不注册 ORCA writer。
  - 转换为简单的坐标格式（如 XYZ）时，计算属性（能量、频率）通常**不会**保留，尽管某些格式（如 SDF）可以将其作为属性存储。
