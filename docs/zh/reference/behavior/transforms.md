# 格式转换与信息保留

`format_transform` 方法可在 `FileBatchModelDisk` 和单个文件对象上使用，用于在不同化学文件格式之间转换，并尽量保留目标格式能够表达的结构信息。

参数默认值、返回值与错误策略见 [API 契约](../../developer/contracts/api.md)。

=== "内存渲染"

    ```python
    from molop import AutoParser

    batch = AutoParser("results/*.out", n_jobs=1)
    rendered = batch.format_transform("xyz", frame=-1, write_to_disk=False)
    ```

    ??? example "返回形状"

        ```text
        dict 1
        ```

    返回值是从每个绝对源路径映射到渲染文本的 mapping；当
    `embed_in_one_file=False` 时，值可以是字符串列表。此模式不会创建输出文件。

=== "写入磁盘"

    ```python
    from pathlib import Path
    from molop import AutoParser

    Path("structures").mkdir(parents=True, exist_ok=True)
    batch = AutoParser("results/*.out", n_jobs=1)
    batch.format_transform(
        "xyz",
        output_dir="structures",
        write_to_disk=True,
    )
    ```

    ??? example "生成文件"

        ```text
        structures/water_mp2.xyz
        ```

    批量 Python API 要求 `output_dir` 已存在；CLI 会创建 `--output-dir`。两种 API 都只替换
    源文件的最后一个后缀。

!!! warning
    目标格式无法表达的属性不会被转换保留。需要保留成键数据时使用支持分子图的格式；需要
    保留能量或频率时使用支持 QM 数据的目标格式。

## 核心行为

- **帧选择**：默认情况下，仅转换最后一帧 (`frame=-1`)。可以指定 `frame="all"` 或一组帧 ID 来转换更多帧。
- **合并输出**：如果 `embed_in_one_file=True`（默认），多个帧将合并到一个输出文件中（如果格式支持，如 SDF 或多帧 XYZ）。
- **文件输出**：Python API 需要显式传入 `write_to_disk=True` 才会写盘。写盘时，
  `file_path` 或 batch 的 `output_dir` 用于选择输出位置；如果没有提供路径，则写到
  源文件所在目录。输出文件名仅替换最后一个后缀，例如转换为 XYZ 时，
  `name.hash.log` 会变为 `name.hash.xyz`。`write_to_disk=False` 时，`file_path` 和
  `output_dir` 会被忽略，
  转换只返回渲染字符串或字符串列表。
- **结构层级**：
    - **COORDS (坐标级)**：`xyz`、`gjf` 和 `orcainp` 等格式主要保留原子坐标和元素信息。
    - **GRAPH (图级)**：`sdf`、`smi` 和 `cml` 等格式保留成键信息（分子图）。如果源文件仅包含坐标（例如 `.log` 文件），MolOP 将自动尝试使用其内置算法重建分子图。
- **元数据保留**：
    - `gjf` writer 会保留结构化 Gaussian 指令和关键字。`orcainp` reader 会把 ORCA 关键词、块和几何段解析到 frame 字段；已注册的 canonical writer 可用显式 `keywords`、资源和 block 从带坐标 frame 构建新输入，但不承诺原文无损 round-trip。
    - 转换为简单的坐标格式（如 XYZ）时，计算属性（能量、频率）通常**不会**保留，尽管某些格式（如 SDF）可以将其作为属性存储。
