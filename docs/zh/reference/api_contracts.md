# API 契约

本页固定 `format_transform`、`to_summary_df`、`parallel_execute` 和 CLI chain
的公开契约。

## Frame Selector

MolOP 公开转换与汇总 API 统一使用一个 frame selector 概念：
`int | Sequence[int] | "all"`。

| 入口 | 参数 | 默认值 | 返回映射 | 错误策略 |
| --- | --- | --- | --- | --- |
| Python 转换 | `frameID` | `-1` | 渲染前归一化为 frame index；负整数从末尾计数。 | 序列中非整数抛出 `TypeError`；除 `"all"` 外的字符串抛出 `ValueError`。 |
| Python 汇总 | `frameIDs` | `-1` | 按文件归一化；`"all"` 展开为每个文件的全部 frame。 | selector 类型错误立即抛出；缺失 frame 由 `on_missing_frame` 控制。 |
| CLI | `--frame` | `-1` | 从 `"all"`、单个整数或逗号分隔整数解析，再传给 Python API。 | 非法文本在执行 chain 前抛出 CLI 用法错误。 |

`slice` 不属于公开 frame selector 契约。文件对象自身的 Python 序列切片行为保留。

## `format_transform`

### Python API

| 参数 | 适用入口 | 默认值 | 含义 | 错误策略 |
| --- | --- | --- | --- | --- |
| `format` | frame、file、batch | 必填 | 目标 writer format id，如 `xyz`、`sdf`、`gjf`、`cml`。 | 不支持的格式抛出 registry writer 错误；writer 校验错误向上传播。 |
| `file_path` | frame、file | `None` | 仅在 `write_to_disk=True` 时作为输出文件路径。 | 不写盘时忽略；目录路径触发断言；memory-only 对象写盘但无路径时抛出 `ValueError`。 |
| `output_dir` | batch | `None` | 仅在 `write_to_disk=True` 时作为输出目录。 | 不写盘时忽略；Python API 写盘时要求目录已存在，CLI 会自动创建。 |
| `frameID` | file、batch | `-1` | frame selector：`int | Sequence[int] | "all"`。 | selector 类型错误在 writer dispatch 前抛出；越界行为由具体 writer 决定，除非 writer 显式校验。 |
| `embed_in_one_file` | file、batch | `True` | 格式支持时，将选中的多个 frame 合并到一个渲染结果。 | writer 错误向上传播。 |
| `write_to_disk` | frame、file、batch | `False` | 将渲染结果写盘；无显式路径时，disk-backed 对象写到源文件所在目录。 | 写盘错误向上传播；memory-only 对象无 `file_path` 时抛出 `ValueError`。 |
| `n_jobs` | batch | `1` | batch 转换并行度。 | joblib 与 worker 异常通过 `parallel_execute` 向上传播。 |
| `graph_policy` | 通过 `**kwargs` 传入 frame、file、batch | writer 默认 | 传给 registry 的分子图策略。 | 非法策略或缺失图数据按 registry/writer 错误处理。 |
| `**kwargs` | frame、file、batch | 无 | writer 特定选项，例如 Gaussian input 选项。 | 未知或非法选项由所选 writer 处理。 |

返回值：

| 入口 | 返回值 |
| --- | --- |
| `FrameFormatTransformMixin.format_transform(...)` | 单帧渲染字符串 `str`。 |
| `FormatTransformMixin.format_transform(...)` | 合并输出时为 `str`，否则为 `list[str]`。 |
| `BatchFormatTransformMixin.format_transform(...)` | `dict[str, str | list[str]]`，从源路径映射到渲染结果。 |

## `to_summary_df`

| 参数 | 默认值 | 含义 | 错误策略 |
| --- | --- | --- | --- |
| `mode` | `"frame"` | `"frame"` 汇总选中 frame；`"file"` 汇总每个文件。 | 其他值抛出 `ValueError`。 |
| `frameIDs` | `-1` | frame selector：`int | Sequence[int] | "all"`；仅 frame mode 使用。 | selector 类型错误在汇总前抛出。 |
| `n_jobs` | `1` | 汇总提取并行度。 | joblib 与 worker 异常向上传播。 |
| `brief` | `True` | 传给 `to_summary_series`；`False` 请求扩展字段。 | 字段级错误由 summary 实现向上传播。 |
| `flatten_columns` | `False` | 将 MultiIndex 列转为 `General.FrameID` 这样的点分列名。 | 非 MultiIndex 列无影响。 |
| `on_missing_frame` | `"skip"` | `"skip"` 跳过缺失 frame；`"error"` 抛错。 | 非法策略抛出 `ValueError`；`"error"` 下缺失 frame 抛出 `IndexError`。 |
| `**kwargs` | 无 | 额外 summary-series 选项。 | 转发错误向上传播。 |

返回值：`pandas.DataFrame`。如果没有任何选中的 summary series，返回空 DataFrame。

## `parallel_execute`

| 参数 | 默认值 | 含义 | 错误策略 |
| --- | --- | --- | --- |
| `func` | 必填 | 以 `func(diskfile, *args, **kwargs)` 调用。 | `func` 异常向上传播。 |
| `desc` | `""` | 进度条描述。 | 进度条 backend 错误向上传播。 |
| `n_jobs` | `1` | 并行任务数。 | joblib 错误向上传播。 |
| `return_as` | `"list"` | joblib 返回模式：`"list"`、`"generator"` 或 `"generator_unordered"`。 | 非法值按 joblib 错误处理。 |
| `*args` | 无 | 传给 `func` 的额外位置参数。 | callable 错误向上传播。 |
| `_diskfiles_snapshot` | `None` | 内部预计算 batch 快照，用于保持过滤、分组等调用中的结果对齐。 | 面向内部调用；错误快照可能造成调用方层面的对齐错误。 |
| `return_results` | `None` | `True` 总是返回结果；`False` 执行完并返回 `None`；`None` 会把全为 `None` 的 list 结果折叠为 `None`。 | generator 模式不会自动检查，除非显式用 `False` 执行完。 |
| `**kwargs` | 无 | 传给 `func` 的关键字参数。 | callable 错误向上传播。 |

返回值：`Iterable[R] | None`。

## CLI Chain

CLI 入口为 `molop parse PATTERN [parse options] OPERATION ...`。返回
`FileBatchModelDisk` 的操作可以继续链接；终止操作必须位于最后，且会在解析文件前校验。

### Parse 选项

| 参数 | 默认值 | 返回/效果 | 错误策略 |
| --- | --- | --- | --- |
| `PATTERN` | 必填 | 构建初始 `FileBatchModelDisk`。 | 缺失输入或无法匹配输入时按 parser 错误处理。 |
| `--parser-detection` | `"auto"` | 选择 parser detection 模式。 | 未知 parser id 按 parser/registry 错误处理。 |
| `--n-jobs`, `-j` | `-1` | 作为未单独设置 `--n-jobs` 的操作默认并行度。 | joblib 错误向上传播。 |
| `--output-format` | `"text"` | 控制终止结果输出格式。 | 非法选项由 Click 拒绝。 |

### `format-transform`

| 参数 | 默认值 | 返回/效果 | 错误策略 |
| --- | --- | --- | --- |
| `--format` | 必填 | 目标 writer format id；终止操作。 | 缺失值或不支持 writer 时抛出 CLI/registry 错误。 |
| `--output-dir` | `None` | 为兼容 CLI 旧行为，提供后隐式写盘。 | 写盘时执行前创建目录。 |
| `--frame` | `"-1"` | CLI frame selector，转发为 `frameID`。 | 非法 selector 文本抛出 CLI 用法错误。 |
| `--embed / --no-embed` | `--embed` | 控制多 frame 合并。 | writer 错误向上传播。 |
| `--write / --no-write` | 自动 | `--write` 在无 `--output-dir` 时写到源目录；`--no-write` 强制只渲染。 | 写盘错误向上传播。 |
| `--n-jobs` | parse 级 `--n-jobs` | 覆盖 batch transform 并行度。 | joblib 错误向上传播。 |
| 动态 writer 选项 | 无 | 转发给所选 writer。 | 动态参数解析错误抛出 CLI 用法错误；writer 错误向上传播。 |

操作实际写盘时 stdout 静默；只渲染时打印返回的 mapping。

### `to-summary-df`

| 参数 | 默认值 | 返回/效果 | 错误策略 |
| --- | --- | --- | --- |
| `--mode` | `"frame"` | 汇总模式；终止操作。 | 非法选项由 Click 拒绝。 |
| `--frame` | `"-1"` | CLI frame selector，转发为 `frameIDs`。 | 非法 selector 文本抛出 CLI 用法错误。 |
| `--n-jobs` | parse 级 `--n-jobs` | 覆盖汇总并行度。 | joblib 错误向上传播。 |
| `--brief / --full` | `--brief` | 控制紧凑或扩展 summary 字段。 | summary 错误向上传播。 |
| `--flatten-columns / --multi-index-columns` | `--multi-index-columns` | 控制是否输出 CSV/JSON 友好的扁平列名。 | 非 MultiIndex 列无影响。 |
| `--on-missing-frame` | `"skip"` | 缺失 frame 策略。 | `"error"` 下选中 frame 缺失时抛出 `IndexError`。 |
| `--out` | `None` | 将 summary 写到文件，否则输出到 stdout。 | 自动创建父目录；写文件错误向上传播。 |
| `--format` | `"csv"` | `--out` 或 stdout 的序列化格式。 | 非法选项由 Click 拒绝。 |
