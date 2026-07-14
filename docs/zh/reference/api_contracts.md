# API 契约

本页固定 `format_transform`、`to_summary_df`、`parallel_execute` 和 CLI chain
的公开契约。

## `AutoParser`

Python API 可接收单个 `str` 或 `os.PathLike`，也可接收由路径和通配符组成的一维
iterable。每个通配符分别展开，全部结果物化后按规范化绝对路径去重，并以绝对路径顺序
返回一个 `FileBatchModelDisk`；不保留输入分组与输入顺序。

未匹配的通配符不贡献文件，空 iterable 返回空 batch。不存在的普通路径会交给 batch
parser 告警并跳过。嵌套 iterable 或非路径成员会抛出带成员索引的 `TypeError`。所有解析
选项统一应用于展开后的每条路径。

```python
from pathlib import Path

batch = AutoParser(["inputs/a.log", "inputs/set_b/*.log", Path("c.log")])
```

## 计算数据库视图

`ChemFile` 和 `ChemFileFrame` 是权威的计算数据模型。MolOP 不维护平行的
`Parsed*` DTO 对象树，也不提供第二套计算解析入口。数据库录入直接使用相同对象的
序列化视图：

```python
from molop import AutoParser

chem_file = AutoParser(
    "calculation.out",
    capture_source_evidence=True,
    release_file_content=True,
)[0]
file_payload = chem_file.to_unitless_dump_with_unit_keys(exclude_none=True)
frame_payloads = (
    frame.to_unitless_dump_with_unit_keys(exclude_none=True)
    for frame in chem_file
)
```

文件 payload 与帧 payload 有意分离。`ChemFile` dump 只包含 artifact、segment 和其他
文件级 metadata，不内嵌 `frames` 或私有 `_frames_` 集合。调用方应增量消费
`frame_payloads`，自行附加文件 identity 与顺序上下文，再交给自己的 Pydantic 接收模型
逐帧校验。接收模型只选择所需字段，通常对格式专用扩展使用 `extra="ignore"`。传输层可以
临时包装一个文件 payload 和一页帧 payload，但该包装不是 MolOP DTO，也不要求构造一个
巨型 JSON 文档。

开启 source evidence 后，每个 frame payload 包含三个不同作用域的序号：`frame_id` 是
当前 `ChemFile` 中的追加序号，`segment_frame_index` 是原文 segment 内序号，
`file_frame_index` 是完整原文件中按 locator 顺序展平后的稳定序号。即使
`only_last_frame=True`，后两者仍保留原文序号；未开启 evidence 时不输出这些 source
序号。数据库持久化原文顺序时必须使用 `file_frame_index`，不得使用可能因筛选而重置的
`frame_id`。

每个正常 parser 生命周期生成的文件 payload 都包含固定的
`schema_version="molop-calculation-export-v1"`、规范 `source_format` 和文件级
`parser_provenance`。
provenance 记录具体 parser、MolOP、MolGR、RDKit 版本，以及本次解析实际使用的 parser、
MolOP 和 MolGR 配置快照；`effective_config_sha256` 是该快照的严格 canonical JSON
SHA-256。provenance 不复制到 frame payload。

`to_unitless_dump_with_unit_keys()` 会递归地把 quantity 转为 magnitude，并将规范化
单位写入对应 key，格式为 `field (unit)`。嵌套 MolOP 模型、数组、序列和 mapping
都由 `BaseDataClassWithUnit` 的通用实现处理。其结果只是 dump，不是另一套需要独立
校验、转换或维护生命周期的数据模型。

单位 label 是该公开序列化视图的一部分。MolOP 根据 Pint 规范单位名使用自身的确定性
语法生成 label，不依赖用于展示的 `str(unit)`：乘积项排序后用 `*` 连接，幂使用 `^`，
多项分母使用括号。例如稳定 label 为 `hartree`、`kilocalorie/mole`、
`bohr^2*unified_atomic_mass_unit` 和 `calorie/(kelvin*mole)`。

数组默认转为 JSON-safe list。数据库录入可选择保留数值 ndarray 的独立副本，以便写入
确定性 NPY 或其他二进制 sidecar，而无需构造巨型 JSON：

```python
frame_payload = frame.to_unitless_dump_with_unit_keys(
    exclude_none=True,
    array_mode="ndarray",
)
```

`ndarray` 模式下，非数值数组仍转为 list；返回的数值数组不会与模型内部数组共享内存。
调用方用引用替换这些数组或选择编码前，返回值有意不保证可 JSON 序列化。其他关键字参数
会在递归转换前原样传给 Pydantic `model_dump()`，因此 `include`、`exclude`、
`exclude_none`、alias 以及其他字段筛选规则保持原有语义。

export 路径会在 Pydantic 生成快照前物化惰性的公开 topology。因此首个 frame payload
已经包含 `bonds`、`formal_charges` 和 `formal_num_radicals`，未修改 frame 的连续两次
dump 相等。bond 端点和逐原子 topology 数组与导出的 `atoms`、`coords` 使用相同的原文
原子顺序。该准备动作仅由此 export 方法触发；普通 `model_dump()` 保持既有的惰性行为与
成本。

可信重建还会输出 map-free `topology_v3000_molblock`、
`source_to_topology_atom_permutation` 和重建配置 provenance。当前 MolGR 后端经元素与坐标
核验为保序时，permutation 明确记录 identity；发生失败或无法无歧义确认原子顺序时不猜测
映射，并以 `parse_presence["topology"]="parse_failed"` 和稳定 diagnostic code 报告。

Source identity、location、semantics 和 evidence 都是现有模型上的可选字段。默认关闭
evidence 捕获；数据库
录入或审计场景可通过 `capture_source_evidence=True` 显式开启。支持的 parser 会在正常
解析与 byte offset 计算中统一使用 `source_encoding` 指定的严格解码器，默认值为
`utf-8`。支持的 parser 会在正常
解析路径中、释放保留原文之前填充这些字段，不会在解析结束后再运行一套 scientific
extractor。与格式无关的协议位于 `molop.io.base_models.source`：`DecodedSource` 保留
精确解码偏移，`LocatedTextBlock` 和 `LocatedSourceSegment` 表达 parser 提供的边界，
`SourceSpan` 保存 byte/character/line 范围。artifact 字段属于 `BaseChemFile`，位置字段
属于 `BaseChemFileFrame`，坐标来源字段属于 `BaseCoordsFrame`，因此坐标、输入和计算文件
都使用同一套生命周期。

每个 `SourceSegmentEvidence` 还包含 portable `protocol` 和 `task_requests`。
`protocol` 是通用 `model_chemistry` 的纯 mapping 投影，`task_requests` 是通用
`QMTaskRequest` 的纯 mapping 列表，不暴露 Gaussian/ORCA 专用 semantic model。证据来自
segment metadata，因此即使某段没有可解析 frame，也能保留其请求协议。正式数据库导入
必须设置 `capture_source_evidence=True`；此时 `SourceSpan` 的 byte、character 和 line
半开区间以及对应 SHA-256 都属于录入核验字段。

文件分段只有一个事实源。每个专用 file parser 必须声明规范 `format_id`，并实现
`_quick_check_file_format()` 和 `_locate_segments()`；locator 直接返回原文上的
segment/frame 半开区间。有 frame 的 segment 中，这些 frame 区间必须覆盖全部非空白
字符；纯空白间隙和零帧 segment 仍然合法。文件级与段级 metadata 可分别通过
`_parse_artifact_metadata()` 和 `_parse_segment_metadata()` 提供，
两者均为可选钩子。不存在 `_split_file()` 或 `_parse_metadata()` fallback，也不得先重建
frame 文本再反向搜索 span。

MolOP 负责输出解析事实与可选证据。存储系统仍需独立核验 artifact bytes，自行选择数组
持久化编码，完成数据库 identity canonicalization、admission/QC，以及 Reaction 或
ReactionPath 创建。

## Frame Selector

MolOP 公开转换与汇总 API 统一使用一个 frame selector 概念：
`int | Sequence[int] | "all"`。

| 入口 | 参数 | 默认值 | 返回映射 | 错误策略 |
| --- | --- | --- | --- | --- |
| Python 转换 | `frame` | `-1` | 渲染前归一化为 frame index；负整数从末尾计数。 | 序列中非整数抛出 `TypeError`；除 `"all"` 外的字符串抛出 `ValueError`。 |
| Python 汇总 | `frame` | `-1` | 按文件归一化；`"all"` 展开为每个文件的全部 frame。 | selector 类型错误立即抛出；缺失 frame 由 `on_missing_frame` 控制。 |
| CLI | `--frame` | `-1` | 从 `"all"`、单个整数或逗号分隔整数解析，再传给 Python API。 | 非法文本在执行 chain 前抛出 CLI 用法错误。 |

`slice` 不属于公开 frame selector 契约。文件对象自身的 Python 序列切片行为保留。

## `format_transform`

### Python API

| 参数 | 适用入口 | 默认值 | 含义 | 错误策略 |
| --- | --- | --- | --- | --- |
| `format` | frame、file、batch | 必填 | 目标 writer format id，如 `xyz`、`sdf`、`gjf`、`cml`。 | 不支持的格式抛出 registry writer 错误；writer 校验错误向上传播。 |
| `file_path` | frame、file | `None` | 仅在 `write_to_disk=True` 时作为输出文件路径。 | 不写盘时忽略；目录路径触发断言；memory-only 对象写盘但无路径时抛出 `ValueError`。 |
| `output_dir` | batch | `None` | 仅在 `write_to_disk=True` 时作为输出目录。 | 不写盘时忽略；Python API 写盘时要求目录已存在，CLI 会自动创建。 |
| `frame` | file、batch | `-1` | frame selector：`int | Sequence[int] | "all"`。 | selector 类型错误在 writer dispatch 前抛出；越界行为由具体 writer 决定，除非 writer 显式校验。 |
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
| `frame` | `-1` | frame selector：`int | Sequence[int] | "all"`；仅 frame mode 使用。 | selector 类型错误在汇总前抛出。 |
| `n_jobs` | `1` | 汇总提取并行度。 | joblib 与 worker 异常向上传播。 |
| `brief` | `True` | 传给 `to_summary_series`；`False` 请求扩展字段。 | 字段级错误由 summary 实现向上传播。 |
| `flatten_columns` | `False` | 将三层 MultiIndex 列转为 `General.FrameID`、`Energy.total_energy.hartree` 这样的点分列名；空单位层会跳过。 | 非 MultiIndex 列无影响。 |
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
| `--frame` | `"-1"` | CLI frame selector，转发为 `frame`。 | 非法 selector 文本抛出 CLI 用法错误。 |
| `--embed / --no-embed` | `--embed` | 控制多 frame 合并。 | writer 错误向上传播。 |
| `--write / --no-write` | 自动 | `--write` 在无 `--output-dir` 时写到源目录；`--no-write` 强制只渲染。 | 写盘错误向上传播。 |
| `--n-jobs` | parse 级 `--n-jobs` | 覆盖 batch transform 并行度。 | joblib 错误向上传播。 |
| 动态 writer 选项 | 无 | 转发给所选 writer。 | 动态参数解析错误抛出 CLI 用法错误；writer 错误向上传播。 |

操作实际写盘时 stdout 静默；只渲染时打印返回的 mapping。

### `to-summary-df`

| 参数 | 默认值 | 返回/效果 | 错误策略 |
| --- | --- | --- | --- |
| `--mode` | `"frame"` | 汇总模式；终止操作。 | 非法选项由 Click 拒绝。 |
| `--frame` | `"-1"` | CLI frame selector，转发为 `frame`。 | 非法 selector 文本抛出 CLI 用法错误。 |
| `--n-jobs` | parse 级 `--n-jobs` | 覆盖汇总并行度。 | joblib 错误向上传播。 |
| `--brief / --full` | `--brief` | 控制紧凑或扩展 summary 字段。 | summary 错误向上传播。 |
| `--flatten-columns / --multi-index-columns` | `--flatten-columns` | 控制是否输出 CSV/JSON 友好的扁平列名。 | 非 MultiIndex 列无影响。 |
| `--on-missing-frame` | `"skip"` | 缺失 frame 策略。 | `"error"` 下选中 frame 缺失时抛出 `IndexError`。 |
| `--out` | `None` | 将 summary 写到文件，否则输出到 stdout。 | 自动创建父目录；写文件错误向上传播。 |
| `--format` | `"csv"` | `--out` 或 stdout 的序列化格式。 | 非法选项由 Click 拒绝。 |
