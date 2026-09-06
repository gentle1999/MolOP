# API 契约

本页按调用路径整理 batch、frame 选择、格式转换、汇总和 CLI 的公开契约。Parser 的实现
规则见[完整 Parser 契约](../developer/parser-contract.md)。

## 快速索引

| 需求 | API | 结果 |
| --- | --- | --- |
| 解析单个文件 | [`AutoFileParser(...)`](#autofileparser) | file 级对象 |
| 解析已加载的文本或 bytes | [`AutoMemoryParser(...)`](#automemoryparser) | 内存 file 级对象 |
| 解析文件和 glob | [`AutoParser(...)`](#autoparser) | `FileBatchModelDisk` |
| 选择 frame | [Frame 选择](#frame-selector) | 规范化后的 frame 索引 |
| 格式转换 | [`format_transform(...)`](#format_transform) | 渲染文本或路径映射 |
| 绘制轨迹 | [`draw_animation(...)`](#draw_animation) | GIF 或动画 SVG |
| 导出 TS 前后体 | [`save_pre_post_ts(...)`](#save_pre_post_ts) | 前后体路径映射 |
| 生成表格 | [`to_summary_df(...)`](#to_summary_df) | `pandas.DataFrame` |
| 对每个文件运行函数 | [`parallel_execute(...)`](#parallel_execute) | 结果或 `None` |
| 使用命令行 | [CLI 命令参考](../command_line_interface.md) | 末端操作结果 |

## `AutoFileParser`

```python
from molop import AutoFileParser

parsed_file = AutoFileParser("water_mp2.out")
```

`AutoFileParser` 只接受一个已存在的文件路径，会探测 reader 候选并直接返回 file 级对象。
它不会展开 glob、构造 batch parser 或调度 worker 进程，适合由调用方自行控制进程或任务模型。
默认 `parser_detection="auto"`，也可以显式传入格式 ID。文件不存在、不支持的格式、格式不匹配
和解析失败都会抛出异常。

## `AutoMemoryParser`

```python
from molop import AutoBytesParser, AutoTextParser

parsed_text = AutoTextParser("1\nwater\nH 0.0 0.0 0.0\n", parser_detection="xyz")
parsed_bytes = AutoBytesParser(raw_bytes, parser_detection="xyz")
```

`AutoMemoryParser` 接受调用方已经持有的文本、bytes-like 数据或文本/二进制流，返回内存中的
file 级对象。`AutoTextParser` 和 `AutoBytesParser` 用于显式区分输入类型。源数据没有文件扩展名
时应使用 `parser_detection`；默认的 `"auto"` 会尝试已注册的内存 reader 及其内容探测。二进制
解析使用 `source_encoding` 解码；设置 `capture_source_evidence=True` 时仍保留精确源字节信息。

## `AutoParser`

```python
from molop import AutoParser

batch = AutoParser("water_mp2.out", n_jobs=1)
```

输入可以是单个路径、`os.PathLike`，或由路径和 glob 组成的一维 iterable。MolOP 负责展开
glob，将匹配项转为绝对路径、去重，并按规范化路径排序；不保留输入分组和输入顺序。

未匹配的 glob 不产生文件；空 iterable 返回空 batch。不存在的普通路径交给 batch parser，
随后告警并跳过。嵌套 iterable 或非路径成员会抛出带输入索引的 `TypeError`。解析选项应用于
每个展开后的路径。

## Frame 选择 {#frame-selector}

Python 转换和汇总 API 使用同一种选择器：`int | Sequence[int] | "all"`。

| 值 | 含义 |
| --- | --- |
| `-1` | 最后一个 frame；负整数从末尾计数。 |
| `0`、`[0, 2]` | 单个 frame 或显式 frame 列表。 |
| `"all"` | 全部 frame。 |

除 `"all"` 外的字符串抛出 `ValueError`；序列中出现非整数抛出 `TypeError`。汇总时，
`on_missing_frame="skip"` 丢弃越界索引，`"error"` 抛出 `IndexError`。`slice` 不属于本
API 契约。

## `format_transform`

按目标输出选择调用层级：

=== "Python"

    ```python
    rendered = batch.format_transform(
        "xyz",
        frame=-1,
        embed_in_one_file=True,
        write_to_disk=False,
        n_jobs=1,
    )
    print(type(rendered).__name__, len(rendered))
    ```

=== "CLI"

    ```bash
    molop -q parse "water_mp2.out" \
      format-transform --format xyz --frame -1 --no-write
    ```

??? example "返回形状"

    ```text
    dict 1
    ```

    对单个输入文件，Python 调用返回一个包含源路径和渲染文本的 `dict`；CLI 在终端输出渲染
    文本，不创建 `.xyz` 文件。需要处理 glob 时，把示例路径替换为 `results/*.out`。

| 参数 | 适用层级 | 默认值 | 契约 |
| --- | --- | --- | --- |
| `format` | frame、file、batch | 必填 | 已注册的 writer ID，例如 `xyz`、`sdf`、`gjf`、`cml`。 |
| `frame` | file、batch | `-1` | `int`、整数序列或 `"all"`。 |
| `embed_in_one_file` | file、batch | `True` | writer 支持时，将所选 frame 合并到一个输出。 |
| `write_to_disk` | frame、file、batch | `False` | 写入文件；未提供路径时，磁盘对象使用由源文件推导的路径。 |
| `file_path` / `output_dir` | frame/file / batch | `None` | 只在写盘时生效。内存对象写盘必须提供 `file_path`；Python batch 的 `output_dir` 必须已存在。 |
| `n_jobs` | batch | `1` | batch 转换使用的 worker 数。 |
| `**kwargs` | frame、file、batch | - | writer 专用选项，例如 Gaussian route 或 link-0 参数。 |

返回值分别为：单 frame 返回 `str`，file 返回 `str` 或 `list[str]`，batch 返回
`dict[str, str | list[str]]`。不支持的格式、writer 校验错误和写盘错误会传递给调用方。
Batch 转换对单个文件失败时记录 warning，并为该文件返回空字符串或空列表；出现空值时应
检查日志。

writer domain、frame 组合和格式专用行为见[转换行为](transform_behavior.md)。

## `draw_animation`

```python
animation = parsed_file.draw_animation(
    image_format="gif",
    file_path="trajectory.gif",
    duration=120,
)
```

`draw_animation` 会渲染所有能够得到 RDKit 分子图的 frame；无法绘制的 frame 会被跳过。默认
legend 保留 frame ID，并在可用时加入 TS 标记和总能量。`image_format` 支持 `"gif"` 与
`"svg"`；`duration` 可以是一个正整数，也可以按原始 frame 提供序列。发生跳帧时，
`duration`、`legends` 和 highlight 列表等按 frame 对齐的序列会同步过滤。

## `save_pre_post_ts`

```python
file_exports = parsed_file.save_pre_post_ts("ts-endpoints", format="sdf")
batch_exports = batch.save_pre_post_ts("ts-endpoints-batch", format="sdf", n_jobs=2)
```

文件级返回 `{frame_id: (pre_path, post_path)}`；batch 级返回值再增加源文件路径这一层。混合
batch 中不支持该操作的文件会记录 warning 并跳过。唯一 source stem 保持原名；重复 stem 会追加
稳定的源路径摘要，避免覆盖其他文件的结果。导出的前后体是基于几何的候选，不是优化后的真实
反应物或产物。前后体推断默认在每个振动方向上采样 `steps=8` 个振幅，范围为
`min_ratio=0.6` 至 `max_ratio=1.4`；每侧选择出现频次最高的拓扑，再将碎片数更多的候选
判为前体。`sampling_method` 默认为 `"harmonic_potential"`，对振幅平方等距采样，从而在
简谐近似下使势能（虚频模式取势能
降低量的绝对值）等距。该参数同样适用于 `additional_pre_post_ts`、`to_diff_rdmol` 和
`save_pre_post_ts`。

标准前后体结果不会被额外采样替代。保留标准前后体的三维构象后，可调用
`additional_pre_post_ts(pre, post)` 生成额外采样的一对表示；成键或断键相关原子固定在原前后体坐标，
其余原子沿虚频振动模式重新扫描。如果不存在键变化，则直接返回传入的前后体。

## `to_summary_df`

```python
summary = batch.to_summary_df(
    mode="frame",
    frame=-1,
    brief=False,
    flatten_columns=True,
    on_missing_frame="skip",
)
print(summary.shape)
```

??? example "输出"

    ```text
    (1, 23)
    ```

| 参数 | 默认值 | 契约 |
| --- | --- | --- |
| `mode` | `"frame"` | 汇总所选 frame；使用 `"file"` 时每个文件一行。 |
| `frame` | `-1` | frame 选择器；仅在 frame 模式使用。 |
| `brief` | `True` | `True` 使用紧凑字段；`False` 请求展开后的结果字段。 |
| `flatten_columns` | `False` | 将三级列索引转换为 `Energy.total_energy.hartree` 等名称。 |
| `on_missing_frame` | `"skip"` | 丢弃缺失 frame，或用 `"error"` 抛出 `IndexError`。 |
| `n_jobs` | `1` | batch 汇总使用的 worker 数。 |
| `**kwargs` | - | 转发给 `to_summary_series` 的其他选项。 |

返回 `pandas.DataFrame`；没有可选行时返回空 DataFrame。[批量汇总指南](../guides/batch.md)
和[批量汇总 Notebook](../examples/02-batch-summary-filter-select.ipynb) 展示共享样例的
完整表格及其全部生成列。

## `parallel_execute`

```python
def describe(parsed_file):
    return parsed_file.filename, parsed_file.detected_format_id

rows = batch.parallel_execute(
    describe,
    n_jobs=1,
    return_results=True,
)
print(list(rows))
```

??? example "输出形状"

    ```text
    [('water_mp2.out', 'orcaout')]
    ```

MolOP 对 batch 中每个文件调用一次 `func(diskfile, *args, **kwargs)`。`n_jobs` 控制执行并行
度，`desc` 设置进度说明，`return_as` 接受 Joblib 的 `"list"`、`"generator"` 和
`"generator_unordered"`。仅执行副作用操作时设置 `return_results=False`；否则返回 iterable。
回调函数或 Joblib 抛出的异常会传递给调用方。

`_diskfiles_snapshot` 是 MolOP 内部操作用于保持结果对齐的参数，普通应用代码不应使用。

## CLI 与序列化

- [CLI 命令参考](../command_line_interface.md)说明解析、筛选、汇总和转换的命令链。
- [Source evidence 与序列化](serialization.md)说明如何把文件级 metadata 与 frame 记录导出给
  数据库或审计系统。

## 错误与边界

| 情况 | 行为 |
| --- | --- |
| 输入或 parser 不支持 | 文件跳过并记录 parser diagnostics；显式非法 selector 或 option 直接抛错。 |
| Frame selector 非法 | 非法字符串抛 `ValueError`；序列成员非整数抛 `TypeError`。 |
| 汇总 frame 缺失 | 默认跳过；`on_missing_frame="error"` 抛 `IndexError`。 |
| writer 不支持或校验失败 | 对象级转换直接传递错误。 |
| Batch 转换单文件失败 | 记录 warning，并为该源返回空字符串或空列表。 |

## 相关页面

- [解析文件](../guides/parsing.md)
- [批量汇总](../guides/batch.md)
- [转换行为](transform_behavior.md)
- [Source evidence 与序列化](serialization.md)
- [格式支持](format_support.md)
- [完整 Parser 契约](../developer/parser-contract.md)
