# CLI 命令参考

MolOP 的业务入口是：

```text
molop [全局选项] parse PATTERN [解析选项] 操作 [操作参数] ...
```

先看任务配方：[CLI 初体验](getting_started/cli.md)和
[CLI 常用任务](guides/cli-recipes.md)。本页用于查参数和链式规则。

以下命令以[共享 ORCA 水分子样例](../assets/examples/water_mp2.out)为输入。将文件下载到当前
目录后即可直接运行；使用自己的文件时替换 `water_mp2.out`。

## 全局选项

| 选项 | 作用 |
| --- | --- |
| `-v, --verbose` | 显示更详细日志 |
| `-q, --quiet` | 抑制进度和非结果输出，适合脚本 |
| `--max-jobs N` | 为解析和操作设置进程级 worker 上限 |
| `--version` | 显示版本 |
| `-h, --help` | 显示帮助 |

## `parse` 选项

| 选项 | 默认值 | 作用 |
| --- | --- | --- |
| `PATTERN` | 必填 | 单路径或 glob |
| `--input PATH_OR_GLOB` | 无 | 追加路径或 glob；可重复传入 |
| `--parser-detection` | `auto` | 自动检测或格式 ID，如 `g16log`, `orcaout`, `xtbout` |
| `-j, --n-jobs` | `-1` | 解析进程数 |
| `--output-format` | `text` | 最终终端输出使用 `text` 或 `json` |
| `--total-charge` | 原文值 | 覆盖所有输入的分子电荷 |
| `--total-multiplicity` | 原文值 | 覆盖所有输入的自旋多重度 |
| `--only-extract-structure` | 关闭 | 跳过能量、频率等非结构结果 |
| `--only-last-frame` | 关闭 | 每个文件只保留最后一帧 |
| `--capture-source-evidence` | 关闭 | 捕获格式支持的原文位置、摘要与 parser provenance |
| `--source-encoding` | `utf-8` | 严格文本解码和原文偏移所用编码 |
| `--release-file-content / --keep-file-content` | 释放 | 解析后释放或保留原文文本 |
| `--force-unit-transform / --no-force-unit-transform` | 配置值 | 覆盖本次解析的单位转换策略 |
| `--graph-reconstruction-backend` | 配置值 | 本次解析使用 `cpp` 或 `python` 图重建后端 |
| `--make-dative-bonds / --no-make-dative-bonds` | 配置值 | 覆盖本次解析的配位键重建策略 |
| `--report` | 关闭 | 不返回 batch，改为保留每个输入的解析结果 |

`--n-jobs -1` 是默认的自动设置。MolOP 会使用当前进程可用 CPU 和可选的
`molopconfig.max_jobs` 上限；排查 parser 或原生库问题时传入 `--n-jobs 1`。parse 和每个操作
命令都独立默认使用 `-1`。使用全局 `--max-jobs N` 限制整条命令；某个操作需要更低上限时，在
该操作上显式传入正整数 `--n-jobs`。

电荷和多重度覆盖会作用于每个匹配输入；混合电荷或多重度的数据应拆成不同命令。
`--only-extract-structure` 是快速路径，会按约定省略非结构科学结果。单位和拓扑参数只覆盖本次
命令；未指定时继承 `molopconfig`，不会修改进程级全局配置。

可组合多个输入范围，无需依赖一个过宽的 glob：

```bash
molop -q parse "gaussian/*.log" \
  --input "orca/*.out" \
  --input selected.xyz \
  --only-last-frame \
  to-summary-df --full --out summary.csv
```

??? example "生成文件"

    ```text
    Summary written to summary.csv
    ```

## 逐文件解析报告

使用 `--report` 审计成功、缺失、不支持、空结果和解析失败的输入。该选项是终止解析模式，不能
再连接操作命令。

```bash
molop -q parse "results/*.log" \
  --input required.out \
  --report --output-format json > parse-report.json
```

??? example "生成文件"

    ```text
    parse-report.json
    ```

JSON 对象包含 `summary` 计数和保持输入顺序的 `outcomes` 数组。每条结果包含输入序号、绝对
路径、状态、检测格式、parser warning，以及可序列化的失败信息。默认文本输出是制表符分隔的
诊断表。

## 可继续链接的操作

| 操作 | 必要参数 | 结果状态 |
| --- | --- | --- |
| `filter-state` | `--state ts|error|opt|normal|thermal|no-img` | 新 batch |
| `filter-value` | `--target charge|multiplicity|format --value VALUE` | 新 batch |
| `filter-by-codec` | `--codec-id FORMAT_ID` | 新 batch |
| `sample` | `--n N [--seed SEED]` | 新 batch |

这些操作之后可以继续筛选，也可以接一个最终操作。

## 最终操作

| 操作 | 主要参数 | 输出 |
| --- | --- | --- |
| `to-summary-df` | `--frame`, `--full`, `--out`, `--format` | CSV/JSON 或终端表 |
| `format-transform` | `--format`, `--output-dir`, `--frame`, `--write/--no-write` | 渲染文本或文件 |
| `draw-grid-image` | `--out` | PNG/SVG 图像 |
| `groupby` | 分组参数 | 分组路径 |
| `copy-to` | 目标目录 | 复制源文件 |
| `move-to` | 目标目录 | 移动源文件 |

最终操作必须位于操作链末尾。

## Summary 参数

```bash
molop parse "water_mp2.out" to-summary-df --help
```

??? example "帮助摘要"

    ```text
    Usage: molop parse PATTERN to-summary-df [OPTIONS]
    --mode [file|frame]       Summary mode.  [default: frame]
    --frame TEXT              Frame selection: all, int, or csv ints.
    --brief / --full          Use compact or full summary fields.
    --out FILE                Output file.
    --format [csv|json]       Output file format.  [default: csv]
    ```

常用选项：

- `--mode frame|file`：frame 或 file 级摘要。
- `--frame -1|all|0,2`：最后、全部或指定 frame。
- `--brief/--full`：基础列或扩展科学结果。
- `--flatten-columns/--multi-index-columns`：扁平列或三层列。
- `--on-missing-frame skip|error`：缺失 frame 的处理。
- `--out PATH --format csv|json`：写盘位置与格式。

执行：

```bash
molop -q parse "water_mp2.out" --parser-detection orcaout --n-jobs 1 \
  to-summary-df --full --format json --out summary.json
```

输出对象包含以下字段（完整 JSON 还会包含其他已解析字段）：

??? example "JSON 输出"

    ```json
    {
      "Calc Parameter.Software": "ORCA",
      "Status.IsNormal": true,
      "Energy.total_energy.hartree": -74.9993745981
    }
    ```

    生成文件：`summary.json`。

## Transform 参数

```bash
molop parse "input.out" format-transform --help
```

??? example "帮助摘要"

    ```text
    Usage: molop parse PATTERN format-transform [OPTIONS] [EXTRA_ARGS]...
    --format TEXT           Target writer format id.  [required]
    --output-dir DIRECTORY  Directory for generated files.
    --frame TEXT            Frame selection: all, int, or csv ints.
    --write / --no-write    Write generated files.
    ```

- `--format FORMAT_ID` 必填。
- `--frame -1|all|0,2` 选择 frame。
- `--embed/--no-embed` 控制多 frame 合并。
- `--output-dir DIR` 指定目录并默认写盘。
- `--write` 无输出目录时写到源目录。
- `--no-write` 即使给出输出目录也只输出渲染结果。

writer 特定参数跟在静态参数之后，例如：

```bash
molop parse "input.out" \
  format-transform --format orcainp --output-dir next \
  --keywords "B3LYP def2-SVP Opt" --nprocs 8 --maxcore 2000
```

??? example "生成文件"

    ```text
    next/input.inp
    ```

## 查看精确 help

```bash
molop --help
molop parse --help
molop parse PATTERN filter-state --help
molop parse PATTERN to-summary-df --help
molop parse PATTERN format-transform --help
```

??? example "帮助命令的共同形状"

    ```text
    molop [OPTIONS] parse PATTERN OPERATION [OPTIONS]
    ```

内部 plan 校验、终止操作约束和动态补全协议见
[CLI 进阶契约](advanced/cli-contract.md)。
