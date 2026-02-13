# 命令行接口

MolOP 提供了一个基于 [Typer](https://typer.tiangolo.com/) 构建的现代命令行接口。主要入口是 `molop
<命令>`，它为与分子数据的交互提供了一种结构化且用户友好的方式。

为保证仓库内可复现，本文示例统一使用 `uv run molop ...`。如果你已在环境中安装 MolOP，可直接使用 `molop ...`。

## 全局选项

MolOP 支持几个可用于任何命令的全局选项：

- `-v`, `--verbose`: 启用详细输出以进行详细日志记录。
- `-q`, `--quiet`: 启用安静模式。强烈建议在脚本或文档中使用此模式以获得稳定、干净的输出。
- `--version`: 显示版本并退出。

## 发现命令

你可以通过运行以下命令查看所有可用的命令和全局选项：

```bash
uv run molop --help
```

```text

 Usage: molop [OPTIONS] COMMAND [ARGS]...

 MolOP: Molecule OPerator CLI

╭─ Options ────────────────────────────────────────────────────────────────────╮
│ --verbose             -v        Enable verbose output.                       │
│ --quiet               -q        Enable quiet mode.                           │
│ --version                       Show the version and exit.                   │
│ --install-completion            Install completion for the current shell.    │
│ --show-completion               Show completion for the current shell, to    │
│                                 copy it or customize the installation.       │
│ --help                -h        Show this message and exit.                  │
╰──────────────────────────────────────────────────────────────────────────────╯
╭─ Commands ───────────────────────────────────────────────────────────────────╮
│ chain             Legacy chain compatibility mode. Forwards all arguments to │
│                   the Fire-based CLI.                                        │
│ summary           Generate a summary of the molecules in the given files.    │
│ visualize         Visualize molecules in a grid image.                       │
│ transform         Transform molecular files to another format.               │
│ filter-by-codec   Filter files by their detected codec ID.                   │
│ sample            Randomly sample N files from the matched files.            │
│ stats             Show statistics of the molecules in the given files.       │
│ groupby           Group files by a key and output the groups.                │
╰──────────────────────────────────────────────────────────────────────────────╯
```

## Typer 命令

### `summary`

生成给定文件中分子的摘要，通常输出为 CSV 或 JSON 文件。

**何时使用：** 当你需要将一批文件中的关键属性（如能量、SMILES 等）提取到结构化表格中时，请使用此命令。

```bash
uv run molop summary --help
```

```text

 Usage: molop summary [OPTIONS] PATTERN

 Generate a summary of the molecules in the given files.

╭─ Arguments ──────────────────────────────────────────────────────────────────╮
│ *    pattern      TEXT  File pattern to match. [required]                    │
╰──────────────────────────────────────────────────────────────────────────────╯
╭─ Options ────────────────────────────────────────────────────────────────────╮
│ --out               -o      PATH     Output file path.                       │
│ --format            -f      TEXT     Output format (csv or json).            │
│                                      [default: csv]                          │
│ --mode              -m      TEXT     Summary mode (file or frame).           │
│                                      [default: frame]                        │
│ --frame                     TEXT     Frame selection (e.g., 'all', '-1',     │
│                                      '1,2').                                 │
│                                      [default: -1]                           │
│ --parser-detection          TEXT     Parser detection mode. [default: auto]  │
│ --n-jobs            -j      INTEGER  Number of parallel jobs. [default: -1]  │
│ --help              -h               Show this message and exit.             │
╰──────────────────────────────────────────────────────────────────────────────╯
```

**示例：**

```bash
uv run molop -q summary "tests/test_files/g16log/2-TS1-Opt.log" --out .sisyphus/tmp/doc_summary.csv --format csv --mode frame --frame -1 --n-jobs 1
```

```text
Summary written to .sisyphus/tmp/doc_summary.csv
```

### `transform`

将分子文件从一种格式转换为另一种格式（例如，将 Gaussian log 转换为 SDF）。

**何时使用：** 用于批量格式转换。

```bash
uv run molop transform --help
```

```text

 Usage: molop transform [OPTIONS] PATTERN

 Transform molecular files to another format.

╭─ Arguments ──────────────────────────────────────────────────────────────────╮
│ *    pattern      TEXT  File pattern to match. [required]                    │
╰──────────────────────────────────────────────────────────────────────────────╯
╭─ Options ────────────────────────────────────────────────────────────────────╮
│ *  --to                                  TEXT     Target format (xyz, sdf,   │
│                                                   cml, gjf, smi, orcainp).   │
│                                                   [required]                 │
│ *  --output-dir                          PATH     Directory to save          │
│                                                   transformed files.         │
│                                                   [required]                 │
│    --frame                               TEXT     Frame selection (e.g.,     │
│                                                   'all', '-1', '1,2').       │
│                                                   [default: -1]              │
│    --embed                 --no-embed             Whether to embed multiple  │
│                                                   frames in one file.        │
│                                                   [default: embed]           │
│    --parser-detection                    TEXT     Parser detection mode.     │
│                                                   [default: auto]            │
│    --n-jobs            -j                INTEGER  Number of parallel jobs.   │
│                                                   [default: -1]              │
│    --help              -h                         Show this message and      │
│                                                   exit.                      │
╰──────────────────────────────────────────────────────────────────────────────╯
```

#### 示例

```bash
uv run molop -q transform "tests/test_files/orca/h2_grad_orca.inp" --to sdf
--output-dir .sisyphus/tmp/doc_transform_out --frame -1 --embed
--parser-detection orcainp --n-jobs 1
```

```text
Transformation completed. Files saved to .sisyphus/tmp/doc_transform_out
```

### `visualize`

在网格图像（PNG 或 SVG）中可视化分子。

**何时使用：** 用于快速直观地检查一批分子的结构。

```bash
uv run molop visualize --help
```

```text

 Usage: molop visualize [OPTIONS] PATTERN

 Visualize molecules in a grid image.

╭─ Arguments ──────────────────────────────────────────────────────────────────╮
│ *    pattern      TEXT  File pattern to match. [required]                    │
╰──────────────────────────────────────────────────────────────────────────────╯
╭─ Options ────────────────────────────────────────────────────────────────────╮
│ *  --out                       PATH     Output file path (e.g., grid.png or  │
│                                         grid.svg).                           │
│                                         [required]                           │
│    --max-mols                  INTEGER  Maximum number of molecules to       │
│                                         visualize.                           │
│                                         [default: 16]                        │
│    --mols-per-row              INTEGER  Number of molecules per row.         │
│                                         [default: 4]                         │
│    --sub-img-width             INTEGER  Width of each sub-image.             │
│                                         [default: 200]                       │
│    --sub-img-height            INTEGER  Height of each sub-image.            │
│                                         [default: 200]                       │
│    --n-jobs            -j      INTEGER  Number of parallel jobs.             │
│                                         [default: -1]                        │
│    --parser-detection          TEXT     Parser detection mode.               │
│                                         [default: auto]                      │
│    --help              -h               Show this message and exit.          │
╰──────────────────────────────────────────────────────────────────────────────╯
```

**示例：**

```bash
uv run molop -q visualize "tests/test_files/g16log/2-TS1-Opt.log" --out .sisyphus/tmp/doc_grid.svg --max-mols 4 --mols-per-row 2 --sub-img-width 200 --sub-img-height 200 --n-jobs 1 --parser-detection auto
```

```text
Visualization saved to .sisyphus/tmp/doc_grid.svg
```

### `filter-by-codec`

根据检测到的编解码器 ID（例如 `g16log`, `xyz`）过滤文件。

**何时使用：** 用于分离不同格式的文件，或确保仅处理特定格式。

```bash
uv run molop filter-by-codec --help
```

```text

 Usage: molop filter-by-codec [OPTIONS] PATTERN

 Filter files by their detected codec ID.

╭─ Arguments ──────────────────────────────────────────────────────────────────╮
│ *    pattern      TEXT  File pattern to match. [required]                    │
╰──────────────────────────────────────────────────────────────────────────────╯
╭─ Options ────────────────────────────────────────────────────────────────────╮
│ *  --codec-id                  TEXT     Codec ID to filter by. [required]    │
│    --negate                             Negate the filter.                   │
│    --on-missing                TEXT     Behavior on missing codec ID (keep,  │
│                                         drop, error).                        │
│                                         [default: drop]                      │
│    --format            -f      TEXT     Output format (text or json).        │
│                                         [default: text]                      │
│    --parser-detection          TEXT     Parser detection mode.               │
│                                         [default: auto]                      │
│    --n-jobs            -j      INTEGER  Number of parallel jobs.             │
│                                         [default: -1]                        │
│    --help              -h               Show this message and exit.          │
╰──────────────────────────────────────────────────────────────────────────────╯
```

**示例：**

```bash
uv run molop -q filter-by-codec "tests/test_files/g16log/2-TS1-Opt.log" --codec-id g16log --format text --n-jobs 1
```

```text
/Users/tmj/Documents/proj/MolOP/tests/test_files/g16log/2-TS1-Opt.log
```

### `sample`

从匹配的文件中随机抽取 N 个文件。

**何时使用：** 用于选取具有代表性的文件子集进行测试或快速检查。

```bash
uv run molop sample --help
```

```text

 Usage: molop sample [OPTIONS] PATTERN

 Randomly sample N files from the matched files.

╭─ Arguments ──────────────────────────────────────────────────────────────────╮
│ *    pattern      TEXT  File pattern to match. [required]                    │
╰──────────────────────────────────────────────────────────────────────────────╯
╭─ Options ────────────────────────────────────────────────────────────────────╮
│ *  --n                         INTEGER  Number of files to sample.           │
│                                         [required]                           │
│    --seed                      INTEGER  Random seed.                         │
│    --format            -f      TEXT     Output format (text or json).        │
│                                         [default: text]                      │
│    --parser-detection          TEXT     Parser detection mode.               │
│                                         [default: auto]                      │
│    --n-jobs            -j      INTEGER  Number of parallel jobs.             │
│                                         [default: -1]                        │
│    --help              -h               Show this message and exit.          │
╰──────────────────────────────────────────────────────────────────────────────╯
```

**示例：**

```bash
uv run molop -q sample "tests/test_files/mix_format/*.gjf" --n 2 --seed 1 --format text --n-jobs 1
```

```text
/Users/tmj/Documents/proj/MolOP/tests/test_files/mix_format/S_Ph_Ni_TS.gjf
/Users/tmj/Documents/proj/MolOP/tests/test_files/mix_format/TS_Zy0fwX_ll_ad_14-19_15-16_optts_g16.gjf
```

### `stats`

显示给定文件中分子的统计信息。

**何时使用：** 用于获取数据集的高层概览（例如格式分布、计算状态）。

```bash
uv run molop stats --help
```

```text

 Usage: molop stats [OPTIONS] PATTERN

 Show statistics of the molecules in the given files.

╭─ Arguments ──────────────────────────────────────────────────────────────────╮
│ *    pattern      TEXT  File pattern to match. [required]                    │
╰──────────────────────────────────────────────────────────────────────────────╯
╭─ Options ────────────────────────────────────────────────────────────────────╮
│ --format            -f      TEXT     Output format (text or json).           │
│                                      [default: text]                         │
│ --parser-detection          TEXT     Parser detection mode. [default: auto]  │
│ --n-jobs            -j      INTEGER  Number of parallel jobs. [default: -1]  │
│ --help              -h               Show this message and exit.             │
╰──────────────────────────────────────────────────────────────────────────────╯
```

**示例：**

```bash
uv run molop -q stats "tests/test_files/g16log/2-TS1-Opt.log"
```

```text
Total files: 1

By format:
  g16log: 1

By state (last frame):
  normal: 1
```

### `groupby`

按键（例如 `detected_format_id`, `state`）对文件进行分组并输出分组结果。

**何时使用：** 用于将文件组织成类别以进行进一步处理。

```bash
uv run molop groupby --help
```

```text

 Usage: molop groupby [OPTIONS] PATTERN

 Group files by a key and output the groups.

╭─ Arguments ──────────────────────────────────────────────────────────────────╮
│ *    pattern      TEXT  File pattern to match. [required]                    │
╰──────────────────────────────────────────────────────────────────────────────╯
╭─ Options ────────────────────────────────────────────────────────────────────╮
│ --key                       TEXT     Key to group by (detected_format_id,    │
│                                      file_format, state).                    │
│                                      [default: detected_format_id]           │
│ --format            -f      TEXT     Output format (text or json).           │
│                                      [default: json]                         │
│ --parser-detection          TEXT     Parser detection mode. [default: auto]  │
│ --n-jobs            -j      INTEGER  Number of parallel jobs. [default: -1]  │
│ --help              -h               Show this message and exit.             │
╰──────────────────────────────────────────────────────────────────────────────╯
```

**示例：**

```bash
uv run molop -q groupby "tests/test_files/mix_format/*.gjf" --key detected_format_id --format json --n-jobs 1
```

```text
{
  "gjf": [
    "/Users/tmj/Documents/proj/MolOP/tests/test_files/mix_format/S_Ph_Ni_TS.gjf",
    "/Users/tmj/Documents/proj/MolOP/tests/test_files/mix_format/TS_Zy0fwX_ll_ad_14-19_15-16_optts_g16.gjf",
    "/Users/tmj/Documents/proj/MolOP/tests/test_files/mix_format/dsgdb9nsd_000007-6.gjf",
    "/Users/tmj/Documents/proj/MolOP/tests/test_files/mix_format/irc.gjf",
    "/Users/tmj/Documents/proj/MolOP/tests/test_files/mix_format/opt_point_charge_xtb.gjf"
  ]
}
```

### `chain`

旧版链式兼容模式。将所有参数转发给基于 Fire 的 CLI。

```bash
uv run molop chain --help
```

```text

 Usage: molop chain [OPTIONS]

 Legacy chain compatibility mode. Forwards all arguments to the Fire-based CLI.

╭─ Options ────────────────────────────────────────────────────────────────────╮
│ --help  -h        Show this message and exit.                                │
╰──────────────────────────────────────────────────────────────────────────────╯
```

---

## 旧版链式语法（高级）

旧版链式语法允许你在单个命令行中串联多个 MolOP 操作。这由 [Google Fire](https://github.com/google/python-fire) 提供支持。

### 操作指南

- **分隔符：** 使用单个连字符 `-` 来分隔链中的不同命令。
- **执行：** 每个链 **必须** 以 `end` 命令结尾，以触发执行并执行必要的清理。
- **引号：** 始终对通配符模式（例如 `"*.log"`）加引号，以防止 Shell 在它们到达 MolOP 之前对其进行展开。

### 链式示例

使用旧版链式语法生成 Gaussian 日志文件的 CSV 摘要。

```bash
uv run molop chain read "tests/test_files/g16log/2-TS1-Opt.log" - summary --output_path ".sisyphus/tmp/chain_ts_summary.csv" - end
```

```text
MolOP parsing with single process 100% ━━━━━━ 1/1  [ 0:00:00 < 0:00:… , ? it/s ]
MolOP processing frame summary with 1 jobs 100% ━━━━ 1/1  [ 0:0… < 0:0… , ?    ]
                                                                          it/s
INFO - Summary saved to .sisyphus/tmp/chain_ts_summary.csv
```

## 重要注意事项

- **Shell 展开：** 在 `zsh`（macOS 默认）等 Shell 中，未加引号且不匹配任何文件的通配符会导致错误（`no matches found`）。请始终使用引号：`"*.log"`。
- **安静模式：** 建议使用 `-q` 或 `--quiet` 以在示例和自动化脚本中获得稳定、干净的输出。
