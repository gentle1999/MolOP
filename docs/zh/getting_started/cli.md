# CLI 初体验

不编写 Python 脚本，直接汇总、筛选和转换一批文件。

## 查看命令

```bash
molop --help
molop parse --help
```

??? example "帮助摘要"

    ```text
    Usage: molop [OPTIONS] COMMAND [ARGS]...
    Commands: completion, parse
    Usage: molop parse [OPTIONS] PATTERN COMMAND1 [ARGS]...
    Commands: copy-to, draw-grid-image, filter-by-codec, filter-state,
              filter-value, format-transform, groupby, move-to, sample,
              to-summary-df
    ```

所有业务操作都从 `molop parse PATTERN` 开始。`PATTERN` 可以是单个路径或 shell glob；一条命令
需要覆盖多组输入时，可重复添加 `--input PATH_OR_GLOB`。

## 导出摘要 CSV

```bash
molop -q parse "results/*.out" \
  to-summary-df --full --out summary.csv
```

默认只汇总每个文件的最后一帧。加入 `--frame all` 可导出全部 frame。

对随文档提供的样例：

```bash
molop -q parse "water_mp2.out" --n-jobs 1 \
  to-summary-df --full --out summary.csv
```

??? example "输出"

    ```text
    Summary written to summary.csv
    ```

## 先筛选再导出

```bash
molop -q parse "results/*.log" \
  filter-state --state normal \
  filter-state --state opt \
  to-summary-df --full --out optimized.csv
```

操作从左到右执行。`filter-state` 保留 batch 状态，可继续连接下一个筛选或一个最终输出操作。

??? example "生成文件"

    ```text
    optimized.csv
    ```

## 导出最终结构

```bash
mkdir -p structures
molop -q parse "results/*.out" \
  filter-state --state normal \
  format-transform --format xyz --output-dir structures
```

指定 `--output-dir` 时会创建目录并写盘。要在终端预览而不写文件，使用 `--no-write`。

对 `water_mp2.out`，输出目录包含：

??? example "生成文件"

    ```text
    structures/water_mp2.xyz
    ```

## 处理扩展名歧义

```bash
molop parse "results/*.out" --parser-detection orcaout \
  to-summary-df --out summary.csv
```

??? example "生成文件"

    ```text
    summary.csv
    ```

## 下一步

- [CLI 常用任务](../guides/cli-recipes.md)
- [过滤与选择](../guides/filtering.md)
- [CLI 完整参考](../reference/cli.md)
