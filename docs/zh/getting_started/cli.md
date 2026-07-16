# CLI 初体验

不编写 Python 脚本，直接汇总、筛选和转换一批文件。

## 查看命令

```bash
molop --help
molop parse --help
```

所有业务操作都从 `molop parse PATTERN` 开始。`PATTERN` 可以是单个路径或 shell glob。

## 导出摘要 CSV

```bash
molop -q parse "results/*.out" \
  to-summary-df --full --out summary.csv
```

默认只汇总每个文件的最后一帧。加入 `--frame all` 可导出全部 frame。

## 先筛选再导出

```bash
molop -q parse "results/*.log" \
  filter-state --state normal \
  filter-state --state opt \
  to-summary-df --full --out optimized.csv
```

操作从左到右执行。`filter-state` 保留 batch 状态，可继续连接下一个筛选或一个最终输出操作。

## 导出最终结构

```bash
mkdir -p structures
molop -q parse "results/*.out" \
  filter-state --state normal \
  format-transform --format xyz --output-dir structures
```

指定 `--output-dir` 时会创建目录并写盘。要在终端预览而不写文件，使用 `--no-write`。

## 处理扩展名歧义

```bash
molop parse "results/*.out" --parser-detection orcaout \
  to-summary-df --out summary.csv
```

## 下一步

- [CLI 常用任务](../guides/cli-recipes.md)
- [过滤与选择](../guides/filtering.md)
- [CLI 完整参考](../command_line_interface.md)
