# CLI 常用任务

以下命令面向已经安装 MolOP 的用户，均从 `molop parse PATTERN` 开始。

## 汇总最后一帧

```bash
molop -q parse "results/*.log" \
  to-summary-df --full --out summary.csv
```

结果：当前目录生成 `summary.csv`，每个成功解析文件一行，列名类似
`Status.IsNormal` 和 `Energy.total_energy.hartree`。

## 汇总全部 frame

```bash
molop -q parse "results/*.log" \
  to-summary-df --frame all --full --out trajectory.csv
```

结果：`trajectory.csv` 中每个 frame 一行，`General.FrameID` 标识 frame。

## 只保留正常优化结果

```bash
molop -q parse "results/*.log" \
  filter-state --state normal \
  filter-state --state opt \
  to-summary-df --full --out optimized.csv
```

结果：`optimized.csv` 只包含两个筛选都通过的文件。

## 选择 ORCA 输出

```bash
molop -q parse "results/*" \
  filter-by-codec --codec-id orcaout \
  to-summary-df --out orca.csv
```

`filter-by-codec` 使用实际 reader ID，不依赖扩展名字符串。

## 导出最终 XYZ

```bash
molop -q parse "results/*.out" \
  filter-state --state normal \
  format-transform --format xyz --output-dir structures
```

结果：`structures/` 下按源文件主名生成 `.xyz` 文件。输出目录不存在时 CLI 会创建。

## 预览而不写盘

```bash
molop -q parse "water_mp2.out" \
  format-transform --format xyz --no-write
```

终端输出包含源路径和渲染后的 XYZ 文本；不会生成 `.xyz` 文件。

## 输出 JSON

```bash
molop -q parse "water_mp2.out" --output-format json \
  filter-state --state normal \
  sample --n 1 --seed 1
```

输出是路径数组：

```json
[
  "/absolute/path/to/water_mp2.out"
]
```

## 规则

- 操作按书写顺序执行。
- `filter-state`、`filter-value`、`filter-by-codec` 和 `sample` 可继续链接。
- `to-summary-df`、`format-transform`、`copy-to`、`move-to`、`groupby` 和
  `draw-grid-image` 是最终操作，之后不能再接业务操作。
- 使用 `molop parse --help` 查看操作列表，使用具体操作的 `--help` 查看参数。

完整约束见 [CLI 命令参考](../command_line_interface.md)。
