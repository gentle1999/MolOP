# 命令行接口

MolOP 只暴露一个业务命令：`molop parse`。

稳定的 CLI chain 契约，包括参数默认值和错误策略，见 [API 契约](reference/api_contracts.md)。

该命令先把文件解析成 `FileBatchModelDisk` 批状态，然后执行经过检查的操作链。返回值仍是
`FileBatchModelDisk` 的操作可以继续链接；返回其他结果的操作必须放在最后一步，且这一关系会在实际解析文件前检查。

为保证仓库内可复现，本文示例统一使用 `uv run molop ...`。如果已安装 MolOP，可直接使用
`molop ...`。

## 全局选项

- `-v`, `--verbose`: 启用详细输出。
- `-q`, `--quiet`: 启用安静模式。
- `--version`: 显示版本并退出。

## 查看命令

```bash
uv run molop --help
uv run molop parse --help
```

## Shell 补全

为当前 shell 安装补全：

```bash
molop completion install
```

可通过 `--shell bash`、`--shell zsh` 或 `--shell fish` 指定 shell。
安装后需要重启 shell，或手动 source 更新后的 rc 文件。

只查看生成的补全脚本、不安装：

```bash
molop completion show --shell bash
```

## 解析并链接操作

按检测到的 codec 过滤并输出路径：

```bash
uv run molop -q parse "tests/test_files/orca/single_point_inputs/h2_grad_orca.inp" \
  --parser-detection orcainp \
  --n-jobs 1 \
  --output-format json \
  filter-by-codec --codec-id orcainp
```

转换为其他格式：

```bash
uv run molop -q parse "tests/test_files/orca/single_point_inputs/h2_grad_orca.inp" \
  --parser-detection orcainp \
  --n-jobs 1 \
  format-transform --format xyz --output-dir .tmp/molop_xyz
```

为保持 CLI 兼容性，`--output-dir` 会隐式写盘。若不指定 `--output-dir`，可用
`--write` 写到每个源文件所在目录；也可用 `--no-write` 强制只渲染、不落盘，即使
命令中提供了输出目录。`format-transform` 写盘时不会把渲染后的文件内容打印到
stdout；生成的文件就是该操作的结果。

生成摘要表：

```bash
uv run molop -q parse "tests/test_files/orca/single_point_inputs/h2_grad_orca.inp" \
  --parser-detection orcainp \
  --n-jobs 1 \
  to-summary-df --out .tmp/molop_summary.csv
```

`to-summary-df` 默认汇总最后一帧（`--frame -1`）。使用 `--frame all`
可汇总所有帧，`--full` 可输出扩展字段。摘要 CSV/JSON 默认使用
`General.FrameID`、`Energy.total_energy.hartree` 这样的扁平列名；如需保留
三层 `(group, field, unit)` MultiIndex 列，使用 `--multi-index-columns`。

## 操作参考

可继续链接的操作：

| 操作 | 说明 |
| --- | --- |
| `filter-state` | 按计算状态过滤。 |
| `filter-value` | 按 charge、multiplicity 或文件格式过滤。 |
| `filter-by-codec` | 按检测到的 reader codec id 过滤。 |
| `sample` | 从当前 batch 中随机抽样。 |

终止操作：

| 操作 | 说明 |
| --- | --- |
| `format-transform` | 转换文件格式。`--output-dir` 隐式写盘；无输出目录时 `--write` 写到源目录；`--no-write` 禁用落盘。 |
| `to-summary-df` | 生成摘要表。 |
| `draw-grid-image` | 渲染分子网格图。 |
| `groupby` | 分组并输出路径。 |
| `copy-to` | 将当前 batch 文件复制到目录。 |
| `move-to` | 将当前 batch 文件移动到目录。 |

终止操作必须是操作链的最后一步。

## 格式相关动态参数

`format-transform` 在常规参数之后还可以接收 writer 特定参数：

```bash
uv run molop parse "input.log" \
  format-transform --format gjf --output-dir out \
  --route-section "#p B3LYP/6-31G(d) opt" \
  --link0-commands "%nprocshared=8"
```

这些动态参数来自已注册 writer 的元数据。安装 shell completion 后，补全可以根据当前
`--format` 给出可用参数名，并显示参数的简短含义。

静态参数面可通过 help 查看：

```bash
uv run molop parse PATTERN format-transform --help
```
