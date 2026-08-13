# CLI 链式契约

本页解释 `molop parse` 的链式规则和终止操作边界。日常使用从[CLI 初体验](../getting_started/cli.md)
开始；参数完整列表见[CLI 命令参考](../command_line_interface.md)。

## 命令形状

```text
molop parse PATTERN [解析选项] 操作 [操作选项] ...
                                      -> 终止操作
```

`PATTERN` 可以是文件路径或 glob；可重复传入 `--input PATH_OR_GLOB` 追加输入范围。解析后，
返回 batch 的操作可以继续连接；返回摘要、渲染文本或路径 mapping 的操作必须位于末尾。CLI 会
在读取文件前检查这条规则。`--report` 也是终止解析模式，不能再连接操作命令。

以下示例使用随文档提供的 [water_mp2.out](../../assets/examples/water_mp2.out)。

## 常见链

=== "按 codec 筛选"

    ```bash
    molop -q parse "water_mp2.out" \
      --parser-detection orcaout \
      --n-jobs 1 \
      --output-format json \
      filter-by-codec --codec-id orcaout
    ```

    ??? example "输出形状"

        ```json
        [
          "<water_mp2.out 的绝对路径>"
        ]
        ```

=== "汇总"

    ```bash
    molop -q parse "water_mp2.out" --n-jobs 1 \
      filter-state --state normal \
      to-summary-df --full --out summary.csv
    ```

    ??? example "终端输出和生成文件"

        ```text
        Summary written to summary.csv
        summary.csv
        ```

=== "转换"

    ```bash
    molop -q parse "water_mp2.out" --n-jobs 1 \
      format-transform --format xyz --output-dir converted
    ```

    ??? example "生成文件"

        ```text
        converted/water_mp2.xyz
        ```

## 操作分类

| 类别 | 操作 | 返回状态 |
| --- | --- | --- |
| 可继续链接 | `filter-state`、`filter-value`、`filter-by-codec`、`sample` | 新 batch |
| 终止 | `to-summary-df` | CSV、JSON 或表格结果 |
| 终止 | `format-transform` | 渲染文本或生成文件 |
| 终止 | `draw-grid-image` | SVG/PNG 图像 |
| 终止 | `groupby`、`copy-to`、`move-to` | 路径分组或文件操作结果 |

终止操作不能再跟随其他操作。例如，下面的链会在解析前被拒绝：

```bash
molop parse "water_mp2.out" \
  to-summary-df --out summary.csv \
  filter-state --state normal
```

??? example "错误形状"

    ```text
    Error: to-summary-df returns non-FileBatchModelDisk and must be the last operation.
    ```

错误内容取决于 CLI 的具体校验分支；修复方式是把筛选放到 `to-summary-df` 之前。

## Writer 动态参数

`format-transform` 的通用参数之后可以接目标 writer 的参数。参数名由当前注册的 writer 提供：

```bash
molop parse "water_mp2.out" --n-jobs 1 \
  format-transform --format gjf --output-dir gaussian_inputs \
  --route-section "#p B3LYP/6-31G(d) opt" \
  --link0-commands "%nprocshared=8"
```

??? example "生成文件"

    ```text
    gaussian_inputs/water_mp2.gjf
    ```

查看当前安装版本的静态参数和 writer 动态参数：

```bash
molop parse PATTERN format-transform --help
```

??? example "帮助输出形状"

    ```text
    Usage: molop parse PATTERN format-transform [OPTIONS] [EXTRA_ARGS]...
    --format TEXT           Target writer format id.  [required]
    --output-dir DIRECTORY  Directory for generated files.
    --frame TEXT            Frame selection: all, int, or csv ints.
    ```

## Shell 补全

只查看补全脚本，不修改 shell 配置：

```bash
molop completion show --shell bash
```

??? example "输出形状"

    ```text
    # bash completion script
    <由当前 Click 命令树生成的脚本>
    ```

需要安装时运行 `molop completion install --shell bash`，然后重新加载 shell。安装命令会修改当前
用户的 shell 配置文件，适合明确的交互式环境；自动化环境使用 `completion show` 取得脚本内容。

## 相关页面

- [CLI 初体验](../getting_started/cli.md)
- [CLI 命令参考](../command_line_interface.md)
- [API 契约](../reference/api_contracts.md)
