# 批量汇总

把一批计算输出整理成可筛选的 pandas DataFrame，并导出 CSV 或 JSON。

## 最短示例

```python
from molop import AutoParser

batch = AutoParser("water_mp2.out", n_jobs=1)
summary = batch.to_summary_df(
    frame=-1,
    brief=False,
    flatten_columns=True,
)
summary.to_csv("summary.csv", index=False)
print(summary.shape)
```

??? example "输出"

    ```text
    (1, 23)
    ```

把 `water_mp2.out` 换成 `results/*.out` 即可处理一批文件；每个成功解析文件一行，无法解析的
文件不会产生行。

## 批处理并行

解析默认使用自动并行。Python API 中的 batch 操作默认使用一个 worker；需要让汇总、筛选或
转换也使用自动并行上限时，显式传入 `n_jobs=-1`：

```python
from molop import AutoParser

batch = AutoParser("results/*.out")
normal = batch.filter_state("normal", n_jobs=-1)
summary = normal.to_summary_df(
    frame=-1,
    brief=False,
    flatten_columns=True,
    n_jobs=-1,
)
summary.to_csv("normal.csv", index=False)
```

`n_jobs=-1` 受 `molopconfig.max_jobs` 限制；需要更低上限时，在创建 batch 前修改全局配置。
排查 parser 或原生库问题时使用 `n_jobs=1`。CLI 会把 parse 级 `--n-jobs` 传给后续操作。

需要对每个文件执行自定义操作时，使用 `parallel_execute` 并显式返回结果：

```python
def file_format(parsed_file):
    return parsed_file.filename, parsed_file.detected_format_id

rows = batch.parallel_execute(file_format, n_jobs=-1, return_results=True)
print(list(rows))
```

对随文档提供的样例，结果为：

??? example "输出"

    ```text
    [('water_mp2.out', 'orcaout')]
    ```

## 用共享样例核对结果

```python
import pandas as pd
from IPython.display import display
from molop import AutoParser

batch = AutoParser("water_mp2.out", n_jobs=1)
summary = batch.to_summary_df(brief=False, flatten_columns=True)
with pd.option_context("display.max_columns", None, "display.max_rows", None, "display.width", 240):
    display(summary)
```

### 完整表格输出

<!-- notebook-output: examples/02-batch-summary-filter-select.ipynb#batch-summary -->

表格来自 [Notebook 02](../examples/02-batch-summary-filter-select.ipynb) 的已执行输出，包含
`summary` 的全部 23 列，没有临时切片。CI 先执行 Notebook，再由 MkDocs 嵌入保存的 HTML。

## 最后一帧与全部 frame

```python
final_rows = batch.to_summary_df(frame=-1)
all_rows = batch.to_summary_df(frame="all")
selected_rows = batch.to_summary_df(frame=[0, -1])
print(final_rows.shape, all_rows.shape, selected_rows.shape)
```

??? example "输出"

    ```text
    (1, 19) (1, 19) (2, 19)
    ```

`on_missing_frame="skip"` 默认跳过不存在的索引；需要把缺失 frame 当作错误时使用
`on_missing_frame="error"`。

## 简要与完整字段

```python
brief = batch.to_summary_df(brief=True, flatten_columns=True)
full = batch.to_summary_df(brief=False, flatten_columns=True)
print(brief.columns[:7].tolist())
```

简要摘要稳定包含文件、分子、计算参数和状态列。共享样例的前七列为：

??? example "列名"

    ```text
    ['DiskStorage.FilePath', 'DiskStorage.FileFormat', 'General.Charge',
     'General.Multiplicity', 'General.CanonicalSMILES', 'General.NumAtoms',
     'General.FrameID']
    ```

完整摘要按实际存在结果增加 `Energy.*`、`Thermal.*`、`Vibration.*` 等列。

## 扁平列与 MultiIndex

```python
flat = batch.to_summary_df(flatten_columns=True)
print(flat["General.Charge"].tolist())

multi = batch.to_summary_df(flatten_columns=False)
print(multi[("General", "Charge", "")].tolist())
```

??? example "输出"

    ```text
    [0]
    [0]
    ```

CSV/JSON 和普通分析脚本优先使用扁平列；需要明确分组、字段和单位三层语义时保留 MultiIndex。

## 分组、抽样与预览

这些操作仍然返回或使用 batch，适合在正式导出前检查数据：

```python
groups = batch.groupby(lambda parsed_file: parsed_file.detected_format_id, n_jobs=-1)
print({key: len(group) for key, group in groups.items()})

sample = batch.sample(n=1, seed=1)
print(sample.file_names)

grid = batch.draw_grid_image(maxMols=16, useSVG=True, n_jobs=1)
print(type(grid).__name__, grid.lstrip()[:4])
```

对随文档提供的样例，文本输出为：

??? example "输出"

    ```text
    {'orcaout': 1}
    ['water_mp2.out']
    str <svg
    ```

`useSVG=True` 时，`grid` 是 SVG 字符串，可直接嵌入 Notebook 或写入 `.svg` 文件。MolOP 默认
使用 `rdkit-dof` 进行景深绘制；设置 `molopconfig.use_dof_effect_drawer = False` 可切换到标准
RDKit 绘制器。设置 `useSVG=False` 时返回栅格图像对象。CLI 等价操作是
`draw-grid-image --out structures.svg` 或 `draw-grid-image --out structures.png`。

## 导出 JSON

```python
summary.to_json("summary.json", orient="records", indent=2)
print("summary.json")
```

??? example "生成文件"

    ```text
    summary.json
    ```

## 大批量建议

- 先用 5-20 个代表性文件和 `n_jobs=1` 检查字段。
- 确认格式后再提高 `n_jobs`；机械硬盘或超大文件不一定随进程数线性加速。
- 只要最终结果时使用 `only_last_frame=True`，但不要用于优化轨迹分析。
- 用输入路径集合与 `batch.file_paths` 的差集记录未解析文件。
- 对不同理论水平分组后再比较 `total_energy`。

## CLI 等价操作

```bash
molop -q parse "results/*.out" \
  to-summary-df --full --flatten-columns --out summary.csv
```

??? example "终端输出"

    ```text
    Summary written to summary.csv
    ```

## 下一步

- [过滤与选择](filtering.md)
- [导出能量和热化学 CSV](../tutorials/energy-csv.md)
- [按科学性质查找字段](../reference/model_fields.md)
- [CLI 常用任务](cli-recipes.md)
