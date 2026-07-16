# 批量汇总

把一批计算输出整理成可筛选的 pandas DataFrame，并导出 CSV 或 JSON。

## 最短示例

```python
from molop import AutoParser

batch = AutoParser("results/*.out", n_jobs=4)
summary = batch.to_summary_df(
    frame=-1,
    brief=False,
    flatten_columns=True,
)
summary.to_csv("summary.csv", index=False)
print(summary.shape)
```

若 `results/` 中有 100 个成功解析的单任务输出，结果通常是 `(100, 列数)`；无法解析的文件
不会产生行。

## 用共享样例核对结果

```python
from molop import AutoParser

batch = AutoParser("water_mp2.out", n_jobs=1)
summary = batch.to_summary_df(brief=False, flatten_columns=True)
print(summary[["Status.IsNormal", "Energy.total_energy.hartree"]])
```

输出：

```text
   Status.IsNormal  Energy.total_energy.hartree
0             True                    -74.999375
```

## 最后一帧与全部 frame

```python
final_rows = batch.to_summary_df(frame=-1)
all_rows = batch.to_summary_df(frame="all")
selected_rows = batch.to_summary_df(frame=[0, -1])
```

`on_missing_frame="skip"` 默认跳过不存在的索引；需要把缺失 frame 当作错误时使用
`on_missing_frame="error"`。

## 简要与完整字段

```python
brief = batch.to_summary_df(brief=True, flatten_columns=True)
full = batch.to_summary_df(brief=False, flatten_columns=True)
print(brief.columns.tolist())
```

简要摘要稳定包含文件、分子、计算参数和状态列，例如：

```text
['DiskStorage.FilePath', ..., 'Status.IsNormal', 'Status.IsTS',
 'Status.IsOptimized']
```

完整摘要按实际存在结果增加 `Energy.*`、`Thermal.*`、`Vibration.*` 等列。

## 扁平列与 MultiIndex

```python
flat = batch.to_summary_df(flatten_columns=True)
print(flat["General.Charge"])

multi = batch.to_summary_df(flatten_columns=False)
print(multi[("General", "Charge", "")])
```

CSV/JSON 和普通分析脚本优先使用扁平列；需要明确分组、字段和单位三层语义时保留 MultiIndex。

## 导出 JSON

```python
summary.to_json("summary.json", orient="records", indent=2)
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

## 下一步

- [过滤与选择](filtering.md)
- [导出能量和热化学 CSV](../tutorials/energy-csv.md)
- [按科学性质查找字段](../reference/model_fields.md)
