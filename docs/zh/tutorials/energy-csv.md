# 导出能量和热化学 CSV

从每个文件的最终 frame 提取能量、零点能、焓和吉布斯自由能。

## Python

```python
from molop import AutoParser

batch = AutoParser("results/*.log")
full = batch.to_summary_df(
    frame=-1,
    brief=False,
    flatten_columns=True,
)

wanted = [
    "DiskStorage.FilePath",
    "Calc Parameter.Method",
    "Energy.total_energy.hartree",
    "Thermal.ZPVE.kilocalorie / mole",
    "Thermal.H_T.kilocalorie / mole",
    "Thermal.G_T.kilocalorie / mole",
]
energy_table = full.reindex(columns=wanted)
energy_table.to_csv("energies.csv", index=False)
print(energy_table.head().to_string(index=False))
```

## 输出

```text
energies.csv
```

表头包含请求的六列；没有频率/热化学区段的单点任务在 `Thermal.*` 列中显示空值，而不是零。

## 只保留有热化学的文件

```python
thermal_batch = batch.filter_state("thermal")
thermal_table = thermal_batch.to_summary_df(
    brief=False,
    flatten_columns=True,
).reindex(columns=wanted)
print(len(batch), len(thermal_batch))
```

输出两个计数，第二个是至少一个 frame 包含热化学结果的文件数。

## 注意

字段单位直接写在扁平列名中。温度、压力、标准态和理论水平仍需结合原始计算设置解释。
