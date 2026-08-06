# 导出能量和热化学 CSV

从每个文件的最终 frame 提取能量、零点能、焓和吉布斯自由能。

## Python

```python
from molop import AutoParser

batch = AutoParser("water_mp2.out", n_jobs=1)
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
print(energy_table.drop(columns="DiskStorage.FilePath").head().to_string(index=False))
```

## 输出

??? example "CSV 预览"

    ```text
    Calc Parameter.Method  Energy.total_energy.hartree  Thermal.ZPVE.kilocalorie / mole  Thermal.H_T.kilocalorie / mole  Thermal.G_T.kilocalorie / mole
                      MP2                  -74.999375                              NaN                             NaN                             NaN
    ```

文件表头包含请求的六列。以上数值来自随文档提供的 `water_mp2.out` 样例；该单点任务没有频率或
热化学区段，因此 `energies.csv` 中的 `Thermal.*` 单元格为空，而不是零。

## 只保留有热化学的文件

```python
thermal_batch = batch.filter_state("thermal")
thermal_table = thermal_batch.to_summary_df(
    brief=False,
    flatten_columns=True,
).reindex(columns=wanted)
print(len(batch), len(thermal_batch))
```

??? example "输出"

    ```text
    1 0
    ```

两个计数分别是成功解析文件数，以及至少一个 frame 包含热化学结果的文件数。将
`water_mp2.out` 替换为 `results/*.log`，即可处理一个结果目录。

## 注意

字段单位直接写在扁平列名中。温度、压力、标准态和理论水平仍需结合原始计算设置解释。
