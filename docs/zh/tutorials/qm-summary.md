# 汇总 Gaussian、ORCA 与 xTB

把不同程序的输出放入同一个 batch，用统一字段导出最终状态和能量。

## 输入

```text
results/
  gaussian_opt.log
  orca_sp.out
  xtb.out
```

`.out` 扩展名可能同时匹配 ORCA 和 xTB，自动模式会继续探测文件内容。

## Python

```python
from molop import AutoParser

batch = AutoParser("results/*", n_jobs=3)
summary = batch.to_summary_df(
    frame=-1,
    brief=False,
    flatten_columns=True,
)

columns = [
    "Calc Parameter.Software",
    "Calc Parameter.Method",
    "Status.IsNormal",
    "Energy.total_energy.hartree",
]
result = summary.reindex(columns=columns)
print(result.to_string(index=False))
result.to_csv("qm_summary.csv", index=False)
```

## 输出

`qm_summary.csv` 每个成功解析文件一行。对随文档提供的样例，选中列打印为：

??? example "输出"

    ```text
    Calc Parameter.Software Calc Parameter.Method Status.IsNormal Energy.total_energy.hartree
                       ORCA                   MP2            True                  -74.999375
    ```

上面的稳定输出使用共享的单文件 ORCA 样例。处理混合目录时，每个成功解析文件占一行，缺失
字段保留为空值；文件路径等环境相关列应单独处理，不应把绝对路径当作固定输出。不要直接比较
不同方法或不同哈密顿量的 `total_energy`。

不准备混合目录时，可用下面的代码复现上述行：

```python
sample = AutoParser("water_mp2.out", n_jobs=1)
sample_summary = sample.to_summary_df(
    frame=-1,
    brief=False,
    flatten_columns=True,
)
sample_columns = [
    "Calc Parameter.Software",
    "Calc Parameter.Method",
    "Status.IsNormal",
    "Energy.total_energy.hartree",
]
print(sample_summary.reindex(columns=sample_columns).to_string(index=False))
```

??? example "输出"

    ```text
    Calc Parameter.Software Calc Parameter.Method Status.IsNormal Energy.total_energy.hartree
                       ORCA                   MP2            True                  -74.999375
    ```

## 按程序检查

```python
for parsed_file in batch:
    print(parsed_file.filename, parsed_file.detected_format_id)
```

典型输出：

??? example "输出"

    ```text
    gaussian_opt.log g16log
    orca_sp.out orcaout
    xtb.out xtbout
    ```

格式的精确能力见 [Gaussian log](../reference/formats/g16log.md)、
[ORCA output](../reference/formats/orcaout.md) 和 [xTB output](../reference/formats/xtbout.md)。
