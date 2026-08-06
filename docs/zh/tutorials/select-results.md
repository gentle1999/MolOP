# 筛选优化结果与过渡态

分别导出稳定优化结构和过渡态候选。

## Python

```python
from molop import AutoParser

batch = AutoParser("water_mp2.out", n_jobs=1)
normal = batch.filter_state("normal")
optimized = normal.filter_state("opt")
transition_states = normal.filter_state("ts")

print({
    "parsed": len(batch),
    "normal": len(normal),
    "optimized": len(optimized),
    "transition_states": len(transition_states),
})

optimized.to_summary_df(
    brief=False, flatten_columns=True
).to_csv("optimized.csv", index=False)

transition_states.to_summary_df(
    brief=False, flatten_columns=True
).to_csv("transition_states.csv", index=False)
```

## 输出

??? example "筛选计数"

    ```text
    {'parsed': 1, 'normal': 1, 'optimized': 0, 'transition_states': 0}
    ```

生成的文件为：

??? example "生成文件"

    ```text
    optimized.csv
    transition_states.csv
    ```

以上计数来自随文档提供的单点样例；实际计数以及两个 CSV 是否包含数据取决于输入文件。过渡态
CSV 的 `Vibration.num_imaginary` 可用于复核虚频数；若列缺失，说明选中 frame 没有结构化频率结果。

筛选优化任务或过渡态任务时，将 `water_mp2.out` 替换为 `results/*.log`。

## CLI

```bash
molop -q parse "results/*.log" \
  filter-state --state normal \
  filter-state --state ts \
  to-summary-df --full --out transition_states.csv
```

??? example "生成文件"

    ```text
    transition_states.csv
    ```

## 科学复核

`is_TS` 是程序化候选判断，不替代虚频方向、连接关系和反应路径检查。对数据库录入或自动反应
工作流，应继续核对振动模式及前后体。
