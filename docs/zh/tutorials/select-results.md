# 筛选优化结果与过渡态

分别导出稳定优化结构和过渡态候选。

## Python

```python
from molop import AutoParser

batch = AutoParser("results/*.log")
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

```text
{'parsed': N, 'normal': A, 'optimized': B, 'transition_states': C}
optimized.csv
transition_states.csv
```

`N/A/B/C` 是当前目录真实计数。过渡态 CSV 的 `Vibration.num_imaginary` 可用于复核虚频数；
若列缺失，说明选中 frame 没有结构化频率结果。

## CLI

```bash
molop -q parse "results/*.log" \
  filter-state --state normal \
  filter-state --state ts \
  to-summary-df --full --out transition_states.csv
```

## 科学复核

`is_TS` 是程序化候选判断，不替代虚频方向、连接关系和反应路径检查。对数据库录入或自动反应
工作流，应继续核对振动模式及前后体。
