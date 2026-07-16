# Select optimized results and transition states

Export stable optimized structures and transition-state candidates separately.

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

## Output

```text
{'parsed': N, 'normal': A, 'optimized': B, 'transition_states': C}
optimized.csv
transition_states.csv
```

`N/A/B/C` are counts from the current dataset. Use `Vibration.num_imaginary` in the transition-state
CSV to review imaginary-mode counts. A missing column means selected frames had no structured
frequency result.

## CLI

```bash
molop -q parse "results/*.log" \
  filter-state --state normal \
  filter-state --state ts \
  to-summary-df --full --out transition_states.csv
```

## Scientific review

`is_TS` is a programmatic candidate decision, not a replacement for checking the imaginary mode,
connectivity change, and reaction path. Database ingestion and automated reaction workflows should
review the displacement and endpoint chemistry.
