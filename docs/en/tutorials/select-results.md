# Select optimized results and transition states

Export stable optimized structures and transition-state candidates separately.

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

## Output

??? example "Selection counts"

    ```text
    {'parsed': 1, 'normal': 1, 'optimized': 0, 'transition_states': 0}
    ```

The generated files are:

??? example "Created files"

    ```text
    optimized.csv
    transition_states.csv
    ```

The counts above are for the bundled single-point sample; counts and whether the two CSV files have
data depend on the input set. Use `Vibration.num_imaginary` in the transition-state CSV to review
imaginary-mode counts. A missing column means selected frames had no structured frequency result.

Replace `water_mp2.out` with `results/*.log` when selecting optimization or transition-state jobs.

## CLI

```bash
molop -q parse "results/*.log" \
  filter-state --state normal \
  filter-state --state ts \
  to-summary-df --full --out transition_states.csv
```

??? example "Created file"

    ```text
    transition_states.csv
    ```

## Scientific review

`is_TS` is a programmatic candidate decision, not a replacement for checking the imaginary mode,
connectivity change, and reaction path. Database ingestion and automated reaction workflows should
review the displacement and endpoint chemistry.
