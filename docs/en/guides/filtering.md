# Filter and select

Select normally terminated, optimized, transition-state, or charge- and format-specific files from
a batch.

## Filter by status

```python
from molop import AutoParser

batch = AutoParser("results/*")
normal = batch.filter_state("normal")
optimized = normal.filter_state("opt")
transition_states = normal.filter_state("ts")

print(len(batch), len(normal), len(optimized), len(transition_states))
```

The four integers are total, normal, optimized, and transition-state file counts. Filtering returns
a new batch and does not mutate the original.

| State | Files retained |
| --- | --- |
| `normal` | All frames show no error and the format provides normal-status evidence |
| `error` | At least one frame shows an error |
| `opt` | Optimization status exists and the closest optimized frame is converged |
| `ts` | At least one frame has a true `is_TS` value |
| `thermal` | At least one frame has thermochemistry |
| `no-img` | Frequency data exists and all relevant frames have zero imaginary modes |

Invert a selection with:

```python
not_normal = batch.filter_state("normal", negate=True)
```

## Check the shared example

```python
batch = AutoParser("water_mp2.out", n_jobs=1)
print(len(batch), len(batch.filter_state("normal")))
```

Output:

```text
1 1
```

## Charge, multiplicity, extension, and codec

```python
neutral = batch.filter_value("charge", 0)
doublets = batch.filter_value("multiplicity", 2)
logs = batch.filter_value("format", ".log")
```

`format` compares the source extension. Select by the detected reader with a codec ID:

```python
orca_outputs = batch.filter_by_codec_id("orcaout")
gaussian_logs = batch.filter_by_codec_id("g16log")
```

## Custom scientific condition

```python
def has_one_imaginary_frequency(parsed_file):
    frame = parsed_file[-1]
    return bool(frame.vibrations and frame.vibrations.num_imaginary == 1)

one_imaginary = batch.filter_custom(has_one_imaginary_frequency)
```

## CLI equivalents

```bash
molop -q parse "results/*.log" \
  filter-state --state normal \
  filter-state --state ts \
  to-summary-df --full --out transition_states.csv
```

```bash
molop -q parse "results/*" \
  filter-by-codec --codec-id orcaout \
  filter-value --target charge --value 0 \
  to-summary-df --out neutral_orca.csv
```

## Note

A status can be `None`, meaning the source did not provide enough evidence. In strict scientific
selection, keep unknown separate from explicit `False` and inspect representative raw outputs.

## Next steps

- [Select optimized results and transition states](../tutorials/select-results.md)
- [Convert and export](conversion.md)
- [CLI task recipes](cli-recipes.md)
