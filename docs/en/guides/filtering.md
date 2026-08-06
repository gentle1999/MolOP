# Filter and select

Select normally terminated, optimized, transition-state, or charge- and format-specific files from
a batch.

## Filter by status

```python
from molop import AutoParser

batch = AutoParser("water_mp2.out", n_jobs=1)
normal = batch.filter_state("normal")
optimized = normal.filter_state("opt")
transition_states = normal.filter_state("ts")

print(len(batch), len(normal), len(optimized), len(transition_states))
```

??? example "Output"

    ```text
    1 1 0 0
    ```

The four integers are total, normal, optimized, and transition-state file counts. Replace the input
with your own glob for a real batch; filtering returns a new batch and does not mutate the original.

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

??? example "Output"

    ```text
    1 1
    ```

## Charge, multiplicity, extension, and codec

```python
neutral = batch.filter_value("charge", 0)
doublets = batch.filter_value("multiplicity", 2)
logs = batch.filter_value("format", ".log")
print(len(neutral), len(doublets), len(logs))
```

`format` compares the source extension. Select by the detected reader with a codec ID:

```python
orca_outputs = batch.filter_by_codec_id("orcaout")
gaussian_logs = batch.filter_by_codec_id("g16log")
print(len(orca_outputs), len(gaussian_logs))
```

??? example "Output"

    ```text
    1 0 0
    1 0
    ```

## Custom scientific condition

```python
def has_one_imaginary_frequency(parsed_file):
    frame = parsed_file[-1]
    return bool(frame.vibrations and frame.vibrations.num_imaginary == 1)

one_imaginary = batch.filter_custom(has_one_imaginary_frequency)
print(len(one_imaginary))
```

??? example "Output"

    ```text
    0
    ```

## Inspect transition-state candidates

`filter_state("ts")` selects files with at least one frame whose `is_TS` value is true. Locate that
frame explicitly before inspecting its imaginary modes; the final frame is not necessarily the TS
frame in a multi-segment file:

```python
ts_batch = batch.filter_state("ts", n_jobs=-1)

for parsed_file in ts_batch:
    ts_frames = [frame for frame in parsed_file if frame.is_TS]
    for frame in ts_frames:
        imaginary = frame.vibrations.num_imaginary if frame.vibrations else None
        print(parsed_file.filename, frame.frame_id, imaginary)
```

??? example "Output shape (input-dependent)"

    ```text
    <source filename> <frame_id> <imaginary-mode count or None>
    ```

The printed count is data-dependent. A true `is_TS` value is a candidate decision; review the
imaginary-mode displacement and the reaction endpoints before using it as a confirmed transition
state.

## CLI equivalents

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

```bash
molop -q parse "results/*" \
  filter-by-codec --codec-id orcaout \
  filter-value --target charge --value 0 \
  to-summary-df --out neutral_orca.csv
```

??? example "Created file"

    ```text
    neutral_orca.csv
    ```

## Note

A status can be `None`, meaning the source did not provide enough evidence. In strict scientific
selection, keep unknown separate from explicit `False` and inspect representative raw outputs.

## Next steps

- [Select optimized results and transition states](../tutorials/select-results.md)
- [Transition-state analysis](../tutorials/transition-states.md)
- [Convert and export](conversion.md)
- [CLI task recipes](cli-recipes.md)
