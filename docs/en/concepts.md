# Understand Batch, File, and Frame

Using MolOP requires only three levels.

```text
FileBatchModelDisk
  ├─ file 0
  │   ├─ frame 0
  │   └─ frame 1
  └─ file 1
      └─ frame 0
```

## Batch

`AutoParser(...)` returns a batch. It stores parsed files and provides summary, filter, group, and
batch conversion operations.

```python
batch = AutoParser("results/*.log")
print(len(batch), batch.file_names)
```

Output is a count and a filename list, for example:

```text
2 ['reactant.log', 'product.log']
```

## File

`batch[0]` selects one source file. A file can hold one single-point frame or an optimization
trajectory, frequency step, and multiple calculation segments.

```python
parsed_file = batch[0]
print(parsed_file.filename, len(parsed_file), parsed_file.detected_format_id)
```

Example output:

```text
reactant.log 12 g16log
```

## Frame

A frame is one structure/result snapshot. Most tasks read the final frame:

```python
final = parsed_file[-1]
print(final.frame_id, final.charge, final.multiplicity)
```

Iterate all frames for trajectories and do not assume every frame has the same result containers.

## The final frame is not always the best structure

Error termination, Link1, multiple segments, and frequency tasks can make the final frame differ
from the converged optimization frame. Use `filter_state("opt")` and
`parsed_file.closest_optimized_frame`, then review the source when correctness matters.

## Learn more

Registry, codec, source lifecycle, and internal model boundaries are developer concepts. See the
[Architecture overview](developer/overview.md) and [API contracts](reference/api_contracts.md).
