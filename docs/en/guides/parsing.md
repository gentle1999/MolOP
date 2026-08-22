# Parse files

Use `AutoParser` with one file, a glob pattern, or a mixed list of paths.

For a single file without batch setup or worker scheduling, use `AutoFileParser`:

```python
from molop import AutoFileParser

parsed_file = AutoFileParser("water_mp2.out")
print(parsed_file.filename, len(parsed_file), parsed_file.detected_format_id)
```

The format is detected from the extension and source content, and the return
value is the file-level model itself.

## Shortest example

```python
from molop import AutoParser

batch = AutoParser("water_mp2.out", n_jobs=1)
for parsed_file in batch:
    print(parsed_file.filename, len(parsed_file), parsed_file.detected_format_id)
```

The returned `FileBatchModelDisk` contains files that parsed successfully and produced at least one
frame.
For the shared example, output is:

??? example "Output"

    ```text
    water_mp2.out 1 orcaout
    ```

## Input forms

```python
one = AutoParser("calculation.log")
many = AutoParser("results/*.out")
mixed = AutoParser([
    "gaussian/*.log",
    "orca/job.out",
    "structures/*.xyz",
])
```

MolOP expands globs itself. Inputs are converted to absolute paths, deduplicated, and sorted.

## Select frames

```python
parsed_file = batch[0]
first = parsed_file[0]
final = parsed_file[-1]

print([(frame.frame_id, len(frame.atoms)) for frame in parsed_file])
```

??? example "Output"

    ```text
    [(0, 3)]
    ```

When only final results matter and memory use is important:

```python
batch = AutoParser("results/*.log", only_last_frame=True)
```

`only_last_frame=True` changes which frames the parser retains. It is not equivalent to parsing the
full trajectory and then indexing `[-1]`; leave it off when you need the optimization history.

## Automatic and explicit format detection

The default `parser_detection="auto"` chooses reader candidates from the extension and probes the
content. For ambiguous extensions, provide a format ID:

```python
orca = AutoParser("job.out", parser_detection="orcaout")
xtb = AutoParser("xtb.out", parser_detection="xtbout")
fchk = AutoParser("molecule.fchk", parser_detection="g16fchk")
```

See the [format overview](../reference/format_support.md) for format IDs.

## Parallel and structure-only parsing

```python
batch = AutoParser("results/*")

structures = AutoParser(
    "results/*.log",
    n_jobs=-1,
    only_extract_structure=True,
)
```

`AutoParser` and all batch operations default to `n_jobs=-1`. MolOP resolves this to the smaller of
the process-aware CPU limit and the global `molopconfig.max_jobs` ceiling. With the default
`max_jobs=None`, the limit follows observable scheduler/container quotas and CPU affinity where the
platform exposes them. Set a positive `max_jobs` for a process-wide ceiling. Pass a positive `n_jobs`
to request a specific worker count, or use `n_jobs=1` for a small batch or debugging. The CLI parse
and operation commands follow the same independent `-1` default.

Parsing itself does not use the MolGR safety budget: `AutoParser` keeps the full `effective_max_jobs`
limit. The stricter two-thirds budget applies only after a graph-dependent operation enters MolGR.

`only_extract_structure=True` skips non-structural results and is unsuitable for extracting energy
or thermochemistry.

## Find failed inputs

```python
from pathlib import Path
from molop import AutoParser

inputs = sorted(Path("results").glob("*.out"))
batch = AutoParser(inputs, n_jobs=1)
parsed = {Path(path).resolve() for path in batch.file_paths}
failed = [path for path in inputs if path.resolve() not in parsed]

for path in failed:
    print("not added to batch:", path)
```

When every input succeeds this prints nothing. A missing result produces one line such as:

??? example "Output when a file is omitted"

    ```text
    not added to batch: results/truncated.out
    ```

Missing files, unsupported formats, and files without usable frames are omitted and reported in the
MolOP log. Retry with `n_jobs=1`, then check extensions, encoding, and format IDs in
[Troubleshooting](troubleshooting.md).

## Learn more

- [Read calculation results](results.md)
- [Batch summaries](batch.md)
- [Format overview](../reference/format_support.md)
- [Source evidence and serialization](../reference/api_contracts.md)
