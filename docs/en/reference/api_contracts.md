# API Contracts

Use this page as a compact reference for the public batch, frame-selection,
conversion, summary, and CLI contracts. Parser implementation rules belong in
the [parser contract](../developer/parser-contract.md).

## At a glance

| Need | API | Result |
| --- | --- | --- |
| Parse one file | [`AutoFileParser(...)`](#autofileparser) | File-level model |
| Parse loaded text or bytes | [`AutoMemoryParser(...)`](#automemoryparser) | Memory file-level model |
| Parse files and globs | [`AutoParser(...)`](#autoparser) | `FileBatchModelDisk` |
| Select frames | [Frame selector](#frame-selector) | Normalized frame indices |
| Convert formats | [`format_transform(...)`](#format_transform) | Rendered text or path mapping |
| Render trajectories | [`draw_animation(...)`](#draw_animation) | GIF or animated SVG |
| Export TS endpoints | [`save_pre_post_ts(...)`](#save_pre_post_ts) | Endpoint path mapping |
| Build a table | [`to_summary_df(...)`](#to_summary_df) | `pandas.DataFrame` |
| Run one function per file | [`parallel_execute(...)`](#parallel_execute) | Results or `None` |
| Compose from the shell | [CLI command reference](../command_line_interface.md) | A terminal operation result |

## `AutoFileParser`

```python
from molop import AutoFileParser

parsed_file = AutoFileParser("water_mp2.out")
```

`AutoFileParser` accepts one existing file path, probes reader candidates, and
returns the file-level model directly. It does not expand globs, construct a
batch parser, or schedule worker processes. This makes it suitable for callers
that own the surrounding process or task model. `parser_detection="auto"` is
the default; explicit format IDs are also accepted. Missing files, unsupported
formats, format mismatches, and parse failures raise exceptions.

## `AutoMemoryParser`

```python
from molop import AutoBytesParser, AutoTextParser

parsed_text = AutoTextParser("1\nwater\nH 0.0 0.0 0.0\n", parser_detection="xyz")
parsed_bytes = AutoBytesParser(raw_bytes, parser_detection="xyz")
```

`AutoMemoryParser` accepts already-loaded text, bytes-like data, or a text/binary
stream and returns a memory-backed file model. `AutoTextParser` and
`AutoBytesParser` make the input type explicit. Use `parser_detection` when the
source has no filename extension; `"auto"` tries the registered in-memory
readers and their content probes. Binary parsing decodes with `source_encoding`
and preserves exact source bytes when `capture_source_evidence=True`.

## `AutoParser`

```python
from molop import AutoParser

batch = AutoParser("water_mp2.out", n_jobs=1)
```

The input is one path, `os.PathLike`, or a one-dimensional iterable of paths and
glob patterns. MolOP expands globs, converts matches to absolute paths,
deduplicates them, and sorts them by normalized path. Input grouping and input
order are not retained.

Unmatched globs add no files. An empty iterable returns an empty batch. A
missing literal path is passed to the batch parser and then warned about and
skipped. Nested iterables and non-path members raise `TypeError` with the input
index. Parse options apply to every expanded path.

## Frame selector {#frame-selector}

The same selector is used by Python conversion and summary APIs:
`int | Sequence[int] | "all"`.

| Value | Meaning |
| --- | --- |
| `-1` | Last frame; negative integers count from the end. |
| `0`, `[0, 2]` | One frame or an explicit list of frames. |
| `"all"` | Every frame. |

Strings other than `"all"` raise `ValueError`; a sequence containing a
non-integer raises `TypeError`. For summaries, `on_missing_frame="skip"` drops
out-of-range indices and `"error"` raises `IndexError`. `slice` is not part of
this API contract.

## `format_transform`

Use the object level that matches the desired output:

=== "Python"

    ```python
    rendered = batch.format_transform(
        "xyz",
        frame=-1,
        embed_in_one_file=True,
        write_to_disk=False,
        n_jobs=1,
    )
    print(type(rendered).__name__, len(rendered))
    ```

=== "CLI"

    ```bash
    molop -q parse "water_mp2.out" \
      format-transform --format xyz --frame -1 --no-write
    ```

??? example "Return shape"

    ```text
    dict 1
    ```

    For one input file, the Python call returns a `dict` containing the source path and rendered
    text; the CLI prints the rendered text without creating an `.xyz` file. Replace the sample path
    with `results/*.out` when processing a glob.

| Parameter | Applies to | Default | Contract |
| --- | --- | --- | --- |
| `format` | frame, file, batch | required | Registered writer id, such as `xyz`, `sdf`, `gjf`, or `cml`. |
| `frame` | file, batch | `-1` | `int`, integer sequence, or `"all"`. |
| `embed_in_one_file` | file, batch | `True` | Combine selected frames when the writer supports it. |
| `write_to_disk` | frame, file, batch | `False` | Write output; disk-backed objects use a source-derived path when no path is supplied. |
| `file_path` / `output_dir` | frame/file / batch | `None` | Destination used only when writing. A memory-only object needs `file_path`; a batch `output_dir` must already exist in the Python API. |
| `n_jobs` | batch | `1` | Number of workers for the batch operation. |
| `**kwargs` | frame, file, batch | - | Writer-specific options, such as Gaussian route or link-0 settings. |

Returns `str` for one frame, `str` or `list[str]` for a file, and
`dict[str, str | list[str]]` for a batch. Unsupported formats, writer
validation errors, and write errors propagate to the caller. Batch conversion
logs a warning and returns an empty value for a file whose conversion fails;
inspect the log when a mapping contains an empty string or list.

See [transform behavior](transform_behavior.md) for writer domains, frame
composition, and format-specific details.

## `draw_animation`

```python
animation = parsed_file.draw_animation(
    image_format="gif",
    file_path="trajectory.gif",
    duration=120,
)
```

`draw_animation` renders every frame with a usable RDKit graph. Invalid frames are skipped, and
default legends retain frame IDs, TS status, and total energy when available. `image_format` accepts
`"gif"` or `"svg"`; `duration` may be one positive integer or one value per original frame. When
frames are skipped, frame-aligned sequences such as `duration`, `legends`, and highlight lists are
filtered to the retained frames.

## `save_pre_post_ts`

```python
file_exports = parsed_file.save_pre_post_ts("ts-endpoints", format="sdf")
batch_exports = batch.save_pre_post_ts("ts-endpoints-batch", format="sdf", n_jobs=2)
```

The file-level method returns `{frame_id: (pre_path, post_path)}` and the batch-level method adds
the source path as its outer key. Only calculation files with this operation are exported; other
files in a mixed batch are skipped with a warning. Unique source stems are retained, while duplicate
stems receive a stable source-path digest to prevent overwriting another file's results. The
endpoints are geometry-based candidates, not optimized reactant or product structures. Endpoint
inference samples `steps=7` amplitudes per displacement side from `min_ratio=0.75` through
`max_ratio=1.75`, selects each side's most frequent topology, and maps the candidate with more
fragments to the precursor.

## `to_summary_df`

```python
summary = batch.to_summary_df(
    mode="frame",
    frame=-1,
    brief=False,
    flatten_columns=True,
    on_missing_frame="skip",
)
print(summary.shape)
```

??? example "Output"

    ```text
    (1, 23)
    ```

| Parameter | Default | Contract |
| --- | --- | --- |
| `mode` | `"frame"` | Summarize selected frames or one row per file with `"file"`. |
| `frame` | `-1` | Frame selector; used in frame mode. |
| `brief` | `True` | Compact fields when `True`; request expanded result fields with `False`. |
| `flatten_columns` | `False` | Convert the three-level column index to names such as `Energy.total_energy.hartree`. |
| `on_missing_frame` | `"skip"` | Drop unavailable frames or raise `IndexError` with `"error"`. |
| `n_jobs` | `1` | Number of workers for batch summary extraction. |
| `**kwargs` | - | Extra options forwarded to `to_summary_series`. |

The return value is a `pandas.DataFrame`; no selected rows produces an empty
DataFrame. The [batch guide](../guides/batch.md) and [batch summary
Notebook](../examples/02-batch-summary-filter-select.ipynb) show the complete
table, including all columns generated by the shared example.

## `parallel_execute`

```python
def describe(parsed_file):
    return parsed_file.filename, parsed_file.detected_format_id

rows = batch.parallel_execute(
    describe,
    n_jobs=1,
    return_results=True,
)
print(list(rows))
```

??? example "Output shape"

    ```text
    [('water_mp2.out', 'orcaout')]
    ```

MolOP invokes `func(diskfile, *args, **kwargs)` once for each file in the batch.
`n_jobs` controls execution, `desc` labels progress, and `return_as` accepts
Joblib's `"list"`, `"generator"`, or `"generator_unordered"` modes. Set
`return_results=False` for side-effect-only work; otherwise results are
returned as an iterable. Exceptions raised by the callable or Joblib propagate.

`_diskfiles_snapshot` is an internal alignment hook for MolOP operations and is
not part of normal application code.

## CLI and serialization

- See the [CLI command reference](../command_line_interface.md) for parsing, filtering, summary,
  and conversion chains.
- See [Source evidence and serialization](serialization.md) for exporting file metadata and frame
  records to a database or audit system.

## Errors and boundaries

| Situation | Behavior |
| --- | --- |
| Unsupported input or parser | File is skipped with parser diagnostics; an explicit invalid selector or option raises. |
| Invalid frame selector | `ValueError` for invalid selector strings; `TypeError` for non-integer sequence members. |
| Missing summary frame | Skip by default; `on_missing_frame="error"` raises `IndexError`. |
| Unsupported writer or writer validation failure | Error propagates for object-level conversion. |
| Batch conversion failure | Warning plus an empty string/list for that source entry. |

## Related pages

- [Parsing files](../guides/parsing.md)
- [Batch summaries](../guides/batch.md)
- [Transform behavior](transform_behavior.md)
- [Source evidence and serialization](serialization.md)
- [Format support](format_support.md)
- [Parser contract](../developer/parser-contract.md)
