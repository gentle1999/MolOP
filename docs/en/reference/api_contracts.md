# API Contracts

This page fixes the public contract for the transform, summary, parallel
execution, and CLI chain surfaces.

## Frame Selector

MolOP uses one frame selector concept for public transform and summary APIs:
`int | Sequence[int] | "all"`.

| Surface | Parameter | Default | Return mapping | Error strategy |
| --- | --- | --- | --- | --- |
| Python transform | `frame` | `-1` | Normalized to frame indices before rendering. Negative integers count from the end. | Non-integer sequence items raise `TypeError`; strings other than `"all"` raise `ValueError`. |
| Python summary | `frame` | `-1` | Normalized per file. `"all"` expands to every frame in each file. | Same selector errors; missing frames are controlled by `on_missing_frame`. |
| CLI | `--frame` | `-1` | Parsed from `"all"`, one integer, or comma-separated integers, then forwarded to the Python API. | Invalid text raises a CLI usage error before executing the chain. |

`slice` is not part of the public frame selector contract. Sequence slicing on
file objects remains normal Python collection behavior.

## `format_transform`

### Python API

| Parameter | Applies to | Default | Meaning | Error strategy |
| --- | --- | --- | --- | --- |
| `format` | frame, file, batch | Required | Target writer format id such as `xyz`, `sdf`, `gjf`, or `cml`. | Unsupported formats raise the registry writer error. Writer-specific validation errors propagate. |
| `file_path` | frame, file | `None` | Output file path used only when `write_to_disk=True`. | Ignored when not writing. Directory paths assert. Memory-only objects need this path when writing. |
| `output_dir` | batch | `None` | Output directory used only when `write_to_disk=True`. | Ignored when not writing. Python API expects an existing directory when writing; CLI creates it. |
| `frame` | file, batch | `-1` | Frame selector: `int | Sequence[int] | "all"`. | Selector type errors are raised before writer dispatch. Out-of-range handling is writer-dependent unless a writer validates. |
| `embed_in_one_file` | file, batch | `True` | Combine selected frames into one rendered output when the format supports it. | Writer errors propagate. |
| `write_to_disk` | frame, file, batch | `False` | Write rendered content to disk. Without an explicit path, disk-backed objects write beside the source file. | Write errors propagate. Memory-only objects without `file_path` raise `ValueError`. |
| `n_jobs` | batch | `1` | Parallel jobs for batch conversion. | Joblib and worker exceptions propagate through `parallel_execute`. |
| `graph_policy` | frame, file, batch via `**kwargs` | Writer default | Molecular graph policy passed to the registry. | Invalid policies or missing graph data follow registry/writer errors. |
| `**kwargs` | frame, file, batch | None | Writer-specific options, for example Gaussian input options. | Unknown or invalid options follow the selected writer implementation. |

Returns:

| Surface | Return value |
| --- | --- |
| `FrameFormatTransformMixin.format_transform(...)` | `str` rendered for one frame. |
| `FormatTransformMixin.format_transform(...)` | `str` when embedded, otherwise `list[str]`. |
| `BatchFormatTransformMixin.format_transform(...)` | `dict[str, str | list[str]]` mapping source path to rendered content. |

## `to_summary_df`

| Parameter | Default | Meaning | Error strategy |
| --- | --- | --- | --- |
| `mode` | `"frame"` | `"frame"` summarizes selected frames; `"file"` summarizes each file. | Any other value raises `ValueError`. |
| `frame` | `-1` | Frame selector: `int | Sequence[int] | "all"`. Used only in frame mode. | Selector type errors are raised before summary generation. |
| `n_jobs` | `1` | Parallel jobs for summary extraction. | Joblib and worker exceptions propagate. |
| `brief` | `True` | Forwarded to `to_summary_series`; `False` requests expanded fields. | Field-specific errors propagate from the summary implementation. |
| `flatten_columns` | `False` | Convert three-level MultiIndex columns to dot-separated strings such as `General.FrameID` and `Energy.total_energy.hartree`; empty unit levels are skipped. | No effect for non-MultiIndex columns. |
| `on_missing_frame` | `"skip"` | `"skip"` drops missing frame indices; `"error"` raises. | Invalid policy raises `ValueError`; missing frames with `"error"` raise `IndexError`. |
| `**kwargs` | None | Extra summary-series options. | Forwarded errors propagate. |

Return value: `pandas.DataFrame`. If no selected series exist, returns an empty
DataFrame.

## `parallel_execute`

| Parameter | Default | Meaning | Error strategy |
| --- | --- | --- | --- |
| `func` | Required | Callable invoked as `func(diskfile, *args, **kwargs)`. | Exceptions from `func` propagate. |
| `desc` | `""` | Progress-bar description. | Progress backend errors propagate. |
| `n_jobs` | `1` | Number of parallel jobs. | Joblib errors propagate. |
| `return_as` | `"list"` | Joblib return mode: `"list"`, `"generator"`, or `"generator_unordered"`. | Invalid values follow joblib errors. |
| `*args` | None | Positional arguments passed after each disk file. | Callable errors propagate. |
| `_diskfiles_snapshot` | `None` | Internal precomputed batch snapshot used to keep filtering/grouping aligned with generated results. | Intended for internal callers; wrong snapshots can produce caller-level alignment errors. |
| `return_results` | `None` | `True` always returns results; `False` exhausts execution and returns `None`; `None` collapses list results containing only `None` to `None`. | Generator modes are not auto-inspected unless explicitly exhausted with `False`. |
| `**kwargs` | None | Keyword arguments passed to `func`. | Callable errors propagate. |

Return value: `Iterable[R] | None`.

## CLI Chain

The CLI exposes `molop parse PATTERN [parse options] OPERATION ...`.
Operations returning a `FileBatchModelDisk` can be followed by more operations.
Terminal operations must be last and are validated before parsing files.

### Parse Options

| Parameter | Default | Return/effect | Error strategy |
| --- | --- | --- | --- |
| `PATTERN` | Required | Builds the initial `FileBatchModelDisk`. | Missing or unmatched inputs follow parser errors. |
| `--parser-detection` | `"auto"` | Selects parser detection mode. | Unknown parser ids follow parser/registry errors. |
| `--n-jobs`, `-j` | `-1` | Default parallelism for operations that do not override `--n-jobs`. | Joblib errors propagate. |
| `--output-format` | `"text"` | Controls terminal result rendering. | Invalid choices are rejected by Click. |

### `format-transform`

| Parameter | Default | Return/effect | Error strategy |
| --- | --- | --- | --- |
| `--format` | Required | Target writer format id. Terminal operation. | Missing value or unsupported writer raises CLI/registry errors. |
| `--output-dir` | `None` | Implies writing for CLI compatibility. | Directory is created before execution when writing. |
| `--frame` | `"-1"` | CLI frame selector forwarded as `frame`. | Invalid selector text raises CLI usage error. |
| `--embed / --no-embed` | `--embed` | Controls multi-frame embedding. | Writer errors propagate. |
| `--write / --no-write` | Auto | `--write` writes beside sources without `--output-dir`; `--no-write` forces render-only behavior. | Write errors propagate. |
| `--n-jobs` | Parse-level `--n-jobs` | Overrides batch transform parallelism. | Joblib errors propagate. |
| Dynamic writer options | None | Forwarded to the selected writer. | Dynamic parse errors raise CLI usage errors; writer errors propagate. |

When the operation writes files, stdout is suppressed. Render-only conversion
prints the returned mapping.

### `to-summary-df`

| Parameter | Default | Return/effect | Error strategy |
| --- | --- | --- | --- |
| `--mode` | `"frame"` | Summary mode. Terminal operation. | Invalid choices are rejected by Click. |
| `--frame` | `"-1"` | CLI frame selector forwarded as `frame`. | Invalid selector text raises CLI usage error. |
| `--n-jobs` | Parse-level `--n-jobs` | Overrides summary parallelism. | Joblib errors propagate. |
| `--brief / --full` | `--brief` | Controls compact versus expanded summary fields. | Summary errors propagate. |
| `--flatten-columns / --multi-index-columns` | `--flatten-columns` | Controls CSV/JSON-friendly flat columns. | No effect on non-MultiIndex columns. |
| `--on-missing-frame` | `"skip"` | Missing frame policy. | `"error"` raises `IndexError` when selected frames are absent. |
| `--out` | `None` | Writes summary to a file instead of stdout. | Parent directory is created; write errors propagate. |
| `--format` | `"csv"` | Output serialization for `--out` or stdout. | Invalid choices are rejected by Click. |
