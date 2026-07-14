# API Contracts

This page fixes the public contract for the transform, summary, parallel
execution, and CLI chain surfaces.

## `AutoParser`

The Python API accepts a single `str` or `os.PathLike`, or a one-dimensional
iterable containing paths and glob patterns. Every pattern is expanded, the
result is materialized and de-duplicated by normalized absolute path, and one
`FileBatchModelDisk` is returned in absolute-path order. Input grouping and
input order are not retained.

Unmatched patterns contribute no files. Empty iterables return an empty batch.
Missing literal paths are passed to the batch parser, which warns and skips
them. Nested iterables and non-path members raise `TypeError` with the member
index. All parse options apply to every expanded path.

```python
from pathlib import Path

batch = AutoParser(["inputs/a.log", "inputs/set_b/*.log", Path("c.log")])
```

## Calculation Database View

`ChemFile` and `ChemFileFrame` are the authoritative calculation models. MolOP
does not maintain a parallel `Parsed*` DTO tree or a second calculation parsing
entry point. Database ingestion uses a serialization view of the same objects:

```python
from molop import AutoParser

chem_file = AutoParser(
    "calculation.out",
    capture_source_evidence=True,
    release_file_content=True,
)[0]
file_payload = chem_file.to_unitless_dump_with_unit_keys(exclude_none=True)
frame_payloads = (
    frame.to_unitless_dump_with_unit_keys(exclude_none=True)
    for frame in chem_file
)
```

File and frame payloads are deliberately separate. A `ChemFile` dump contains
artifact, segment, and other file metadata; it does not embed `frames` or the
private `_frames_` collection. Callers consume `frame_payloads` incrementally,
attach their own file identity and ordering context, and validate each payload
with their own Pydantic receiving model. Such a model should select the fields
it needs and normally use `extra="ignore"` for format-specific additions. A
transport may wrap a file payload and a page of frame payloads, but that wrapper
is not a MolOP DTO or a requirement to build one giant JSON document.

With source evidence enabled, each frame payload has three indices with
different scopes. `frame_id` is the append-local index in the current
`ChemFile`, `segment_frame_index` is the original index inside its segment, and
`file_frame_index` is the stable locator-order ordinal across the complete
source artifact. The latter two retain their source values under
`only_last_frame=True`; source indices are omitted when evidence capture is
disabled. Persisted source ordering must use `file_frame_index`, not a
`frame_id` that can be reset by selection.

Every normal parser lifecycle adds the fixed
`schema_version="molop-calculation-export-v1"`, canonical `source_format`, and
file-scoped `parser_provenance`. Provenance records the concrete parser, MolOP,
MolGR, and RDKit versions plus snapshots of the effective parser, MolOP, and MolGR
configuration. `effective_config_sha256` is the strict canonical-JSON SHA-256
of that snapshot. Provenance is not copied into frame payloads.

`to_unitless_dump_with_unit_keys()` recursively converts quantities to their
magnitudes and writes the normalized unit into the corresponding key as
`field (unit)`. Nested MolOP models, arrays, sequences, and mappings are handled
by the common `BaseDataClassWithUnit` implementation. The result is a dump, not
a second model with its own validation, conversion, or lifecycle.

Unit labels are part of this public serialization view. They are assembled from
canonical Pint unit names with MolOP's own deterministic grammar rather than
Pint's display-oriented `str(unit)` output. Products are sorted and joined with
`*`, powers use `^`, and a multi-term denominator is parenthesized. For example,
the stable labels are `hartree`, `kilocalorie/mole`,
`bohr^2*unified_atomic_mass_unit`, and `calorie/(kelvin*mole)`.

Arrays use JSON-safe lists by default. Database ingestion can retain independent
copies of numeric ndarrays and write deterministic NPY or other binary sidecars
without constructing a giant JSON document:

```python
frame_payload = frame.to_unitless_dump_with_unit_keys(
    exclude_none=True,
    array_mode="ndarray",
)
```

In `ndarray` mode, non-numeric arrays still become lists; the returned numeric
arrays never alias the model's internal arrays. The output is intentionally not
JSON serializable until the caller replaces those arrays with references or
chooses an encoding. All other keyword arguments are passed to Pydantic
`model_dump()` before recursive conversion, so `include`, `exclude`,
`exclude_none`, aliases, and the other model-dump selection rules keep their
normal semantics.

The export path materializes lazy public topology before Pydantic takes its
snapshot. Consequently the first frame payload already contains `bonds`,
`formal_charges`, and `formal_num_radicals`, and repeated dumps of an unchanged
frame are equal. Bond endpoints and per-atom topology arrays use the same source
atom order as the exported `atoms` and `coords`. This preparation is limited to
this export method; ordinary `model_dump()` retains its existing lazy behavior
and cost.

Trusted reconstruction also exports a map-free `topology_v3000_molblock`,
`source_to_topology_atom_permutation`, and reconstruction-configuration
provenance. When element and coordinate checks prove that a MolGR backend
preserved order, the permutation is explicitly the identity. Failed or
ambiguous reconstruction does not guess a mapping; it reports
`parse_presence["topology"]="parse_failed"` with a stable diagnostic code.

Source identity, locations, semantics, and evidence are optional fields on the
existing models.
Evidence capture is disabled by default and can be enabled with
`capture_source_evidence=True` for database ingestion or audit workflows.
`source_encoding` selects the strict decoder used for both parsing and byte
offset calculation and defaults to `utf-8`.
Supported parsers populate those fields during their normal parsing path before
retained source text is released; they do not run a separate scientific
extractor afterward. The format-independent protocol lives in
`molop.io.base_models.source`: `DecodedSource` preserves exact decoded offsets,
`LocatedTextBlock` and `LocatedSourceSegment` describe parser-owned boundaries,
and `SourceSpan` stores byte/character/line ranges. Artifact fields live on
`BaseChemFile`, location fields live on `BaseChemFileFrame`, and
coordinate-source fields live on `BaseCoordsFrame`, so coordinate, input, and
calculation formats share the same lifecycle.

Each `SourceSegmentEvidence` also carries portable `protocol` and
`task_requests` values. `protocol` is the plain-mapping projection of common
`model_chemistry`, and `task_requests` is a list of plain-mapping common
`QMTaskRequest` values; no Gaussian- or ORCA-specific semantic model leaks into
the evidence contract. Because this evidence comes from segment metadata, a
segment with no parseable frame still retains its requested protocol. Formal
database ingestion must set `capture_source_evidence=True`; the resulting
half-open byte, character, and line ranges and their SHA-256 values are
admission evidence.

File segmentation has one source of truth. Every format-specific file parser
declares a canonical `format_id` and implements `_quick_check_file_format()`
and `_locate_segments()`; the locator returns half-open segment and frame ranges
over the original text. In a segment with frames, those ranges must cover every
non-whitespace character; pure-whitespace gaps and zero-frame segments remain
valid. Optional
artifact- and segment-scoped metadata use `_parse_artifact_metadata()` and
`_parse_segment_metadata()`. There is no `_split_file()` or `_parse_metadata()`
fallback, and parsers must not rebuild frame text and search backward for spans.

MolOP reports parsed facts and optional evidence. Storage systems remain
responsible for independently verifying artifact bytes, choosing persisted
array encodings, canonicalizing database identity, applying admission/QC
policy, and creating reaction or path records.

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
