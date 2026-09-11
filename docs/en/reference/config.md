# molop.config

This page provides the API reference for the `molop.config` module.

The module exposes the process-wide `molopconfig` object. It controls progress
display, the maximum parallel-job limit, logging, and structure-recovery
options. Prefer explicit function arguments for per-call behavior; change the
global configuration only when the policy should apply to subsequent calls.

## Configuration overview

The `MolOPConfig` fields and defaults are listed below. `from molop import molopconfig`
returns a process-wide shared object; changing it affects subsequent calls. Creating an
independent `MolOPConfig(...)` applies its progress, file-logging, and native-logging
settings during construction.

| Option | Type | Default | Description |
| --- | --- | --- | --- |
| `show_progress_bar` | `bool` | `True` | Show batch progress bars; also controls MolOP's console log handler. |
| `max_jobs` | `int \| None` | `None` | Process-wide parallel-job ceiling; `None` uses the CPUs available to the process. |
| `graph_reconstruction_backend` | `"cpp" \| "python"` | `"cpp"` | Backend used for molecular-graph reconstruction. |
| `reconstruction_failure_policy` | `"raise" \| "return_suspicious"` | `"raise"` | Raise on reconstruction failure, or retain a fallback graph marked suspicious. |
| `prewarm_topologies` | `bool` | `False` | Prewarm graph-dependent topologies in the parent process; otherwise reconstruct lazily in the worker using the result. |
| `make_dative_bonds` | `bool` | `True` | Create dative bonds during graph reconstruction. |
| `make_stereochemistry` | `bool` | `True` | Assign stereochemistry during graph reconstruction. |
| `force_unit_transform` | `bool` | `False` | Force unit conversion. |
| `parallel_max_size` | `int` | `8 * 1024**2` | Size threshold used by parallel dispatch and joblib data transfer, in bytes (8 MiB). |
| `max_recursion_depth` | `int` | `3000` | Python recursion depth requested by `set_max_recursion_depth()`; import does not change the interpreter limit. |
| `log_to_file` | `bool` | `False` | Enable MolOP file logging. |
| `log_file_path` | `str` | `"molop.log"` | File-log path; effective when `log_to_file` is enabled or `enable_file_logging()` is called. |
| `suppress_rdkit_logs` | `bool` | `True` | Suppress RDKit native diagnostics. |
| `suppress_openbabel_logs` | `bool` | `True` | Suppress Open Babel native diagnostics. |
| `use_dof_effect_drawer` | `bool` | `True` | Prefer the `rdkit-dof` depth-of-field drawer; use the standard RDKit drawer when disabled. |

`effective_max_jobs` and `effective_molgr_max_jobs` are read-only limits calculated from
available CPU resources and the configured ceilings.

## Import side effects and explicit initialization

Importing `molop` initializes the global `molopconfig` and suppresses RDKit/Open Babel native diagnostics
by default. It does not create `molop.log` or change the host process's Python recursion limit. `rdkit-dof`
is loaded only when depth-aware drawing or an explicit related setting requires it. Native diagnostics belong
to the native libraries' own output and are not automatically copied to MolOP's file log.

To retain native diagnostics, change the flags and apply the policy to the current process:

```python
from molop import molopconfig

molopconfig.suppress_rdkit_logs = False
molopconfig.suppress_openbabel_logs = False
molopconfig.configure_native_logging()
```

Keep the other flag set to `True` when only one native library should remain quiet. New `MolOPConfig(...)`
instances apply both native-log flags during construction; call `configure_native_logging()` after changing
flags on an existing object.

To enable file logging at the same time as quiet native diagnostics, create a configuration explicitly:

```python
from molop.config import MolOPConfig

config = MolOPConfig(
    log_to_file=True,
    log_file_path="run.log",
    suppress_rdkit_logs=True,
    suppress_openbabel_logs=True,
)
```

`enable_file_logging()` and `disable_file_logging()` own the file-handler lifecycle. Setting
`log_to_file = True` on the global object does not create a handler automatically; call
`enable_file_logging()` and use `disable_file_logging()` to turn it off. `set_log_level()` accepts
`DEBUG`, `INFO`, `WARNING`, `ERROR`, or `CRITICAL`.

Use the methods below to switch the progress bar and MolOP's console handler together:

```python
molopconfig.quiet()    # disable the progress bar and MolOP console logs
molopconfig.verbose()  # restore the progress bar and MolOP console logs
```

`set_max_recursion_depth()` changes the current process only when explicitly called.

Parent-process topology prewarming is controlled by `prewarm_topologies` and defaults to `False`.
With the default, graph-dependent operations reconstruct lazily in the spawn-like `loky` worker that
uses the result, avoiding a separate parent-process warmup step. Set it to `True` when a workflow needs
the deterministic parent-side cache:

```python
from molop import molopconfig

molopconfig.prewarm_topologies = True
```

MolGR reconstruction failures are strict by default. Enable raw fallback retention when the
workflow must preserve invalid structures for review:

```python
molopconfig.reconstruction_failure_policy = "return_suspicious"
```

Retained graphs are marked with `topology_reconstruction_status == "suspicious_fallback"`; they
must not be treated as trusted chemistry.

## Parallelism defaults

`max_jobs` defaults to `None`, which means automatic process-aware detection. MolOP takes the most
restrictive observable value from joblib/loky's scheduler and container limits, `os.cpu_count()`, and
the current process affinity where the operating system exposes it. The direct psutil dependency
enables joblib's affinity fallback on Windows and other supported platforms. Thus
`n_jobs=-1` uses the CPUs available to the current process without exceeding a detected allocation.

Linux supports affinity and cgroup quota detection. Windows affinity is read through psutil and
joblib also applies its safe worker limit for Windows process pools. macOS does not expose a portable
per-process CPU-affinity API, so its automatic value uses the CPUs visible to the process and any
joblib limit. On every platform, use `max_jobs`, `MOLOP_MAX_JOBS`, or `LOKY_MAX_CPU_COUNT` when the
runtime allocation is not observable from the operating system.

Set a positive `max_jobs` to impose a process-wide ceiling:

```python
from molop import molopconfig

molopconfig.max_jobs = 32
```

For deployments that construct the default configuration from the environment, set
`MOLOP_MAX_JOBS=32`. An explicit `MolOPConfig(max_jobs=None)` keeps automatic detection even when the
environment variable is present. Positive `n_jobs` values still request a per-call limit, capped by
`effective_max_jobs`; `n_jobs=1` remains the serial debugging mode. CLI users can set the same
process-wide ceiling for one invocation with `molop --max-jobs 32 parse ...`.

Tasks that may enter MolGR use the separate `effective_molgr_max_jobs` budget. Automatic and explicit
values are capped at `floor(2 * available_cpu_count / 3)` (with a minimum of one serial worker), then
by `max_jobs`. This covers topology reconstruction, graph-dependent summaries, exports, transforms,
and user callbacks such as `filter_custom`. File parsing and other operations that do not invoke
MolGR keep the full `effective_max_jobs` limit.

::: molop.config
