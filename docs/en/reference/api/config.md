# molop.config

This page provides the API reference for the `molop.config` module.

The module exposes the process-wide `molopconfig` object. It controls progress
display, the maximum parallel-job limit, logging, and structure-recovery
options. Prefer explicit function arguments for per-call behavior; change the
global configuration only when the policy should apply to subsequent calls.

## Import side effects and explicit initialization

Importing `molop` does not create `molop.log`, change the host process's Python recursion limit, or
silence RDKit/Open Babel native diagnostics by default. `rdkit-dof` is also loaded only when depth-aware
drawing or an explicit related setting requires it. Create a configuration explicitly when a workflow
wants file logging or quiet native diagnostics:

```python
from molop.config import MolOPConfig

config = MolOPConfig(
    log_to_file=True,
    log_file_path="run.log",
    suppress_rdkit_logs=True,
    suppress_openbabel_logs=True,
)
config.configure_native_logging()
```

`enable_file_logging()` and `disable_file_logging()` own the file-handler lifecycle;
`set_max_recursion_depth()` changes the current process only when the caller explicitly requests it.

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
