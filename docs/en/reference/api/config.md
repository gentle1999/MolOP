# molop.config

This page provides the API reference for the `molop.config` module.

The module exposes the process-wide `molopconfig` object. It controls progress
display, the maximum parallel-job limit, logging, and structure-recovery
options. Prefer explicit function arguments for per-call behavior; change the
global configuration only when the policy should apply to subsequent calls.

## Parallelism defaults

`max_jobs` defaults to `None`, which means automatic process-aware detection. MolOP takes the most
restrictive value from joblib's scheduler/container quota, CPU affinity, and `os.cpu_count()`. Thus
`n_jobs=-1` uses all CPUs available to the current process, without exceeding a scheduler allocation.

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

::: molop.config
