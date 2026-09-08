"""Runtime defaults for native numerical libraries."""

from __future__ import annotations

import os


_NATIVE_THREAD_ENV_VARS = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
)


def configure_native_thread_limits() -> None:
    """Default native numerical libraries to one thread per process.

    MolOP uses process-level parallelism for batch work.  Limiting nested
    native thread pools avoids oversubscription and reduces instability around
    RDKit/Open Babel native calls.  Explicit user settings remain untouched.
    """
    for variable in _NATIVE_THREAD_ENV_VARS:
        os.environ.setdefault(variable, "1")
