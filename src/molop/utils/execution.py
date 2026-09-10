"""Execution boundaries for MolOP's Python and native workloads.

This module owns process identity checks, joblib/loky lifecycle tracking, and
the mutual exclusion boundary around MolGR reconstruction.  Progress rendering
imports these helpers for compatibility, but does not own their state.
"""

from __future__ import annotations

import multiprocessing
import threading
from collections.abc import Callable, Generator
from contextlib import contextmanager
from typing import Any, cast

from joblib.externals import loky
from molgr.process_guard import ensure_current_process


DEFAULT_JOBLIB_BACKEND = "loky"

_parallel_state = threading.Condition()
_active_loky_threads: dict[int, int] = {}
_owned_loky_executor_ids: set[int] = set()
_native_reconstruction_active = False
_native_reconstruction_owner: int | None = None


class NativeReconstructionConcurrencyError(RuntimeError):
    """Raised when native reconstruction would overlap unsafe parallel work."""


def is_loky_worker() -> bool:
    """Return whether the current process is a joblib child worker."""

    process = multiprocessing.current_process()
    return process.name != "MainProcess" and (
        bool(getattr(process, "_identity", ()))
        or process.name.startswith(("LokyProcess", "ForkPoolWorker", "SpawnPoolWorker"))
    )


def _reusable_loky_executor() -> Any | None:
    """Return loky's current reusable executor when joblib exposes one."""

    return getattr(loky.reusable_executor, "_executor", None)


def _remember_owned_loky_executor(previous: Any | None) -> None:
    """Remember an executor created by a MolOP-managed call."""

    current = _reusable_loky_executor()
    if current is not None and (previous is None or id(previous) in _owned_loky_executor_ids):
        _owned_loky_executor_ids.add(id(current))


def _shutdown_reusable_loky_executor() -> None:
    """Drain and discard an idle reusable executor owned by MolOP.

    A reusable executor created by another library is never closed here.  A
    pool with pending work is rejected because native reconstruction cannot
    safely overlap it, regardless of ownership.
    """

    reusable_executor = loky.reusable_executor
    if not hasattr(reusable_executor, "_executor_lock") or not hasattr(
        reusable_executor, "_executor"
    ):
        raise NativeReconstructionConcurrencyError(
            "Cannot verify joblib/loky executor state before MolGR native "
            "reconstruction; this joblib version is unsupported by the native guard."
        )
    executor_lock = cast(Any, reusable_executor._executor_lock)
    executor = cast(Any, reusable_executor._executor)
    if executor_lock is None:
        raise NativeReconstructionConcurrencyError(
            "Cannot verify joblib/loky executor state before MolGR native "
            "reconstruction; this joblib version is unsupported by the native guard."
        )
    with cast(Any, executor_lock):
        if executor is None:
            return
        missing = object()
        pending_work = getattr(executor, "_pending_work_items", missing)
        running_work = getattr(executor, "_running_work_items", missing)
        shutdown = getattr(executor, "shutdown", None)
        if pending_work is missing or running_work is missing or not callable(shutdown):
            raise NativeReconstructionConcurrencyError(
                "Cannot verify joblib/loky executor state before MolGR native "
                "reconstruction; this joblib version is unsupported by the native guard."
            )
        if pending_work or running_work:
            raise NativeReconstructionConcurrencyError(
                "MolGR topology reconstruction cannot start while an external "
                "joblib/loky task is still running; consume or close that "
                "parallel result before prewarming topologies."
            )
        if id(executor) not in _owned_loky_executor_ids:
            return
        shutdown(wait=True, kill_workers=False)
        if getattr(reusable_executor, "_executor", None) is executor:
            reusable_executor._executor = None
            reusable_executor._executor_kwargs = None
        _owned_loky_executor_ids.discard(id(executor))


def _reject_active_child_processes() -> None:
    """Reject unsafe unmanaged children at the native boundary."""

    try:
        children = multiprocessing.active_children()
    except Exception as exc:
        raise NativeReconstructionConcurrencyError(
            "Cannot verify child-process state before MolGR native reconstruction."
        ) from exc
    unsafe_children = [
        child
        for child in children
        if not getattr(child, "name", "").startswith("LokyProcess-")
        and getattr(child, "_start_method", None) not in {"spawn", "forkserver"}
    ]
    if unsafe_children:
        names = ", ".join(
            f"{getattr(child, 'name', type(child).__name__)}" for child in unsafe_children[:4]
        )
        suffix = "..." if len(unsafe_children) > 4 else ""
        raise NativeReconstructionConcurrencyError(
            "MolGR topology reconstruction cannot start while unmanaged child "
            f"processes are active ({names}{suffix}); join or close them first."
        )


@contextmanager
def loky_parallel_guard() -> Generator[None, None, None]:
    """Serialize MolOP loky work against native reconstruction."""

    thread_id = threading.get_ident()
    executor_before = _reusable_loky_executor()
    with _parallel_state:
        if _native_reconstruction_active:
            if _native_reconstruction_owner == thread_id:
                message = (
                    "joblib parallel execution cannot start from inside the MolGR native "
                    "reconstruction guard."
                )
            else:
                message = (
                    "joblib parallel execution cannot start while the MolGR native "
                    "reconstruction guard is active; prewarm before starting parallel work."
                )
            raise NativeReconstructionConcurrencyError(message)
        _active_loky_threads[thread_id] = _active_loky_threads.get(thread_id, 0) + 1
    try:
        yield
    finally:
        with _parallel_state:
            remaining = _active_loky_threads.get(thread_id, 1) - 1
            if remaining > 0:
                _active_loky_threads[thread_id] = remaining
            else:
                _active_loky_threads.pop(thread_id, None)
            _parallel_state.notify_all()
        _remember_owned_loky_executor(executor_before)


@contextmanager
def native_reconstruction_guard(
    *,
    process_validator: Callable[[str], None] | None = None,
) -> Generator[None, None, None]:
    """Acquire the process-wide native reconstruction boundary."""

    validator = process_validator or ensure_current_process
    try:
        validator("MolOP native reconstruction")
    except RuntimeError as exc:
        raise NativeReconstructionConcurrencyError(str(exc)) from exc

    global _native_reconstruction_active, _native_reconstruction_owner
    thread_id = threading.get_ident()
    with _parallel_state:
        if _active_loky_threads:
            raise NativeReconstructionConcurrencyError(
                "MolGR topology reconstruction cannot start while a MolOP-managed "
                "joblib task or result generator is active. Prewarm topologies before "
                "starting parallel execution."
            )
        if _native_reconstruction_active:
            if _native_reconstruction_owner == thread_id:
                message = "MolGR topology reconstruction guard cannot be re-entered."
            else:
                message = (
                    "MolGR topology reconstruction guard is already active in another "
                    "call; concurrent native reconstruction is not supported."
                )
            raise NativeReconstructionConcurrencyError(message)
        _native_reconstruction_active = True
        _native_reconstruction_owner = thread_id
    try:
        _shutdown_reusable_loky_executor()
        _reject_active_child_processes()
        yield
    finally:
        with _parallel_state:
            _native_reconstruction_active = False
            _native_reconstruction_owner = None
            _parallel_state.notify_all()


__all__ = [
    "DEFAULT_JOBLIB_BACKEND",
    "NativeReconstructionConcurrencyError",
    "is_loky_worker",
    "loky_parallel_guard",
    "native_reconstruction_guard",
]
