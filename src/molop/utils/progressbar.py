"""Adaptive progress bars and parallel execution helpers.

Progress rendering is delegated to tqdm (with an auto-selected backend), but
every bar also broadcasts machine-readable :class:`ProgressEvent` notifications.
Third-party tools can capture MolOP progress without parsing terminal output:

- **In-process** — register a listener with
  :func:`register_progress_listener` (or the :func:`progress_listener` context
  manager) and receive every event synchronously, or use the thread-safe
  :class:`ProgressRecorder` to keep a live snapshot of every running bar.
- **Cross-process** — use :func:`progress_jsonl_sink` to append JSON Lines to a
  file that other tools can tail or read.
"""

import importlib
import itertools
import json
import logging
import multiprocessing
import os
import sys
import threading
import time
import warnings
from collections.abc import Callable, Generator, Iterable
from contextlib import contextmanager, suppress
from dataclasses import asdict, dataclass
from pathlib import Path
from types import TracebackType
from typing import TYPE_CHECKING, Any, Literal, NoReturn, TypeAlias, TypeVar, cast, overload

from joblib import Parallel, delayed
from joblib.parallel import get_active_backend
from molgr.process_guard import ensure_current_process
from tqdm import tqdm as st_tqdm


try:
    from tqdm.std import TqdmExperimentalWarning
except ImportError:
    TqdmExperimentalWarning = None


T = TypeVar("T")
R = TypeVar("R")

# ``loky`` starts independent interpreter processes (spawn-like on POSIX) and
# does not inherit MolGR's native runtime state through ``fork``. Keep this
# explicit at every MolOP-managed process boundary instead of mutating
# joblib's process-global context.
DEFAULT_JOBLIB_BACKEND = "loky"


_parallel_state = threading.Condition()
_active_loky_threads: dict[int, int] = {}
_native_reconstruction_active = False
_native_reconstruction_owner: int | None = None


class NativeReconstructionConcurrencyError(RuntimeError):
    """Raised when native reconstruction would overlap Python parallel work."""


def is_loky_worker() -> bool:
    """Return whether the current process is a joblib child worker.

    ``joblib`` names loky workers ``LokyProcess-*`` and its legacy
    multiprocessing backend ``ForkPoolWorker-*``/``SpawnPoolWorker-*``.
    Checking process identity instead of one naming convention keeps the
    native boundary intact for either process backend. Threading backends stay
    in the main process and are intentionally not classified as workers.
    """

    process = multiprocessing.current_process()
    return process.name != "MainProcess" and (
        bool(getattr(process, "_identity", ()))
        or process.name.startswith(("LokyProcess", "ForkPoolWorker", "SpawnPoolWorker"))
    )


def _shutdown_reusable_loky_executor() -> None:
    """Drain and discard joblib's process-global reusable executor.

    joblib's ``Parallel`` context only releases per-call resources; loky keeps
    its reusable worker pool alive by design.  MolGR native reconstruction must
    not start while those workers and queue threads are still alive, so the
    native boundary explicitly drains that pool.  This is isolated here because
    joblib does not expose a public "shutdown reusable executor" API.
    """

    try:
        from joblib.externals.loky import reusable_executor
    except ImportError:
        return

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
        # No MolOP-managed parallel operation can be active here.  Waiting for
        # an idle external loky pool is safe and removes all worker/queue
        # threads before native reconstruction starts.
        shutdown(wait=True, kill_workers=False)
        if getattr(reusable_executor, "_executor", None) is executor:
            reusable_executor._executor = None
            reusable_executor._executor_kwargs = None


def _reject_active_child_processes() -> None:
    """Reject unmanaged Python child processes at the native boundary."""

    try:
        children = multiprocessing.active_children()
    except Exception as exc:
        raise NativeReconstructionConcurrencyError(
            "Cannot verify child-process state before MolGR native reconstruction."
        ) from exc
    if children:
        names = ", ".join(
            f"{getattr(child, 'name', type(child).__name__)}" for child in children[:4]
        )
        suffix = "..." if len(children) > 4 else ""
        raise NativeReconstructionConcurrencyError(
            "MolGR topology reconstruction cannot start while unmanaged child "
            f"processes are active ({names}{suffix}); join or close them first."
        )


@contextmanager
def loky_parallel_guard() -> Generator[None, None, None]:
    """Serialize MolOP loky work against the native reconstruction boundary."""

    thread_id = threading.get_ident()
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


@contextmanager
def native_reconstruction_guard() -> Generator[None, None, None]:
    """Acquire the process-wide native reconstruction boundary."""

    # MolGR records the PID that initialized its native/Open Babel runtime and
    # rejects POSIX fork children. This also allows a fresh spawn/loky worker to
    # perform an isolated single-molecule reconstruction after importing MolGR
    # in that worker.
    try:
        ensure_current_process("MolOP native reconstruction")
    except RuntimeError as exc:
        # Keep MolGR's precise PID-baseline diagnostic while making it visible
        # to callers that deliberately fail closed on native-boundary errors.
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


if TYPE_CHECKING:
    from tqdm.notebook import tqdm_notebook
    from tqdm.rich import tqdm_rich

    TqdmType: TypeAlias = (
        type[st_tqdm[NoReturn]] | type[tqdm_notebook[NoReturn]] | type[tqdm_rich[NoReturn]]
    )
    TqdmObj: TypeAlias = st_tqdm[NoReturn] | tqdm_notebook[NoReturn] | tqdm_rich[NoReturn]
else:
    TqdmType: TypeAlias = Any
    TqdmObj: TypeAlias = Any


_progress_logger = logging.getLogger("molop.utils.progressbar")


@dataclass(frozen=True, slots=True)
class ProgressEvent:
    """A machine-readable progress notification emitted by MolOP progress bars.

    Attributes:
        task_id: Stable identifier of the bar instance that emitted the event.
        event: Kind of event: ``"start"``, ``"update"``, ``"description"`` or
            ``"close"``. Events are emitted regardless of whether the bar is
            rendered (``disable=True`` only suppresses terminal output).
        description: Current bar description (the ``desc`` argument), possibly empty.
        done: Completed units so far.
        total: Total units, or None while unknown (e.g. streamed inputs).
        timestamp: Wall-clock Unix timestamp of the event.
    """

    task_id: int
    event: Literal["start", "update", "description", "close"]
    description: str
    done: int | float
    total: int | float | None
    timestamp: float


ProgressListener: TypeAlias = Callable[[ProgressEvent], None]


class _ProgressListenerRegistration:
    __slots__ = ("listener",)

    def __init__(self, listener: ProgressListener) -> None:
        self.listener = listener


_progress_listeners: list[_ProgressListenerRegistration] = []
_progress_listeners_lock = threading.Lock()
_progress_task_ids = itertools.count(1)


def _same_listener(left: ProgressListener, right: ProgressListener) -> bool:
    if left is right:
        return True
    try:
        comparison = left == right
        return bool(comparison)
    except Exception:
        return False


def register_progress_listener(listener: ProgressListener) -> Callable[[], None]:
    """Register ``listener`` to receive every MolOP progress event.

    Listeners are invoked synchronously from the thread driving the progress
    bar, once per event, and must be cheap; exceptions raised by a listener are
    logged and never interrupt progress. Returns an unregister callable.
    """

    registration = _ProgressListenerRegistration(listener)
    with _progress_listeners_lock:
        _progress_listeners.append(registration)

    def unregister() -> None:
        with _progress_listeners_lock, suppress(ValueError):
            _progress_listeners.remove(registration)

    return unregister


def unregister_progress_listener(listener: ProgressListener) -> None:
    """Remove a previously registered progress listener (no-op if absent)."""

    with _progress_listeners_lock:
        _progress_listeners[:] = [
            registration
            for registration in _progress_listeners
            if not _same_listener(registration.listener, listener)
        ]


@contextmanager
def progress_listener(listener: ProgressListener) -> Generator[ProgressListener, None, None]:
    """Register ``listener`` for the duration of the ``with`` block."""

    unregister = register_progress_listener(listener)
    try:
        yield listener
    finally:
        unregister()


class ProgressRecorder:
    """Thread-safe, snapshot-based tracker of every live MolOP progress bar.

    The recorder is itself a progress listener. Register it with
    :func:`register_progress_listener` or simply use it as a context manager,
    and it maintains an up-to-date view of every running progress bar keyed by
    ``task_id``. Because events are state snapshots (each carries the full
    ``done``/``total``/``description``), the recorded state is always in sync
    with the bars — it can be read from any thread at any time, including
    while the bars are still running.

    Example
    -------
    .. code-block:: python

        from molop.utils.progressbar import ProgressRecorder, parallel_map

        rec = ProgressRecorder()
        with rec:  # registers and unregisters the listener automatically
            results = parallel_map(
                lambda value: value * 2,
                range(3),
                n_jobs=1,
                total=3,
                disable=True,
                return_results=True,
            )

        print(rec.history())  # full event timeline

        # Mid-run, from another thread (e.g. a GUI poller or a dashboard):
        #   for task_id, event in rec.snapshot().items():
        #       widget.set(event.done / event.total)
    """

    def __init__(self, *, keep_history: bool = True) -> None:
        self._lock = threading.Lock()
        self._active: dict[int, ProgressEvent] = {}
        self._history: list[ProgressEvent] = []
        self._keep_history = keep_history
        self._unregisters: list[Callable[[], None]] = []

    def __enter__(self) -> "ProgressRecorder":
        self._unregisters.append(register_progress_listener(self))
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_value: BaseException | None,
        traceback: TracebackType | None,
    ) -> None:
        if self._unregisters:
            self._unregisters.pop()()

    def __call__(self, event: ProgressEvent) -> None:
        with self._lock:
            if event.event == "close":
                self._active.pop(event.task_id, None)
            else:
                self._active[event.task_id] = event
            if self._keep_history:
                self._history.append(event)

    def snapshot(self) -> dict[int, ProgressEvent]:
        """Latest state of every live progress bar, keyed by ``task_id``."""

        with self._lock:
            return dict(self._active)

    def active_tasks(self) -> tuple[int, ...]:
        """Task ids of the currently live progress bars."""

        with self._lock:
            return tuple(self._active)

    def latest(self, task_id: int) -> ProgressEvent | None:
        """Latest event of ``task_id``, or None if that bar is not live."""

        with self._lock:
            return self._active.get(task_id)

    def percent(self, task_id: int) -> float | None:
        """Completion percentage (0-100) of ``task_id``, or None if unknown."""

        with self._lock:
            event = self._active.get(task_id)
        if event is None or not event.total:
            return None
        return event.done / event.total * 100.0

    def history(self) -> tuple[ProgressEvent, ...]:
        """Full recorded event timeline (empty when ``keep_history=False``)."""

        with self._lock:
            return tuple(self._history)

    def clear(self) -> None:
        """Drop all recorded state (active bars and history)."""

        with self._lock:
            self._active.clear()
            self._history.clear()


def _emit(event: ProgressEvent) -> None:
    with _progress_listeners_lock:
        listeners: list[ProgressListener] = []
        for registration in _progress_listeners:
            if not any(_same_listener(registration.listener, listener) for listener in listeners):
                listeners.append(registration.listener)
    for listener in listeners:
        try:
            listener(event)
        except Exception:
            _progress_logger.exception("Progress listener %r failed on %s", listener, event)


def _make_emitter(
    progress_callback: ProgressListener | None,
) -> Callable[[ProgressEvent], None]:
    if progress_callback is None:
        return _emit

    def emit(event: ProgressEvent) -> None:
        _emit(event)
        try:
            progress_callback(event)
        except Exception:
            _progress_logger.exception(
                "Progress callback %r failed on %s", progress_callback, event
            )

    return emit


class JsonLinesProgressListener:
    """Progress listener that appends each event as one JSON line to a file.

    Suitable for cross-process capture: an external tool can tail the file
    (``tail -f``, a reader process, a dashboard) and parse JSON Lines
    (https://jsonlines.org/). The file is flushed after every event.
    """

    def __init__(self, path: os.PathLike[str] | str) -> None:
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._handle = self.path.open("a", encoding="utf-8")
        self._lock = threading.Lock()
        self._closed = False

    def __call__(self, event: ProgressEvent) -> None:
        with self._lock:
            if self._closed:
                return
            self._handle.write(json.dumps(asdict(event)) + "\n")
            self._handle.flush()

    def close(self) -> None:
        with self._lock:
            if not self._closed:
                self._handle.close()
                self._closed = True


@contextmanager
def progress_jsonl_sink(
    path: os.PathLike[str] | str,
) -> Generator[JsonLinesProgressListener, None, None]:
    """Capture MolOP progress as JSON Lines appended to ``path``.

    Every progress event produced inside the ``with`` block is appended to
    ``path``, one JSON object per line, e.g.::

        {
            "task_id": 1,
            "event": "start",
            "description": "MolOP parsing with 4 processes",
            "done": 0,
            "total": 12,
            "timestamp": 1700000000.0,
        }

    The listener is unregistered and the file handle closed on block exit.
    """

    listener = JsonLinesProgressListener(path)
    unregister = register_progress_listener(listener)
    try:
        yield listener
    finally:
        unregister()
        listener.close()


_observed_cls_cache: dict[type, type] = {}


def _observed_cls_for(base: TqdmType) -> TqdmType:
    """Return a cached tqdm-compatible subclass that broadcasts progress events.

    Subclassing (instead of monkeypatching instances) keeps every tqdm
    attribute, method and rendering behavior intact while ``update``,
    ``set_description`` and ``close`` also emit :class:`ProgressEvent` objects.
    """

    cached = _observed_cls_cache.get(base)
    if cached is not None:
        return cached

    class _ObservedBar(base):  # type: ignore[valid-type, misc]
        """tqdm subclass that broadcasts :class:`ProgressEvent` notifications."""

        def __init__(
            self,
            iterable: Iterable[Any] | None = None,
            *args: Any,
            _molop_emitter: Callable[[ProgressEvent], None] | None = None,
            **kwargs: Any,
        ) -> None:
            initial_desc = kwargs.get("desc", args[0] if args else None)
            self._molop_task_id = next(_progress_task_ids)
            self._molop_emitter = _molop_emitter
            self._molop_closed = False
            self._molop_initialized = False
            self._molop_refresh_suppressed = 0
            self._molop_description = str(initial_desc or "")
            self._molop_description_rendered = False
            self._molop_last_state: tuple[str, int | float, int | float | None] | None = None
            super().__init__(iterable, *args, **kwargs)
            self._molop_initialized = True
            self._molop_emit("start")

        def _molop_state(self) -> tuple[str, int | float, int | float | None]:
            description = self._molop_description
            rendered_description = getattr(self, "desc", None)
            if rendered_description is not None:
                rendered_description = str(rendered_description)
                expected_description = description + (
                    ": " if self._molop_description_rendered else ""
                )
                if rendered_description != expected_description:
                    # A caller may mutate ``bar.desc`` directly. Preserve the
                    # same raw-description contract as the setter methods.
                    description = (
                        rendered_description[:-2]
                        if rendered_description.endswith(": ")
                        else rendered_description
                    )
                    self._molop_description = description
                    self._molop_description_rendered = rendered_description.endswith(": ")
            return (
                description,
                cast(int | float, getattr(self, "n", 0)),
                cast(int | float | None, getattr(self, "total", None)),
            )

        def _molop_emit(self, event: Literal["start", "update", "description", "close"]) -> None:
            emitter = getattr(self, "_molop_emitter", None)
            if emitter is None:
                return
            timestamp_fn = getattr(time, "time", None)
            if not callable(timestamp_fn):
                return
            timestamp = cast(Callable[[], float], timestamp_fn)()
            description, done, total = self._molop_state()
            self._molop_last_state = (description, done, total)
            try:
                emitter(
                    ProgressEvent(
                        task_id=self._molop_task_id,
                        event=event,
                        description=description,
                        done=done,
                        total=total,
                        timestamp=timestamp,
                    )
                )
            except Exception:
                # tqdm may invoke __del__ while module globals are being torn
                # down. Progress cleanup must not produce an ignored exception.
                return

        def _molop_emit_if_changed(self, event: Literal["update", "description", "close"]) -> None:
            if getattr(self, "_molop_closed", False):
                return
            state = self._molop_state()
            previous = self._molop_last_state
            if state == previous:
                return
            if event == "update" and previous is not None and state[0] != previous[0]:
                event = "description"
            self._molop_emit(event)

        def refresh(self, *args: Any, **kwargs: Any) -> Any:
            result = super().refresh(*args, **kwargs)
            if getattr(self, "_molop_initialized", False) and not getattr(
                self, "_molop_refresh_suppressed", 0
            ):
                self._molop_emit_if_changed("update")
            return result

        def update(self, n: int | float = 1) -> bool | None:
            if getattr(self, "_molop_closed", False):
                return None
            if self.disable:
                # tqdm does not track state for disabled bars; keep the
                # machine-readable counter accurate for external listeners.
                self.n += n
                self._molop_emit("update")
                return None
            self._molop_refresh_suppressed += 1
            try:
                result = super().update(n)
            finally:
                self._molop_refresh_suppressed -= 1
            self._molop_emit("update")
            return result

        def set_description(self, desc: str | None = None, refresh: bool = True) -> None:
            self._molop_refresh_suppressed += 1
            try:
                super().set_description(desc=desc, refresh=refresh)
            finally:
                self._molop_refresh_suppressed -= 1
            self._molop_description = str(desc or "")
            self._molop_description_rendered = True
            if getattr(self, "_molop_initialized", False) and not getattr(
                self, "_molop_closed", False
            ):
                self._molop_emit("description")

        def set_description_str(self, desc: str | None = None, refresh: bool = True) -> None:
            self._molop_refresh_suppressed += 1
            try:
                super().set_description_str(desc=desc, refresh=refresh)
            finally:
                self._molop_refresh_suppressed -= 1
            self._molop_description = str(desc or "")
            self._molop_description_rendered = False
            if getattr(self, "_molop_initialized", False) and not getattr(
                self, "_molop_closed", False
            ):
                self._molop_emit("description")

        def reset(self, total: int | float | None = None) -> Any:
            if getattr(self, "_molop_closed", False):
                return None
            # tqdm's optional rich backend currently calls ``Progress.reset``
            # without the required task id. Update rich's task explicitly and
            # then use the std tqdm state reset for every backend.
            rich_progress = getattr(self, "_prog", None)
            rich_task_id = getattr(self, "_task_id", None)
            if rich_progress is not None and rich_task_id is not None:
                with suppress(Exception):
                    rich_progress.reset(rich_task_id, total=total)
            result = st_tqdm.reset(self, total=total)
            if getattr(self, "_molop_initialized", False):
                self._molop_emit_if_changed("update")
            return result

        def close(self) -> None:
            if self._molop_closed:
                return
            self._molop_closed = True
            self._molop_refresh_suppressed += 1
            try:
                try:
                    super().close()
                finally:
                    if getattr(self, "_molop_initialized", False):
                        self._molop_emit("close")
            finally:
                self._molop_refresh_suppressed -= 1

        def __iter__(self) -> Generator[Any, None, None]:
            iterable = self.iterable
            if iterable is None:
                raise TypeError("Iterating a progress bar requires an iterable.")
            if getattr(self, "_molop_closed", False):
                yield from iterable
                return
            try:
                if self.disable:
                    for obj in iterable:
                        yield obj
                        self.n += 1
                        self._molop_emit("update")
                else:
                    for obj in iterable:
                        yield obj
                        # tqdm's own __iter__ batches updates for speed, which
                        # would hide intermediate progress from listeners.
                        self.update(1)
            finally:
                self.close()

    _ObservedBar.__name__ = f"_MolOPObserved{base.__name__}"
    _ObservedBar.__qualname__ = _ObservedBar.__name__
    _observed_cls_cache[base] = _ObservedBar
    return _ObservedBar


_best_tqdm_cls: TqdmType | None = None


def _is_notebook() -> bool:
    try:
        ipython_module = sys.modules.get("IPython")
        if ipython_module is None:
            return False
        get_ipython = getattr(ipython_module, "get_ipython", None)
        if get_ipython is None:
            return False
        ipython_instance = get_ipython()
        if ipython_instance is None:
            return False
        if "IPKernelApp" not in ipython_instance.config:
            return False
    except ImportError:
        return False
    except AttributeError:
        return False
    return True


def _detect_best_backend() -> TqdmType:
    if _is_notebook():
        try:
            notebook_module = importlib.import_module("tqdm.notebook")
            nb_tqdm = notebook_module.tqdm_notebook

            if getattr(notebook_module, "IProgress", None) is not None:
                return nb_tqdm
        except ImportError:
            pass
    else:
        try:
            from tqdm.rich import tqdm_rich as rich_tqdm

            return rich_tqdm
        except ImportError:
            pass
    return st_tqdm


@overload
def AdaptiveProgress(
    iterable: Iterable[T],
    *args: Any,
    progress_callback: ProgressListener | None = None,
    **kwargs: Any,
) -> Iterable[T]: ...
@overload
def AdaptiveProgress(
    iterable: None = ...,
    *args: Any,
    progress_callback: ProgressListener | None = None,
    **kwargs: Any,
) -> TqdmObj: ...
def AdaptiveProgress(
    iterable: Iterable[T] | None = None,
    *args: Any,
    progress_callback: ProgressListener | None = None,
    **kwargs: Any,
) -> Iterable[T] | TqdmObj:
    global _best_tqdm_cls
    if _best_tqdm_cls is None:
        _best_tqdm_cls = _detect_best_backend()
    emitter = _make_emitter(progress_callback)
    tqdm_factory = cast(Callable[..., Any], _observed_cls_for(_best_tqdm_cls))
    with warnings.catch_warnings():
        if TqdmExperimentalWarning is not None:
            warnings.filterwarnings("ignore", category=TqdmExperimentalWarning)
        else:
            warnings.filterwarnings("ignore", message=".*rich is experimental/alpha.*")
        try:
            obj = tqdm_factory(iterable, *args, _molop_emitter=emitter, **kwargs)
        except ImportError:
            if _best_tqdm_cls is st_tqdm:
                raise
            _best_tqdm_cls = st_tqdm
            tqdm_factory = cast(Callable[..., Any], _observed_cls_for(st_tqdm))
            obj = tqdm_factory(iterable, *args, _molop_emitter=emitter, **kwargs)
    if iterable is None:
        return cast(TqdmObj, obj)
    return cast(Iterable[T], obj)


def _finalize_parallel_results(
    results: Iterable[R], return_results: bool | None
) -> Iterable[R] | None:
    if return_results is False:
        for _ in results:
            pass
        return None
    if (
        return_results is None
        and isinstance(results, list)
        and results
        and all(result is None for result in results)
    ):
        return None
    return results


def _close_progress_bar(progress: Any) -> None:
    close = getattr(progress, "close", None)
    if callable(close):
        with suppress(Exception):
            close()


def _install_parallel_progress_hook(
    parallel: Any,
    progress: Any,
) -> tuple[Callable[[], None], bool]:
    """Mirror joblib's completed-task counter into an adaptive progress bar."""

    completed_seen = 0
    print_progress = getattr(parallel, "print_progress", None)
    if not callable(print_progress):
        return (lambda: None), False

    def observe() -> None:
        nonlocal completed_seen
        completed = getattr(parallel, "n_completed_tasks", None)
        if not isinstance(completed, int):
            return
        delta = completed - completed_seen
        if delta <= 0:
            return
        completed_seen = completed
        with suppress(Exception):
            progress.update(delta)

    def hooked_print_progress(*args: Any, **kwargs: Any) -> Any:
        result = print_progress(*args, **kwargs)
        observe()
        return result

    with suppress(Exception):
        parallel.print_progress = hooked_print_progress
        return observe, True
    return (lambda: None), False


def _with_default_joblib_backend(joblib_kwargs: dict[str, Any]) -> dict[str, Any]:
    """Select MolOP's safe process backend unless the caller chose another mode."""

    if joblib_kwargs.get("backend") is not None:
        return joblib_kwargs
    joblib_kwargs.pop("backend", None)
    if joblib_kwargs.get("prefer") == "threads" or joblib_kwargs.get("require") == "sharedmem":
        return joblib_kwargs

    # Preserve threading and third-party backends. Normalize only joblib's
    # built-in process choices: the default is made explicit, while the legacy
    # MultiprocessingBackend could otherwise inherit POSIX fork semantics.
    active_backend, _ = get_active_backend()
    if active_backend.__class__.__name__ in {"LokyBackend", "MultiprocessingBackend"}:
        joblib_kwargs["backend"] = DEFAULT_JOBLIB_BACKEND
    return joblib_kwargs


@overload
def parallel_map(
    func: Callable[[T], R],
    iterable: Iterable[T],
    n_jobs: int = -1,
    total: int | None = None,
    desc: str = "Processing",
    disable: bool = False,
    return_as: str | None = None,
    tqdm_kwargs: dict[str, Any] | None = None,
    progress_callback: ProgressListener | None = None,
    *,
    return_results: Literal[False],
    **joblib_kwargs: Any,
) -> None: ...
@overload
def parallel_map(
    func: Callable[[T], R],
    iterable: Iterable[T],
    n_jobs: int = -1,
    total: int | None = None,
    desc: str = "Processing",
    disable: bool = False,
    return_as: str | None = None,
    tqdm_kwargs: dict[str, Any] | None = None,
    progress_callback: ProgressListener | None = None,
    *,
    return_results: Literal[True],
    **joblib_kwargs: Any,
) -> Iterable[R]: ...
@overload
def parallel_map(
    func: Callable[[T], R],
    iterable: Iterable[T],
    n_jobs: int = -1,
    total: int | None = None,
    desc: str = "Processing",
    disable: bool = False,
    return_as: str | None = None,
    tqdm_kwargs: dict[str, Any] | None = None,
    progress_callback: ProgressListener | None = None,
    *,
    return_results: None = None,
    **joblib_kwargs: Any,
) -> Iterable[R] | None: ...
def parallel_map(
    func: Callable[[T], R],
    iterable: Iterable[T],
    n_jobs: int = -1,
    total: int | None = None,
    desc: str = "Processing",
    disable: bool = False,
    return_as: str | None = None,
    tqdm_kwargs: dict[str, Any] | None = None,
    progress_callback: ProgressListener | None = None,
    *,
    return_results: bool | None = None,
    **joblib_kwargs: Any,
) -> Iterable[R] | None:
    """
    A parallel wrapper with type hints and progress bar

    Parameters
    ----------
    func : Callable[[T], R]
        The function to apply to each item in the iterable.
    iterable : Iterable[T]
        The iterable to apply the function to.
    n_jobs : int, optional
        The number of parallel jobs to run. Default is -1, which means using all available cores.
    total : Optional[int], optional
        The total number of items in the iterable. Default is None, which means inferring the total from the iterable.
    desc : str, optional
        The description to display in the progress bar. Default is "Processing".
    disable : bool, optional
        Whether to disable the progress bar. Default is False.
    tqdm_kwargs : Optional[dict[str, Any]], optional
        Additional keyword arguments to pass to tqdm. Default is None.
    progress_callback : Optional[Callable[[ProgressEvent], None]], optional
        Optional callback invoked with a :class:`ProgressEvent` for every
        progress update of the wrapped bar, in addition to the global listener
        registry (see :func:`register_progress_listener`). Default is None.
    return_results : Optional[bool], optional
        When False, exhaust execution and return None. When None, list results
        containing only implicit None values are returned as None.
    **joblib_kwargs : Any
        Additional keyword arguments to pass to joblib.Parallel.

    Returns
    -------
    Iterable[R] | None
        The iterable of results after applying the function to each item, or None
        for side-effect-only calls.
    """
    original_joblib_kwargs = joblib_kwargs
    joblib_kwargs = _with_default_joblib_backend(dict(joblib_kwargs))
    effective_return_as = return_as
    if effective_return_as is None:
        # joblib's legacy multiprocessing backend cannot expose a result
        # generator. Keep the streaming default for loky and other backends,
        # but select a compatible eager result for this explicit backend.
        effective_return_as = (
            ("list" if joblib_kwargs.get("backend") == "multiprocessing" else "generator")
            if original_joblib_kwargs
            else "list"
        )

    tqdm_kwargs = dict(tqdm_kwargs or {})
    # Accept progress_callback either as a first-class argument or as a
    # pass-through key in tqdm_kwargs without double-binding it.
    tqdm_callback = tqdm_kwargs.pop("progress_callback", None)
    if progress_callback is None:
        progress_callback = tqdm_callback
    progress_total = total
    if progress_total is None:
        try:
            progress_total = len(iterable)  # type: ignore[arg-type]
        except (AttributeError, TypeError):
            progress_total = None

    def new_progress() -> Any:
        return AdaptiveProgress(
            None,
            desc=desc,
            total=progress_total,
            disable=disable,
            progress_callback=progress_callback,
            **tqdm_kwargs,
        )

    if n_jobs == 1:
        if effective_return_as in {"generator", "generator_unordered"}:

            def serial_generator() -> Iterable[R]:
                progress = new_progress()
                try:
                    for item in iterable:
                        result = func(item)
                        progress.update(1)
                        yield result
                finally:
                    _close_progress_bar(progress)

            return _finalize_parallel_results(serial_generator(), return_results)
        progress = new_progress()
        try:
            source = iter(iterable)
        except Exception:
            _close_progress_bar(progress)
            raise
        try:
            results: list[R] = []
            for item in source:
                results.append(func(item))
                progress.update(1)
            return _finalize_parallel_results(results, return_results)
        finally:
            _close_progress_bar(progress)

    parallel = Parallel(n_jobs=n_jobs, return_as=effective_return_as, **joblib_kwargs)
    enter = getattr(parallel, "__enter__", None)
    if effective_return_as in {"generator", "generator_unordered"}:
        # Do not leave loky's reusable executor alive while a later native
        # MolGR prewarm runs in this process.  The returned generator owns the
        # coordination guard until its results are fully consumed or abandoned.
        def process_generator() -> Iterable[R]:
            progress = new_progress()
            observe_parallel, has_parallel_hook = _install_parallel_progress_hook(
                parallel, progress
            )
            try:
                with loky_parallel_guard():
                    if callable(enter):
                        with parallel as executor:
                            source = iter(iterable)
                            results = executor(delayed(func)(item) for item in source)
                            if results is None:
                                raise ValueError("The parallel map returned None.")
                            for result in cast(Iterable[R], results):
                                if not has_parallel_hook:
                                    progress.update(1)
                                yield result
                    else:
                        source = iter(iterable)
                        results = parallel(delayed(func)(item) for item in source)
                        if results is None:
                            raise ValueError("The parallel map returned None.")
                        for result in cast(Iterable[R], results):
                            if not has_parallel_hook:
                                progress.update(1)
                            yield result
            finally:
                observe_parallel()
                _close_progress_bar(progress)

        return _finalize_parallel_results(process_generator(), return_results)

    progress = new_progress()
    try:
        source = iter(iterable)
    except Exception:
        _close_progress_bar(progress)
        raise
    observe_parallel, has_parallel_hook = _install_parallel_progress_hook(parallel, progress)
    try:
        with loky_parallel_guard():
            raw_results: Any
            if callable(enter):
                with parallel as executor:
                    raw_results = executor(delayed(func)(item) for item in source)
            else:
                raw_results = parallel(delayed(func)(item) for item in source)
        if raw_results is None:
            raise ValueError("The parallel map returned None.")
        observe_parallel()
        completed_results: list[R]
        if isinstance(raw_results, list):
            completed_results = raw_results
        else:
            completed_results = list(cast(Iterable[R], raw_results))
        if not has_parallel_hook:
            for _ in completed_results:
                progress.update(1)
        return _finalize_parallel_results(completed_results, return_results)
    finally:
        observe_parallel()
        _close_progress_bar(progress)
