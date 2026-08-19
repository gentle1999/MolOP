import importlib
import multiprocessing
import sys
import threading
import warnings
from collections.abc import Callable, Generator, Iterable
from contextlib import contextmanager
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
def AdaptiveProgress(iterable: Iterable[T], *args: Any, **kwargs: Any) -> Iterable[T]: ...
@overload
def AdaptiveProgress(iterable: None = ..., *args: Any, **kwargs: Any) -> TqdmObj: ...
def AdaptiveProgress(
    iterable: Iterable[T] | None = None, *args: Any, **kwargs: Any
) -> Iterable[T] | TqdmObj:
    global _best_tqdm_cls
    if _best_tqdm_cls is None:
        _best_tqdm_cls = _detect_best_backend()
    tqdm_factory = cast(Callable[..., Any], _best_tqdm_cls)
    with warnings.catch_warnings():
        if TqdmExperimentalWarning is not None:
            warnings.filterwarnings("ignore", category=TqdmExperimentalWarning)
        else:
            warnings.filterwarnings("ignore", message=".*rich is experimental/alpha.*")
        try:
            obj = tqdm_factory(iterable, *args, **kwargs)
        except ImportError:
            if _best_tqdm_cls is st_tqdm:
                raise
            _best_tqdm_cls = st_tqdm
            obj = st_tqdm(iterable, *args, **kwargs)
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

    iterator = AdaptiveProgress(
        iterable, desc=desc, total=total, disable=disable, **(tqdm_kwargs or {})
    )
    if n_jobs == 1:
        if effective_return_as in {"generator", "generator_unordered"}:
            return _finalize_parallel_results((func(item) for item in iterator), return_results)
        return _finalize_parallel_results([func(item) for item in iterator], return_results)

    parallel = Parallel(n_jobs=n_jobs, return_as=effective_return_as, **joblib_kwargs)
    enter = getattr(parallel, "__enter__", None)
    if effective_return_as in {"generator", "generator_unordered"}:
        # Do not leave loky's reusable executor alive while a later native
        # MolGR prewarm runs in this process.  The returned generator owns the
        # coordination guard until its results are fully consumed or abandoned.
        def process_generator() -> Iterable[R]:
            with loky_parallel_guard():
                if callable(enter):
                    with parallel as executor:
                        results = executor(delayed(func)(item) for item in iterator)
                        if results is None:
                            raise ValueError("The parallel map returned None.")
                        yield from cast(Iterable[R], results)
                else:
                    results = parallel(delayed(func)(item) for item in iterator)
                    if results is None:
                        raise ValueError("The parallel map returned None.")
                    yield from cast(Iterable[R], results)

        return _finalize_parallel_results(process_generator(), return_results)

    with loky_parallel_guard():
        if callable(enter):
            with parallel as executor:
                results = executor(delayed(func)(item) for item in iterator)
        else:
            results = parallel(delayed(func)(item) for item in iterator)
    if results is None:
        raise ValueError("The parallel map returned None.")
    return _finalize_parallel_results(cast(Iterable[R], results), return_results)
