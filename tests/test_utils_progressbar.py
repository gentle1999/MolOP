import multiprocessing
import sys
import threading
from concurrent.futures import ThreadPoolExecutor
from types import ModuleType, SimpleNamespace
from typing import Any

import pytest
from joblib import parallel_config
from tqdm import tqdm as st_tqdm

import molop.utils.progressbar as progressbar_module
from molop.utils.progressbar import (
    AdaptiveProgress,
    NativeReconstructionConcurrencyError,
    _is_notebook,
    loky_parallel_guard,
    native_reconstruction_guard,
    parallel_map,
)


def test_adaptive_progress_is_iterable_and_yields_input_items(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(progressbar_module, "_best_tqdm_cls", None, raising=False)
    items = ["a", "b", "c"]
    yielded = list(AdaptiveProgress(items, disable=True))
    assert yielded == items


def test_parallel_map_n_jobs_one_returns_expected_results(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(progressbar_module, "_best_tqdm_cls", None, raising=False)
    items = [1, 2, 3, 4]
    results = parallel_map(lambda value: value * value, items, n_jobs=1, disable=True)
    assert results == [1, 4, 9, 16]


def _worker_process_identity(_value: int) -> tuple[str, str]:
    process = multiprocessing.current_process()
    return process.name, type(process).__module__


def test_parallel_map_default_uses_loky_workers() -> None:
    results = parallel_map(
        _worker_process_identity,
        [1, 2],
        n_jobs=2,
        disable=True,
        return_results=True,
    )

    assert all(name.startswith("LokyProcess-") for name, _module in results)
    assert all(module == "joblib.externals.loky.backend.process" for _name, module in results)


def test_parallel_map_suppresses_implicit_none_results_by_default(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(progressbar_module, "_best_tqdm_cls", None, raising=False)
    items = [1, 2, 3]
    seen: list[int] = []

    results = parallel_map(lambda value: seen.append(value), items, n_jobs=1, disable=True)

    assert results is None
    assert seen == items


def test_parallel_map_can_preserve_explicit_none_results(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(progressbar_module, "_best_tqdm_cls", None, raising=False)
    items = [1, 2, 3]

    results = parallel_map(
        lambda _value: None,
        items,
        n_jobs=1,
        disable=True,
        return_results=True,
    )

    assert results == [None, None, None]


def test_parallel_map_return_results_false_exhausts_generator(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(progressbar_module, "_best_tqdm_cls", None, raising=False)
    items = [1, 2, 3]
    seen: list[int] = []

    results = parallel_map(
        lambda value: seen.append(value),
        items,
        n_jobs=1,
        disable=True,
        return_as="generator",
        return_results=False,
    )

    assert results is None
    assert seen == items


def test_parallel_map_n_jobs_one_supports_generator_return_as(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(progressbar_module, "_best_tqdm_cls", None, raising=False)
    items = [1, 2, 3, 4]

    results = parallel_map(
        lambda value: value * value,
        items,
        n_jobs=1,
        disable=True,
        return_as="generator",
    )

    assert list(results) == [1, 4, 9, 16]


def test_is_notebook_returns_false_when_ipython_module_missing(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(progressbar_module, "_best_tqdm_cls", None, raising=False)
    monkeypatch.delitem(sys.modules, "IPython", raising=False)
    assert _is_notebook() is False


def test_is_notebook_returns_false_when_get_ipython_returns_none(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(progressbar_module, "_best_tqdm_cls", None, raising=False)
    monkeypatch.setitem(
        sys.modules,
        "IPython",
        SimpleNamespace(get_ipython=lambda: None),
    )
    assert _is_notebook() is False


def test_is_notebook_returns_false_when_ipkernelapp_missing_from_config(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(progressbar_module, "_best_tqdm_cls", None, raising=False)
    ipython_instance = SimpleNamespace(config={})
    monkeypatch.setitem(
        sys.modules,
        "IPython",
        SimpleNamespace(get_ipython=lambda: ipython_instance),
    )
    assert _is_notebook() is False


def test_is_notebook_returns_true_when_ipkernelapp_present_in_config(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(progressbar_module, "_best_tqdm_cls", None, raising=False)
    ipython_instance = SimpleNamespace(config={"IPKernelApp": True})
    monkeypatch.setitem(
        sys.modules,
        "IPython",
        SimpleNamespace(get_ipython=lambda: ipython_instance),
    )
    assert _is_notebook() is True


def test_adaptive_progress_falls_back_when_notebook_widgets_are_unavailable(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    notebook_module = ModuleType("tqdm.notebook")
    notebook_module.IProgress = None  # type: ignore[attr-defined]
    notebook_module.tqdm_notebook = object  # type: ignore[attr-defined]
    monkeypatch.setitem(sys.modules, "tqdm.notebook", notebook_module)
    monkeypatch.setattr(progressbar_module, "_is_notebook", lambda: True)
    monkeypatch.setattr(progressbar_module, "_best_tqdm_cls", None)

    progress = AdaptiveProgress([1, 2], disable=True)

    assert list(progress) == [1, 2]
    assert progressbar_module._best_tqdm_cls is st_tqdm


def test_adaptive_progress_falls_back_when_notebook_backend_initialization_fails(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    class BrokenNotebookProgress:
        def __init__(self, *_args: Any, **_kwargs: Any) -> None:
            raise ImportError("IProgress not found")

    monkeypatch.setattr(progressbar_module, "_best_tqdm_cls", BrokenNotebookProgress)

    progress = AdaptiveProgress([1, 2], disable=True)

    assert list(progress) == [1, 2]
    assert progressbar_module._best_tqdm_cls is st_tqdm


def test_parallel_map_joblib_kwargs_defaults(monkeypatch: pytest.MonkeyPatch) -> None:
    captured_kwargs = {}

    class StubParallel:
        def __init__(self, **kwargs: Any) -> None:
            nonlocal captured_kwargs
            captured_kwargs = kwargs

        def __call__(self, iterable: Any) -> list[Any]:
            return [None] * len(list(iterable))

    monkeypatch.setattr(progressbar_module, "Parallel", StubParallel)
    monkeypatch.setattr(progressbar_module, "_best_tqdm_cls", None, raising=False)

    items = [1, 2, 3]
    parallel_map(lambda x: x, items, n_jobs=2, disable=True)

    assert captured_kwargs["n_jobs"] == 2
    assert captured_kwargs["backend"] == "loky"
    assert captured_kwargs["return_as"] == "list"


def test_parallel_map_joblib_kwargs_overrides(monkeypatch: pytest.MonkeyPatch) -> None:
    captured_kwargs = {}

    class StubParallel:
        def __init__(self, **kwargs: Any) -> None:
            nonlocal captured_kwargs
            captured_kwargs = kwargs

        def __call__(self, iterable: Any) -> list[Any]:
            return [None] * len(list(iterable))

    monkeypatch.setattr(progressbar_module, "Parallel", StubParallel)
    monkeypatch.setattr(progressbar_module, "_best_tqdm_cls", None, raising=False)

    items = [1, 2, 3]
    parallel_map(lambda x: x, items, n_jobs=2, disable=True, batch_size=10)

    assert captured_kwargs["n_jobs"] == 2
    assert captured_kwargs["backend"] == "loky"
    assert captured_kwargs["batch_size"] == 10
    assert captured_kwargs["return_as"] == "generator"


def test_parallel_map_normalizes_none_backend_to_loky(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured_kwargs: dict[str, Any] = {}

    class StubParallel:
        def __init__(self, **kwargs: Any) -> None:
            captured_kwargs.update(kwargs)

        def __call__(self, iterable: Any) -> list[Any]:
            return [None] * len(list(iterable))

    monkeypatch.setattr(progressbar_module, "Parallel", StubParallel)
    monkeypatch.setattr(progressbar_module, "_best_tqdm_cls", None, raising=False)

    parallel_map(lambda value: value, [1, 2], n_jobs=2, disable=True, backend=None)

    assert captured_kwargs["backend"] == "loky"


@pytest.mark.parametrize(
    ("backend", "expected_backend"),
    [("threading", None), ("multiprocessing", "loky")],
)
def test_parallel_map_preserves_threading_but_normalizes_process_backend_context(
    backend: str,
    expected_backend: str | None,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured_kwargs: dict[str, Any] = {}

    class StubParallel:
        def __init__(self, **kwargs: Any) -> None:
            captured_kwargs.update(kwargs)

        def __call__(self, iterable: Any) -> list[Any]:
            return [None] * len(list(iterable))

    monkeypatch.setattr(progressbar_module, "Parallel", StubParallel)
    monkeypatch.setattr(progressbar_module, "_best_tqdm_cls", None, raising=False)

    with parallel_config(backend=backend):
        parallel_map(lambda value: value, [1, 2], n_jobs=2, disable=True)

    if expected_backend is None:
        assert "backend" not in captured_kwargs
    else:
        assert captured_kwargs["backend"] == expected_backend


def test_native_guard_discards_joblib_reusable_executor() -> None:
    from joblib.externals.loky import reusable_executor

    parallel_map(lambda value: value + 1, [1, 2], n_jobs=2, disable=True, return_results=True)
    assert reusable_executor._executor is not None

    with native_reconstruction_guard():
        pass

    assert reusable_executor._executor is None


def test_worker_detection_covers_joblib_multiprocessing_children(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        progressbar_module.multiprocessing,
        "current_process",
        lambda: SimpleNamespace(name="ForkPoolWorker-1", _identity=(1,)),
    )

    assert progressbar_module.is_loky_worker() is True


def test_worker_detection_does_not_classify_main_process(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        progressbar_module.multiprocessing,
        "current_process",
        lambda: SimpleNamespace(name="MainProcess", _identity=()),
    )

    assert progressbar_module.is_loky_worker() is False


def test_native_guard_rejects_pending_external_loky_work(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from joblib.externals.loky import reusable_executor

    class PendingExecutor:
        _pending_work_items = {"work": object()}
        _running_work_items: set[int] = set()
        shutdown_called = False

        def shutdown(self, **_kwargs: Any) -> None:
            self.shutdown_called = True

    executor = PendingExecutor()
    monkeypatch.setattr(reusable_executor, "_executor", executor)

    with (
        pytest.raises(RuntimeError, match="external joblib/loky task"),
        native_reconstruction_guard(),
    ):
        pass

    assert executor.shutdown_called is False


def test_native_guard_rejects_unmanaged_child_processes(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        progressbar_module.multiprocessing,
        "active_children",
        lambda: [SimpleNamespace(name="external-worker")],
    )

    with (
        pytest.raises(RuntimeError, match="unmanaged child processes"),
        native_reconstruction_guard(),
    ):
        pass


def test_native_guard_delegates_pid_validation_to_molgr(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls: list[str] = []
    monkeypatch.setattr(
        progressbar_module,
        "ensure_current_process",
        lambda api: calls.append(api),
    )

    with native_reconstruction_guard():
        pass

    assert calls == ["MolOP native reconstruction"]


def test_native_guard_surfaces_molgr_fork_rejection(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def reject(_api: str) -> None:
        raise RuntimeError("forked child")

    monkeypatch.setattr(
        progressbar_module,
        "ensure_current_process",
        reject,
    )

    with (
        pytest.raises(NativeReconstructionConcurrencyError, match="forked child"),
        progressbar_module.native_reconstruction_guard(),
    ):
        pass


def test_native_guard_rejects_reentrant_entry() -> None:
    with native_reconstruction_guard():  # noqa: SIM117 - outer guard must remain active
        with (
            pytest.raises(NativeReconstructionConcurrencyError, match="cannot be re-entered"),
            native_reconstruction_guard(),
        ):
            pass


def test_parallel_guard_rejects_nested_parallelism_inside_native() -> None:
    with (
        native_reconstruction_guard(),
        pytest.raises(NativeReconstructionConcurrencyError, match="inside the MolGR"),
        loky_parallel_guard(),
    ):
        pass


def test_parallel_guard_fails_fast_when_native_is_active_in_another_thread() -> None:
    entered = threading.Event()

    def attempt_parallel() -> None:
        entered.set()
        with loky_parallel_guard():
            pass

    with native_reconstruction_guard(), ThreadPoolExecutor(max_workers=1) as executor:
        future = executor.submit(attempt_parallel)
        assert entered.wait(1)
        with pytest.raises(NativeReconstructionConcurrencyError, match="native"):
            future.result(timeout=1)


def test_native_guard_fails_fast_when_another_thread_holds_native_guard() -> None:
    entered = threading.Event()

    def attempt_native() -> None:
        entered.set()
        with native_reconstruction_guard():
            pass

    with native_reconstruction_guard(), ThreadPoolExecutor(max_workers=1) as executor:
        future = executor.submit(attempt_native)
        assert entered.wait(1)
        with pytest.raises(NativeReconstructionConcurrencyError, match="already active"):
            future.result(timeout=1)


def test_native_guard_rejects_any_active_parallel_thread() -> None:
    with (
        loky_parallel_guard(),
        pytest.raises(
            NativeReconstructionConcurrencyError,
            match="MolOP-managed.*task or result generator",
        ),
        native_reconstruction_guard(),
    ):
        pass


def test_parallel_map_multiprocessing_backend_uses_supported_return_mode(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured_kwargs: dict[str, Any] = {}

    class StubParallel:
        def __init__(self, **kwargs: Any) -> None:
            captured_kwargs.update(kwargs)

        def __call__(self, iterable: Any) -> list[Any]:
            return [item() for item in iterable]

    monkeypatch.setattr(progressbar_module, "Parallel", StubParallel)
    monkeypatch.setattr(
        progressbar_module,
        "delayed",
        lambda func: lambda item: lambda: func(item),
    )
    monkeypatch.setattr(progressbar_module, "_best_tqdm_cls", None, raising=False)

    result = parallel_map(
        lambda value: value + 1,
        [1, 2],
        n_jobs=2,
        backend="multiprocessing",
        disable=True,
        return_results=True,
    )

    assert result == [2, 3]
    assert captured_kwargs["return_as"] == "list"
