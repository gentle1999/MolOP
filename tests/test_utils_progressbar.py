import io
import json
import multiprocessing
import sys
import threading
from collections.abc import Callable
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict
from types import ModuleType, SimpleNamespace
from typing import Any

import pytest
from joblib import parallel_config
from tqdm import tqdm as st_tqdm

import molop.utils.progressbar as progressbar_module
from molop.utils.progressbar import (
    AdaptiveProgress,
    JsonLinesProgressListener,
    NativeReconstructionConcurrencyError,
    ProgressEvent,
    ProgressRecorder,
    _is_notebook,
    loky_parallel_guard,
    native_reconstruction_guard,
    parallel_map,
    progress_jsonl_sink,
    progress_listener,
    register_progress_listener,
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


# ---------------------------------------------------------------------------
# Progress event capture (third-party listeners, callbacks, JSONL sink)
# ---------------------------------------------------------------------------


def _collect_events() -> tuple[list[ProgressEvent], Callable[[ProgressEvent], None]]:
    events: list[ProgressEvent] = []

    def listener(event: ProgressEvent) -> None:
        events.append(event)

    return events, listener


def test_progress_events_emitted_for_iterable_bar(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(progressbar_module, "_best_tqdm_cls", None, raising=False)
    events, listener = _collect_events()

    with progress_listener(listener):
        yielded = list(AdaptiveProgress(["a", "b", "c"], disable=True))

    assert yielded == ["a", "b", "c"]
    kinds = [event.event for event in events]
    assert kinds == ["start", "update", "update", "update", "close"]
    updates = [event for event in events if event.event == "update"]
    assert [event.done for event in updates] == [1, 2, 3]
    assert all(event.total == 3 for event in updates)
    assert all(event.task_id == events[0].task_id for event in events)
    assert all(event.description == "" for event in events)


def test_progress_events_emitted_for_active_bar(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(progressbar_module, "_best_tqdm_cls", st_tqdm, raising=False)
    events, listener = _collect_events()

    with progress_listener(listener):
        yielded = list(AdaptiveProgress(["a", "b"], disable=False, file=io.StringIO()))

    assert yielded == ["a", "b"]
    assert [event.event for event in events] == ["start", "update", "update", "close"]
    assert [event.done for event in events if event.event == "update"] == [1, 2]


def test_progress_events_for_manual_bar(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(progressbar_module, "_best_tqdm_cls", None, raising=False)
    events, listener = _collect_events()

    with progress_listener(listener):
        bar = AdaptiveProgress(None, disable=True, desc="manual", total=10)
        bar.update(3)
        bar.total = 12
        bar.update(2)
        bar.set_description("manual-2")
        bar.close()

    kinds = [event.event for event in events]
    assert kinds == ["start", "update", "update", "description", "close"]
    assert events[1].done == 3 and events[1].total == 10
    assert events[2].done == 5 and events[2].total == 12
    assert events[3].description == "manual-2"
    assert events[4].done == 5 and events[4].total == 12


def test_progress_reset_is_supported_by_rich_backend(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    try:
        from tqdm.rich import tqdm_rich
    except ImportError:  # pragma: no cover - rich is an optional tqdm backend
        pytest.skip("rich backend is unavailable")

    monkeypatch.setattr(progressbar_module, "_best_tqdm_cls", tqdm_rich, raising=False)
    events, listener = _collect_events()

    with progress_listener(listener):
        bar = AdaptiveProgress(None, disable=False, total=2)
        bar.update(1)
        bar.reset(total=3)
        bar.close()

    assert [(event.event, event.done, event.total) for event in events] == [
        ("start", 0, 2),
        ("update", 1, 2),
        ("update", 0, 3),
        ("close", 0, 3),
    ]


def test_progress_events_capture_description_str_and_ignore_closed_updates(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(progressbar_module, "_best_tqdm_cls", None, raising=False)
    events, listener = _collect_events()

    with progress_listener(listener):
        bar = AdaptiveProgress(None, disable=True, total=2)
        bar.set_description_str("raw")
        bar.close()
        bar.update(1)
        bar.set_description("after-close")
        bar.refresh()

    assert [(event.event, event.description, event.done) for event in events] == [
        ("start", "", 0),
        ("description", "raw", 0),
        ("close", "raw", 0),
    ]


def test_progress_events_track_disabled_bars(monkeypatch: pytest.MonkeyPatch) -> None:
    # disable=True only suppresses rendering; listeners still see progress.
    monkeypatch.setattr(progressbar_module, "_best_tqdm_cls", None, raising=False)
    events, listener = _collect_events()

    with progress_listener(listener):
        list(AdaptiveProgress(range(5), desc="parse: ", disable=True))

    assert [event.done for event in events if event.event == "update"] == [1, 2, 3, 4, 5]
    assert all(event.description == "parse: " for event in events)


def test_progress_listener_register_and_unregister(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(progressbar_module, "_best_tqdm_cls", None, raising=False)
    events, listener = _collect_events()

    unregister = register_progress_listener(listener)
    bar = AdaptiveProgress(None, disable=True, total=1)
    bar.close()
    unregister()

    bar = AdaptiveProgress(None, disable=True, total=1)
    bar.close()

    assert [event.event for event in events] == ["start", "close"]


def test_progress_listener_is_deduplicated(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(progressbar_module, "_best_tqdm_cls", None, raising=False)
    events, listener = _collect_events()

    first = register_progress_listener(listener)
    second = register_progress_listener(listener)
    try:
        bar = AdaptiveProgress(None, disable=True, total=1)
        bar.close()
    finally:
        first()
        second()

    assert [event.event for event in events] == ["start", "close"]


def test_progress_listener_registration_lifetimes_are_independent(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(progressbar_module, "_best_tqdm_cls", None, raising=False)
    events, listener = _collect_events()

    first = register_progress_listener(listener)
    second = register_progress_listener(listener)
    first()
    bar = AdaptiveProgress(None, disable=True, total=1)
    bar.close()
    second()
    bar = AdaptiveProgress(None, disable=True, total=1)
    bar.close()

    assert [event.event for event in events] == ["start", "close"]


def test_progress_listener_contexts_can_be_nested(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(progressbar_module, "_best_tqdm_cls", None, raising=False)
    events, listener = _collect_events()

    with progress_listener(listener):
        with progress_listener(listener):
            bar = AdaptiveProgress(None, disable=True, total=1)
            bar.close()
        bar = AdaptiveProgress(None, disable=True, total=1)
        bar.close()

    assert [event.event for event in events] == ["start", "close", "start", "close"]


def test_progress_events_close_when_iterable_generator_is_closed(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(progressbar_module, "_best_tqdm_cls", None, raising=False)
    events, listener = _collect_events()

    with progress_listener(listener):
        bar = AdaptiveProgress([1, 2], disable=True, total=2)
        iterator = iter(bar)
        assert next(iterator) == 1
        iterator.close()

    assert [event.event for event in events] == ["start", "close"]


def test_parallel_map_closes_progress_on_serial_failure(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(progressbar_module, "_best_tqdm_cls", None, raising=False)
    recorder = ProgressRecorder()

    def fail_on_second(value: int) -> int:
        if value == 2:
            raise RuntimeError("work failed")
        return value

    with recorder, pytest.raises(RuntimeError, match="work failed"):
        parallel_map(fail_on_second, [1, 2, 3], n_jobs=1, disable=True)

    assert [event.event for event in recorder.history()] == [
        "start",
        "update",
        "close",
    ]
    assert recorder.snapshot() == {}


def test_parallel_map_does_not_start_unconsumed_generator_progress(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(progressbar_module, "_best_tqdm_cls", None, raising=False)
    recorder = ProgressRecorder()

    with recorder:
        results = parallel_map(
            lambda value: value * 2,
            [1, 2, 3],
            n_jobs=1,
            return_as="generator",
            disable=True,
        )
        assert recorder.history() == ()
        results.close()  # type: ignore[union-attr]

    assert recorder.snapshot() == {}
    assert recorder.history() == ()


def test_progress_refresh_emits_changed_total(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(progressbar_module, "_best_tqdm_cls", None, raising=False)
    recorder = ProgressRecorder()

    with recorder:
        bar = AdaptiveProgress(None, disable=True)
        task_id = next(iter(recorder.snapshot()))
        bar.total = 3
        bar.refresh()
        assert recorder.latest(task_id) is not None
        assert recorder.latest(task_id).total == 3  # type: ignore[union-attr]
        bar.close()

    assert [(event.event, event.total) for event in recorder.history()] == [
        ("start", None),
        ("update", 3),
        ("close", 3),
    ]


def test_progress_listener_exception_does_not_break_bar(
    monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
) -> None:
    monkeypatch.setattr(progressbar_module, "_best_tqdm_cls", None, raising=False)

    def broken(_event: ProgressEvent) -> None:
        raise RuntimeError("listener boom")

    with progress_listener(broken):
        yielded = list(AdaptiveProgress([1, 2], disable=True))

    assert yielded == [1, 2]
    assert any(
        record.exc_info is not None and "listener boom" in str(record.exc_info[1])
        for record in caplog.records
    )


def test_progress_callback_receives_events(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(progressbar_module, "_best_tqdm_cls", None, raising=False)
    events, callback = _collect_events()

    list(AdaptiveProgress(["a", "b"], disable=True, progress_callback=callback))

    assert [event.event for event in events] == ["start", "update", "update", "close"]


def test_parallel_map_progress_callback_emits_events(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(progressbar_module, "_best_tqdm_cls", None, raising=False)
    events, callback = _collect_events()

    results = parallel_map(
        lambda value: value * 2,
        [1, 2, 3],
        n_jobs=1,
        disable=True,
        progress_callback=callback,
    )

    assert results == [2, 4, 6]
    kinds = [event.event for event in events]
    assert kinds[0] == "start" and kinds[-1] == "close"
    assert [event.done for event in events if event.event == "update"] == [1, 2, 3]


def test_parallel_map_explicit_progress_callback_wins_over_tqdm_kwargs(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(progressbar_module, "_best_tqdm_cls", None, raising=False)
    explicit_events, explicit_callback = _collect_events()
    legacy_events, legacy_callback = _collect_events()

    parallel_map(
        lambda value: value * 2,
        [1, 2],
        n_jobs=1,
        disable=True,
        progress_callback=explicit_callback,
        tqdm_kwargs={"progress_callback": legacy_callback},
    )

    assert len(explicit_events) == 4
    assert legacy_events == []


def test_parallel_map_progress_tracks_completed_batches(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(progressbar_module, "_best_tqdm_cls", None, raising=False)
    events, callback = _collect_events()

    results = list(
        parallel_map(
            lambda value: value * 2,
            [1, 2, 3, 4],
            n_jobs=2,
            backend="threading",
            batch_size=1,
            disable=True,
            progress_callback=callback,
            return_results=True,
        )
    )

    assert results == [2, 4, 6, 8]
    updates = [event.done for event in events if event.event == "update"]
    assert updates == [1, 2, 3, 4]
    assert events[-1].event == "close"


def test_progress_events_are_json_serializable(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(progressbar_module, "_best_tqdm_cls", None, raising=False)
    events, listener = _collect_events()

    with progress_listener(listener):
        bar = AdaptiveProgress(None, disable=True, desc="task", total=2)
        bar.update(1)
        bar.close()

    for event in events:
        json.dumps(asdict(event))


def test_progress_jsonl_sink_writes_json_lines(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Any
) -> None:
    monkeypatch.setattr(progressbar_module, "_best_tqdm_cls", None, raising=False)
    path = tmp_path / "progress.jsonl"

    with progress_jsonl_sink(path):
        list(AdaptiveProgress(["a", "b"], disable=True))

    lines = path.read_text(encoding="utf-8").strip().splitlines()
    records = [json.loads(line) for line in lines]
    assert [record["event"] for record in records] == [
        "start",
        "update",
        "update",
        "close",
    ]
    assert [record["done"] for record in records if record["event"] == "update"] == [1, 2]
    assert all(
        set(record) == {"task_id", "event", "description", "done", "total", "timestamp"}
        for record in records
    )


def test_progress_jsonl_sink_appends_and_closes(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Any
) -> None:
    monkeypatch.setattr(progressbar_module, "_best_tqdm_cls", None, raising=False)
    path = tmp_path / "progress.jsonl"

    with progress_jsonl_sink(path):
        bar = AdaptiveProgress(None, disable=True, total=1)
        bar.close()

    # A second sink appends to the same file rather than truncating it.
    with progress_jsonl_sink(path):
        bar = AdaptiveProgress(None, disable=True, total=1)
        bar.close()

    lines = path.read_text(encoding="utf-8").strip().splitlines()
    assert len(lines) == 4
    assert all(json.loads(line)["event"] in {"start", "close"} for line in lines)


def test_json_lines_listener_close_is_idempotent(tmp_path: Any) -> None:
    listener = JsonLinesProgressListener(tmp_path / "nested" / "progress.jsonl")
    listener(
        ProgressEvent(task_id=1, event="start", description="", done=0, total=1, timestamp=0.0)
    )
    listener.close()
    listener.close()

    lines = (
        (tmp_path / "nested" / "progress.jsonl").read_text(encoding="utf-8").strip().splitlines()
    )
    assert len(lines) == 1


# ---------------------------------------------------------------------------
# Official ProgressRecorder
# ---------------------------------------------------------------------------


def test_progress_recorder_tracks_active_bar(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(progressbar_module, "_best_tqdm_cls", None, raising=False)
    rec = ProgressRecorder()

    with rec:
        bar = AdaptiveProgress(None, disable=True, desc="task", total=10)
        bar.update(3)
        active = rec.snapshot()
        assert len(active) == 1
        task_id = next(iter(active))
        assert active[task_id].done == 3
        assert active[task_id].total == 10
        assert rec.latest(task_id) is active[task_id]
        assert rec.percent(task_id) == 30.0
        assert rec.active_tasks() == (task_id,)
        bar.close()

    assert rec.snapshot() == {}
    assert rec.active_tasks() == ()
    assert [event.event for event in rec.history()] == ["start", "update", "close"]


def test_progress_recorder_handles_concurrent_bars(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(progressbar_module, "_best_tqdm_cls", None, raising=False)
    rec = ProgressRecorder()

    with rec:
        bar_a = AdaptiveProgress(None, disable=True, desc="a", total=5)
        bar_b = AdaptiveProgress(None, disable=True, desc="b", total=5)
        bar_a.update(2)
        bar_b.update(4)
        snapshot = rec.snapshot()
        assert len(snapshot) == 2
        bar_a.close()
        assert len(rec.snapshot()) == 1
        bar_b.close()

    assert rec.snapshot() == {}
    assert len(rec.history()) == 2 * 3  # start/close per bar, no updates after close


def test_progress_recorder_percent_unknown_total(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(progressbar_module, "_best_tqdm_cls", None, raising=False)
    rec = ProgressRecorder()

    with rec:
        bar = AdaptiveProgress(None, disable=True)  # total unknown
        bar.update(3)
        task_id = next(iter(rec.snapshot()))
        assert rec.percent(task_id) is None
        bar.close()

    assert rec.percent(task_id) is None  # bar no longer live


def test_progress_recorder_keep_history_false(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(progressbar_module, "_best_tqdm_cls", None, raising=False)
    rec = ProgressRecorder(keep_history=False)

    with rec:
        bar = AdaptiveProgress(None, disable=True, total=1)
        bar.update(1)
        bar.close()

    assert rec.history() == ()
    assert rec.snapshot() == {}


def test_progress_recorder_clear(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(progressbar_module, "_best_tqdm_cls", None, raising=False)
    rec = ProgressRecorder()

    with rec:
        bar = AdaptiveProgress(None, disable=True, total=1)
        bar.update(1)
        assert rec.snapshot() != {}
        rec.clear()
        assert rec.snapshot() == {}
        assert rec.history() == ()
        bar.close()


def test_progress_recorder_context_manager_unregisters(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(progressbar_module, "_best_tqdm_cls", None, raising=False)
    rec = ProgressRecorder()

    with rec:
        bar = AdaptiveProgress(None, disable=True, total=1)
        bar.close()
    assert len(rec.history()) == 2

    # Outside the with-block the recorder must no longer receive events.
    bar = AdaptiveProgress(None, disable=True, total=1)
    bar.close()
    assert len(rec.history()) == 2


def test_progress_recorder_snapshot_is_thread_safe(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(progressbar_module, "_best_tqdm_cls", None, raising=False)
    rec = ProgressRecorder()
    stop = threading.Event()
    errors: list[BaseException] = []

    def reader() -> None:
        try:
            while not stop.is_set():
                rec.snapshot()
                rec.history()
                rec.active_tasks()
        except BaseException as exc:  # pragma: no cover - failure path
            errors.append(exc)

    with rec:
        reader_thread = threading.Thread(target=reader, daemon=True)
        reader_thread.start()
        bar = AdaptiveProgress(None, disable=True, total=100)
        for _ in range(100):
            bar.update(1)
        bar.close()
        stop.set()
        reader_thread.join()

    assert errors == []
    assert rec.snapshot() == {}
    assert len(rec.history()) == 1 + 100 + 1
