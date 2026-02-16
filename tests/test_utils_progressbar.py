import sys
from types import SimpleNamespace
from typing import Any

import pytest

import molop.utils.progressbar as progressbar_module
from molop.utils.progressbar import AdaptiveProgress, _is_notebook, parallel_map


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
    assert captured_kwargs["batch_size"] == 10
    assert captured_kwargs["return_as"] == "list"
