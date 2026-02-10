import sys
from types import SimpleNamespace

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
