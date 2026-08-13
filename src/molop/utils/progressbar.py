import importlib
import sys
import warnings
from collections.abc import Callable, Iterable
from typing import TYPE_CHECKING, Any, Literal, NoReturn, TypeAlias, TypeVar, cast, overload

from joblib import Parallel, delayed
from tqdm import tqdm as st_tqdm


try:
    from tqdm.std import TqdmExperimentalWarning
except ImportError:
    TqdmExperimentalWarning = None


T = TypeVar("T")
R = TypeVar("R")

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
    effective_return_as = return_as
    if effective_return_as is None:
        effective_return_as = "generator" if joblib_kwargs else "list"

    iterator = AdaptiveProgress(
        iterable, desc=desc, total=total, disable=disable, **(tqdm_kwargs or {})
    )
    if n_jobs == 1:
        if effective_return_as in {"generator", "generator_unordered"}:
            return _finalize_parallel_results((func(item) for item in iterator), return_results)
        return _finalize_parallel_results([func(item) for item in iterator], return_results)

    results = Parallel(n_jobs=n_jobs, return_as=effective_return_as, **joblib_kwargs)(
        delayed(func)(item) for item in iterator
    )
    if results is None:
        raise ValueError("The parallel map returned None.")
    return _finalize_parallel_results(cast(Iterable[R], results), return_results)
