import sys
import warnings
from collections.abc import Callable, Iterable
from typing import TYPE_CHECKING, Any, NoReturn, TypeAlias, TypeVar, cast, overload

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
            from tqdm.notebook import tqdm_notebook as nb_tqdm

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
        obj = tqdm_factory(iterable, *args, **kwargs)
    if iterable is None:
        return cast(TqdmObj, obj)
    return cast(Iterable[T], obj)


def parallel_map(
    func: Callable[[T], R],
    iterable: Iterable[T],
    n_jobs: int = -1,
    total: int | None = None,
    desc: str = "Processing",
    disable: bool = False,
    tqdm_kwargs: dict[str, Any] | None = None,
    **joblib_kwargs: Any,
) -> list[R]:
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
    **joblib_kwargs : Any
        Additional keyword arguments to pass to joblib.Parallel.

    Returns
    -------
    list[R]
        The list of results after applying the function to each item in the iterable.
    """
    iterator = AdaptiveProgress(
        iterable, desc=desc, total=total, disable=disable, **(tqdm_kwargs or {})
    )
    if n_jobs == 1:
        return [func(item) for item in iterator]

    results = Parallel(n_jobs=n_jobs, return_as="list", **joblib_kwargs)(
        delayed(func)(item) for item in iterator
    )
    if results is None:
        raise ValueError("The parallel map returned None.")
    return cast(list[R], list(results))
