from typing import TYPE_CHECKING, Any, Callable, Iterable, NoReturn, Optional, TypeVar

from joblib import Parallel, delayed
from tqdm import tqdm as st_tqdm

T = TypeVar("T")
R = TypeVar("R")

if TYPE_CHECKING:
    from tqdm.notebook import tqdm_notebook
    from tqdm.rich import tqdm_rich

    TqdmType = (
        type[st_tqdm[NoReturn]]
        | type[tqdm_notebook[NoReturn]]
        | type[tqdm_rich[NoReturn]]
    )
    TqdmObj = st_tqdm[NoReturn] | tqdm_notebook[NoReturn] | tqdm_rich[NoReturn]
else:
    TqdmType = Any
    TqdmObj = Any


def _is_notebook() -> bool:
    try:
        from IPython.core.getipython import get_ipython

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


class AdaptiveProgress:
    _best_tqdm_cls: Optional[TqdmType] = None

    def __new__(cls, iterable=None, *args, **kwargs) -> TqdmObj:
        if cls._best_tqdm_cls is None:
            cls._best_tqdm_cls = cls._detect_best_backend()
        return cls._best_tqdm_cls(iterable, *args, **kwargs)

    @classmethod
    def _detect_best_backend(cls) -> TqdmType:
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


def parallel_map(
    func: Callable[[T], R],
    iterable: Iterable[T],
    n_jobs: int = -1,
    total: Optional[int] = None,
    desc: str = "Processing",
    disable: bool = False,
    tqdm_kwargs: Optional[dict[str, Any]] = None,
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
    return Parallel(n_jobs=n_jobs, return_as="list", **joblib_kwargs)(
        delayed(func)(item) for item in iterator
    )  # pyright: ignore[reportReturnType]
