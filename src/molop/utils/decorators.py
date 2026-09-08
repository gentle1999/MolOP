"""Reusable decorators for public MolOP APIs."""

from __future__ import annotations

import warnings
from collections.abc import Callable
from functools import wraps
from typing import ParamSpec, TypeVar


P = ParamSpec("P")
R = TypeVar("R")


class ExperimentalWarning(UserWarning):
    """Warning emitted when an experimental public API is called."""


def experimental(func: Callable[P, R]) -> Callable[P, R]:
    """Mark a callable as experimental and warn at each call site."""

    @wraps(func)
    def wrapped(*args: P.args, **kwargs: P.kwargs) -> R:
        warnings.warn(
            f"Experimental API {func.__name__!r} is not stable; "
            "its behavior and interface may change without notice.",
            category=ExperimentalWarning,
            stacklevel=2,
        )
        return func(*args, **kwargs)

    setattr(wrapped, "__experimental__", True)  # noqa: B010
    return wrapped


__all__ = ["ExperimentalWarning", "experimental"]
