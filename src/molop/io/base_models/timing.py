"""Shared wall-clock aliases and independent totals for calculation timing."""

from collections.abc import Mapping, Sequence
from math import isclose
from typing import Any

from pint.facets.plain import PlainQuantity


def resolve_wall_time(
    running_time: PlainQuantity | None, wall_time: PlainQuantity | None
) -> PlainQuantity | None:
    """Resolve the legacy running_time name as an alias of wall_time."""
    if (
        running_time is not None
        and wall_time is not None
        and not isclose(
            float(running_time.m_as("second")),
            float(wall_time.m_as("second")),
            rel_tol=1e-12,
            abs_tol=1e-12,
        )
    ):
        raise ValueError("running_time and wall_time must describe the same wall-clock duration")
    return wall_time if wall_time is not None else running_time


def accumulate_calculation_timings(
    segments: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    """Sum like quantities once per segment; incomplete totals remain unknown."""
    result: dict[str, Any] = {}
    for field in ("cpu_time", "wall_time"):
        values = [value for segment in segments if (value := segment.get(field)) is not None]
        result[field] = (
            sum(values[1:], values[0]) if values and len(values) == len(segments) else None
        )
    result["running_time"] = result["wall_time"]
    sources = [segment.get("cpu_time_source") for segment in segments]
    result["cpu_time_source"] = None
    if result["cpu_time"] is not None:
        if "estimated" in sources:
            result["cpu_time_source"] = "estimated"
        elif sources and all(source == "observed" for source in sources):
            result["cpu_time_source"] = "observed"
    counts = [segment.get("cpu_time_num_processes") for segment in segments]
    result["cpu_time_num_processes"] = (
        counts[0] if counts and all(count == counts[0] for count in counts) else None
    )
    return result
