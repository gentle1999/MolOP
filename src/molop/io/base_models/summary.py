from __future__ import annotations

from collections.abc import Iterable
from typing import Any, TypeAlias

import numpy as np
import pandas as pd


SummaryColumn: TypeAlias = tuple[str, str, str]
SummaryDict: TypeAlias = dict[SummaryColumn, Any]


def is_quantity(value: Any) -> bool:
    return hasattr(value, "magnitude") and hasattr(value, "units")


def summary_column(group: str, field: str, unit: str | None = "") -> SummaryColumn:
    return (group, field, "" if unit is None else str(unit))


def summary_value(value: Any) -> Any:
    if is_quantity(value):
        magnitude = value.magnitude
        if isinstance(magnitude, np.ndarray):
            return magnitude.tolist()
        return magnitude
    if isinstance(value, np.ndarray):
        return value.tolist()
    return value


def summary_item(group: str, field: str, value: Any) -> tuple[SummaryColumn, Any] | None:
    if value is None:
        return None
    if is_quantity(value):
        return summary_column(group, field, str(value.units)), summary_value(value)
    return summary_column(group, field), summary_value(value)


def summary_dict_from_fields(obj: Any, group: str, **kwargs: Any) -> SummaryDict:
    summary: SummaryDict = {}
    for field in obj.model_dump(**kwargs):
        item = summary_item(group, field, getattr(obj, field))
        if item is not None:
            column, value = item
            summary[column] = value
    return summary


def normalize_summary_column(column: Any) -> SummaryColumn:
    if isinstance(column, tuple):
        if len(column) == 3:
            return (str(column[0]), str(column[1]), "" if column[2] is None else str(column[2]))
        if len(column) == 2:
            return (str(column[0]), str(column[1]), "")
    return ("", str(column), "")


def normalize_summary_series(series: pd.Series) -> pd.Series:
    normalized = series.copy()
    normalized.index = pd.MultiIndex.from_tuples(
        [normalize_summary_column(column) for column in normalized.index],
        names=["group", "field", "unit"],
    )
    return normalized


def flatten_summary_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Return a copy with dot-separated column names when columns are a MultiIndex."""

    if not isinstance(df.columns, pd.MultiIndex):
        return df
    flattened = df.copy()
    flattened.columns = [
        ".".join(str(level) for level in column if level not in ("", None))
        for column in flattened.columns.to_flat_index()
    ]
    return flattened


def build_summary_df(
    series: Iterable[pd.Series],
    *,
    flatten_columns: bool = False,
) -> pd.DataFrame:
    """Build a summary DataFrame with stable column group ordering."""

    series_list = list(series)
    if not series_list:
        return pd.DataFrame()
    series_list = [normalize_summary_series(item) for item in series_list]
    df: pd.DataFrame = pd.concat(series_list, axis=1).T
    if isinstance(df.columns, pd.MultiIndex):
        top_level_order = df.columns.get_level_values(0).unique()
        df = pd.DataFrame(df.loc[:, top_level_order])
    if flatten_columns:
        return flatten_summary_columns(df)
    return df
