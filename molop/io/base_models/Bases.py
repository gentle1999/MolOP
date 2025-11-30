"""
Author: TMJ
Date: 2025-07-26 19:08:33
LastEditors: TMJ
LastEditTime: 2025-11-30 18:46:56
Description: 请填写简介
"""

from abc import abstractmethod
from typing import Any, Dict, Mapping, Union, overload

import numpy as np
import pandas as pd
from pint._typing import Magnitude, UnitLike
from pint.facets.plain import PlainQuantity
from pydantic import BaseModel, ConfigDict, PrivateAttr, model_validator
from typing_extensions import Self

from molop.config import molopconfig, moloplogger
from molop.unit import unit_transform


class BaseDataClassWithUnit(BaseModel):
    model_config = ConfigDict(arbitrary_types_allowed=True)
    _default_units: Dict[str, UnitLike] = PrivateAttr(default_factory=dict)
    _set_default_units: bool = PrivateAttr(default=False)

    @staticmethod
    @overload
    def __unitless_dump__(item: PlainQuantity) -> Magnitude: ...
    @staticmethod
    @overload
    def __unitless_dump__(item: np.ndarray) -> list: ...
    @staticmethod
    @overload
    def __unitless_dump__(item: Magnitude) -> Magnitude: ...
    @staticmethod
    @overload
    def __unitless_dump__(item: str) -> str: ...
    @staticmethod
    @overload
    def __unitless_dump__(item: list) -> list: ...
    @staticmethod
    @overload
    def __unitless_dump__(item: tuple) -> tuple: ...
    @staticmethod
    @overload
    def __unitless_dump__(item: Mapping[str, Any]) -> Mapping[str, Any]: ...
    @staticmethod
    @overload
    def __unitless_dump__(
        item: "BaseDataClassWithUnit",
    ) -> Union[str, list, tuple, Mapping[str, Any], Magnitude]: ...
    @staticmethod
    def __unitless_dump__(
        item: Union[
            PlainQuantity,
            np.ndarray,
            str,
            list,
            tuple,
            Magnitude,
            Mapping[str, Any],
            "BaseDataClassWithUnit",
        ],
        **kwargs,
    ) -> Union[str, list, tuple, Mapping[str, Any], Magnitude]:
        if isinstance(item, PlainQuantity):
            if isinstance(item.m, np.ndarray):
                return item.m.tolist()
            else:
                return item.m
        if isinstance(item, np.ndarray):
            return item.tolist()
        if isinstance(item, list):
            return [BaseDataClassWithUnit.__unitless_dump__(i, **kwargs) for i in item]
        if isinstance(item, tuple):
            return tuple(
                BaseDataClassWithUnit.__unitless_dump__(i, **kwargs) for i in item
            )
        if isinstance(item, Mapping):
            return {
                k: BaseDataClassWithUnit.__unitless_dump__(v, **kwargs)
                for k, v in item.items()
            }
        if isinstance(item, BaseDataClassWithUnit):
            return item.to_unitless_dump(**kwargs)
        return item

    def to_unitless_dump(self, **kwargs) -> Dict[str, Any]:
        """
        Parameters:
            kwargs:
                follow the same format as
                [`model_dump`](https://docs.pydantic.dev/latest/concepts/serialization/) method.
        """
        return {
            k: BaseDataClassWithUnit.__unitless_dump__(getattr(self, k), **kwargs)
            for k, v in self.model_dump(**kwargs).items()
        }

    @model_validator(mode="after")
    def __unit_transform__(self) -> Self:
        self._add_default_units()
        if self._set_default_units or molopconfig.force_unit_transform:
            self._transform_units(self._default_units)
            moloplogger.debug(f"Data class {self.__class__.__name__} parsed.\n{self}")
        return self

    @abstractmethod
    def _add_default_units(self) -> None:
        raise NotImplementedError

    def _transform_units(self, unit_dict: Mapping[str, UnitLike]) -> None:
        for key, unit in unit_dict.items():
            if hasattr(self, key):
                setattr(self, key, unit_transform(getattr(self, key), unit))

    def to_summary_series(self, **kwargs) -> pd.Series:
        return pd.Series(self.to_summary_dict(**kwargs))

    @abstractmethod
    def to_summary_dict(self, **kwargs) -> Dict[tuple[str, str], Any]:
        raise NotImplementedError
