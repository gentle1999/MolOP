"""
Author: TMJ
Date: 2025-07-26 19:08:33
LastEditors: TMJ
LastEditTime: 2026-02-04 15:47:57
Description: 请填写简介
"""

from abc import abstractmethod
from collections.abc import Mapping
from typing import Any, ClassVar

import numpy as np
import pandas as pd
from pint._typing import UnitLike
from pydantic import BaseModel, ConfigDict, model_validator
from typing_extensions import Self

from molop.config import molopconfig, moloplogger
from molop.unit import unit_transform


class BaseDataClassWithUnit(BaseModel):
    model_config = ConfigDict(arbitrary_types_allowed=True)
    default_units: ClassVar[dict[str, UnitLike]] = {}
    set_default_units: ClassVar[bool] = False

    @staticmethod
    def __unitless_dump__(item: Any, **kwargs) -> Any:
        if hasattr(item, "m") and hasattr(item, "units"):
            if isinstance(item.m, np.ndarray):
                return item.m.tolist()
            else:
                return item.m
        if isinstance(item, np.ndarray):
            return item.tolist()
        if isinstance(item, list):
            return [BaseDataClassWithUnit.__unitless_dump__(i, **kwargs) for i in item]
        if isinstance(item, tuple):
            return tuple(BaseDataClassWithUnit.__unitless_dump__(i, **kwargs) for i in item)
        if isinstance(item, Mapping):
            return {
                k: BaseDataClassWithUnit.__unitless_dump__(v, **kwargs) for k, v in item.items()
            }
        if isinstance(item, BaseDataClassWithUnit):
            return item.to_unitless_dump(**kwargs)
        return item

    def to_unitless_dump(self, **kwargs) -> dict[str, Any]:
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
        if self.set_default_units or molopconfig.force_unit_transform:
            self._transform_units(self.default_units)
            moloplogger.debug(f"Data class {self.__class__.__name__} parsed.\n{self}")
        return self

    def _transform_units(self, unit_dict: Mapping[str, UnitLike]) -> None:
        for key, unit in unit_dict.items():
            if hasattr(self, key):
                setattr(self, key, unit_transform(getattr(self, key), unit))

    def to_summary_series(self, **kwargs) -> pd.Series:
        return pd.Series(self.to_summary_dict(**kwargs))

    @abstractmethod
    def to_summary_dict(self, **kwargs) -> dict[tuple[str, str], Any]:
        raise NotImplementedError

    @classmethod
    def __pydantic_init_subclass__(cls, **kwargs):
        super().__pydantic_init_subclass__(**kwargs)
        final_default_units: dict[str, UnitLike] = {}
        for base in reversed(cls.__mro__):
            if hasattr(base, "default_units"):
                base_default_units = getattr(base, "default_units", {})
                final_default_units.update(base_default_units)
        cls.default_units = final_default_units
