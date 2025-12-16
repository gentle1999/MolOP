"""
Author: TMJ
Date: 2025-07-29 19:04:21
LastEditors: TMJ
LastEditTime: 2025-11-30 18:49:52
Description: 请填写简介
"""

from abc import abstractmethod
from collections.abc import Mapping
from typing import Any, Generic, TypeVar

from pydantic import Field, PrivateAttr

from molop.io.base_models.Bases import BaseDataClassWithUnit
from molop.io.base_models.ChemFileFrame import ChemFileFrame


class BaseFrameParser(BaseDataClassWithUnit, Generic[ChemFileFrame]):
    only_extract_structure: bool = Field(default=False, exclude=True, repr=False)
    _file_frame_class_: type[ChemFileFrame] = PrivateAttr()
    _block: str = ""

    def parse(self, block: str, *, additional_data: dict[str, Any] | None = None) -> ChemFileFrame:
        self._block = block
        temp_dict = {"frame_content": block}
        if additional_data is not None:
            temp_dict.update(additional_data)
        temp_dict.update(self._parse_frame())
        return self._file_frame_class_.model_validate(temp_dict)

    @abstractmethod
    def _parse_frame(self) -> Mapping[str, Any]:
        raise NotImplementedError()

    def _add_default_units(self) -> None: ...

    def to_summary_dict(self, **kwargs) -> dict[tuple[str, str], Any]:
        return {("FrameParser", "only_extract_structure"): self.only_extract_structure}


FrameParser = TypeVar("FrameParser", bound="BaseFrameParser")
