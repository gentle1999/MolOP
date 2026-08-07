"""
Author: TMJ
Date: 2025-07-29 19:04:21
LastEditors: TMJ
LastEditTime: 2026-03-23 10:59:00
Description: 请填写简介
"""

from abc import abstractmethod
from collections.abc import Mapping
from dataclasses import dataclass
from types import MappingProxyType
from typing import Any, Generic, Protocol, TypeVar, cast

from pydantic import Field, PrivateAttr

from molop.io.base_models.Bases import BaseDataClassWithUnit
from molop.io.base_models.ChemFileFrame import BaseChemFileFrame
from molop.io.base_models.summary import SummaryDict, summary_column
from molop.io.codec_types import ParseOptions


FrameT = TypeVar("FrameT", bound=BaseChemFileFrame)


class _HasParseMethod(Protocol):
    only_extract_structure: bool
    capture_source_evidence: bool

    def parse(self, block: str, *, additional_data: dict[str, Any] | None = None) -> Any: ...


@dataclass(frozen=True, slots=True)
class FrameParseContext:
    """File-level metadata supplied to one frame parse."""

    additional_data: Mapping[str, Any]
    parse_options: ParseOptions | None = None

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "additional_data",
            MappingProxyType(dict(self.additional_data)),
        )


class BaseFrameParser(BaseDataClassWithUnit, Generic[FrameT]):
    only_extract_structure: bool = Field(default=False, exclude=True, repr=False)
    capture_source_evidence: bool = Field(default=False, exclude=True, repr=False)
    parse_options: ParseOptions | None = Field(default=None, exclude=True, repr=False)
    _file_frame_class_: type[FrameT] = PrivateAttr()

    def parse(self, block: str, *, additional_data: dict[str, Any] | None = None) -> FrameT:
        context = FrameParseContext(
            additional_data=dict(additional_data) if additional_data is not None else {},
            parse_options=self.parse_options,
        )
        temp_dict = {"frame_content": block}
        temp_dict.update(self._parse_frame(block, context=context))
        temp_dict.update(context.additional_data)
        validation_context = (
            {"force_unit_transform": context.parse_options.force_unit_transform}
            if context.parse_options is not None
            else None
        )
        return cast(
            FrameT,
            self._file_frame_class_.model_validate(temp_dict, context=validation_context),
        )

    @abstractmethod
    def _parse_frame(
        self,
        block: str,
        *,
        context: FrameParseContext,
    ) -> Mapping[str, Any]:
        raise NotImplementedError()

    def to_summary_dict(self, **kwargs) -> SummaryDict:
        return {
            summary_column("FrameParser", "only_extract_structure"): self.only_extract_structure,
            summary_column("FrameParser", "capture_source_evidence"): (
                self.capture_source_evidence
            ),
        }


FrameParser = TypeVar("FrameParser", bound="BaseFrameParser[Any]")
