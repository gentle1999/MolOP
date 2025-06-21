"""
Author: TMJ
Date: 2025-02-16 16:12:34
LastEditors: TMJ
LastEditTime: 2025-02-17 22:54:21
Description: 请填写简介
"""

import re
from typing import Union

from pydantic import BaseModel, ConfigDict, Field, PrivateAttr, model_validator
from typing_extensions import Self


class MolOPPattern(BaseModel):
    model_config = ConfigDict(arbitrary_types_allowed=True)
    start_pattern: Union[str, None] = Field(
        default=None, description="The start pattern of the pattern."
    )
    start_offset: int = Field(default=0, description="The start offset of the pattern.")
    end_pattern: Union[str, None] = Field(
        default=None, description="The end pattern of the pattern."
    )
    end_offset: int = Field(default=0, description="The end offset of the pattern.")

    content_pattern: Union[str, None] = Field(
        default=None, description="The content pattern of the pattern."
    )
    content_repeat: int = Field(
        default=1,
        description="Whether the content pattern should be repeated. -1 means unlimited. "
        "> 0 means the number of times the content pattern should be repeated.",
    )
    description: str = Field(default="", description="The description of the pattern.")

    _start_pattern_compiled: Union[re.Pattern, None] = PrivateAttr(None)
    _end_pattern_compiled: Union[re.Pattern, None] = PrivateAttr(None)
    _content_pattern_compiled: Union[re.Pattern, None] = PrivateAttr(None)

    @model_validator(mode="after")
    def validate_pattern(self) -> Self:
        if self.start_pattern:
            self._start_pattern_compiled = re.compile(self.start_pattern)
        if self.end_pattern:
            self._end_pattern_compiled = re.compile(self.end_pattern)
        if self.content_pattern:
            self._content_pattern_compiled = re.compile(self.content_pattern)
        return self

    @property
    def start_pattern_compiled(self) -> Union[re.Pattern, None]:
        return self._start_pattern_compiled

    @property
    def end_pattern_compiled(self) -> Union[re.Pattern, None]:
        return self._end_pattern_compiled

    @property
    def content_pattern_compiled(self) -> Union[re.Pattern, None]:
        return self._content_pattern_compiled
