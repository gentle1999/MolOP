"""
Author: TMJ
Date: 2025-07-29 16:32:20
LastEditors: TMJ
LastEditTime: 2025-11-26 13:47:31
Description: 请填写简介
"""

from collections.abc import Generator, Iterator
from typing import Any

import regex
from pydantic import BaseModel, ConfigDict, Field, PrivateAttr, model_validator
from typing_extensions import Self


def find_iter_no_regex(
    pattern: str, content: str, start_idx: int = 0, end_idx: int | None = None
) -> Generator[tuple[int, int], None, None]:
    if end_idx is None:
        end_idx = len(content)
    current_idx = start_idx
    pattern_len = len(pattern)

    while True:
        match_start = content.find(pattern, current_idx, end_idx)
        if match_start == -1:
            break
        match_end = match_start + pattern_len
        yield match_start, match_end
        current_idx = match_start + pattern_len


def find_no_regex(
    pattern: str, content: str, start_idx: int = 0, end_idx: int | None = None
) -> tuple[int, int] | None:
    if end_idx is None:
        end_idx = len(content)
    match_start = content.find(pattern, start_idx, end_idx)
    if match_start == -1:
        return None
    match_end = match_start + len(pattern)
    return match_start, match_end


class MolOPPattern(BaseModel):
    model_config = ConfigDict(arbitrary_types_allowed=True)
    start_pattern: str | None = Field(default=None, description="The start pattern of the pattern.")
    start_offset: int = Field(default=0, ge=0, description="The start offset of the pattern.")
    start_regex: bool = Field(
        default=True,
        description="Whether the start pattern should be treated as a regex pattern.",
    )
    end_pattern: str | None = Field(default=None, description="The end pattern of the pattern.")
    end_offset: int = Field(default=0, ge=0, description="The end offset of the pattern.")
    end_regex: bool = Field(
        default=True,
        description="Whether the end pattern should be treated as a regex pattern.",
    )
    content_pattern: str | None = Field(
        default=None, description="The content pattern of the pattern."
    )
    content_repeat: int = Field(
        default=1,
        description="Whether the content pattern should be repeated. 0 means unlimited. "
        "> 0 means the number of times the content pattern should be repeated.",
    )
    description: str = Field(default="", description="The description of the pattern.")

    _start_pattern_compiled: regex.Pattern[str] | None = PrivateAttr(None)
    _end_pattern_compiled: regex.Pattern[str] | None = PrivateAttr(None)
    _content_pattern_compiled: regex.Pattern[str] | None = PrivateAttr(None)

    @model_validator(mode="after")
    def validate_pattern(self) -> Self:
        if self.start_pattern and self.start_regex:
            self._start_pattern_compiled = regex.compile(self.start_pattern, regex.MULTILINE)
        if self.end_pattern and self.end_regex:
            self._end_pattern_compiled = regex.compile(self.end_pattern, regex.MULTILINE)
        if self.content_pattern:
            self._content_pattern_compiled = regex.compile(self.content_pattern, regex.MULTILINE)
        return self

    @property
    def start_pattern_compiled(self) -> regex.Pattern[str] | None:
        return self._start_pattern_compiled

    @property
    def end_pattern_compiled(self) -> regex.Pattern[str] | None:
        return self._end_pattern_compiled

    @property
    def content_pattern_compiled(self) -> regex.Pattern[str] | None:
        return self._content_pattern_compiled

    def locate_content(self, content: str) -> None | tuple[int, int, int, int]:
        """
        Locate the content of the pattern in the given content.

        Parameters:
            content (str): The content to be searched.

        Returns:
            (None | tuple[int, int, int, int]): The start, end pos of start_pattern and start, end pos of end_pattern.
        """
        if self.start_pattern_compiled:
            for idx, match in enumerate(self.start_pattern_compiled.finditer(content)):
                if idx >= self.start_offset:
                    start_index, start_pos = match.start(), match.end()
                    break
            else:
                return None
        elif self.start_pattern:
            for idx, match in enumerate(find_iter_no_regex(self.start_pattern, content)):
                if idx >= self.start_offset:
                    start_index, start_pos = match[0], match[1]
                    break
            else:
                return None
        else:
            start_index, start_pos = 0, 0
        if self.end_pattern_compiled:
            for idx, match in enumerate(
                self.end_pattern_compiled.finditer(content, pos=start_index)
            ):
                if idx >= self.end_offset:
                    end_index, end_pos = match.start(), match.end()
                    break
            else:
                return None
        elif self.end_pattern:
            for idx, match in enumerate(find_iter_no_regex(self.end_pattern, content, start_index)):
                if idx >= self.end_offset:
                    end_index, end_pos = match[0], match[1]
                    break
            else:
                return None
        else:
            end_index, end_pos = start_pos, len(content)
        assert end_index >= start_index, (
            f"end_index should be greater than or equal to start_index, but got {end_index} < {start_index}"
        )
        assert end_pos >= start_pos, (
            f"end_pos should be greater than or equal to start_pos, but got {end_pos} < {start_pos}"
        )
        return start_index, start_pos, end_index, end_pos

    def match_content(self, content: str) -> None | list[tuple[str | Any, ...]]:
        """
        Match the content of the pattern in the given content.

        Parameters:
            content (str): The content to be searched.

        Returns:
            (None | list[tuple[str | Any, ...]]): The matched content.
        """
        if located_content_index := self.locate_content(content):
            start_start, start_end, end_start, end_end = located_content_index
            located_content = content[start_start:end_end]
            total_matches = self.get_matches(located_content)
            return total_matches
        return None

    def get_matches(self, located_content: str) -> None | list[tuple[str | Any, ...]]:
        if not self.content_pattern_compiled:
            return None
        total_matches: list[tuple[str | Any, ...]] = []
        if self.content_repeat == 0:
            total_matches = [
                match.groups() for match in self.content_pattern_compiled.finditer(located_content)
            ]
        if self.content_repeat > 0:
            while (len(total_matches) < self.content_repeat) and (
                match := self.content_pattern_compiled.search(located_content)
            ):
                total_matches.append(match.groups())
                located_content = located_content[match.end() :]
        if self.content_repeat < 0:
            while (len(total_matches) < abs(self.content_repeat)) and (
                match := self.content_pattern_compiled.search(located_content)
            ):
                total_matches.append(match.groups())
                located_content = located_content[: match.start()]
        return total_matches

    def find_iter(self, located_content: str) -> None | Iterator[regex.Match[str]]:
        if not self.content_pattern_compiled:
            return None
        if self.content_repeat == 0:
            return self.content_pattern_compiled.finditer(located_content)
        else:
            return (
                match
                for idx, match in enumerate(self.content_pattern_compiled.finditer(located_content))
                if idx < abs(self.content_repeat)
            )

    def split_content(self, content: str) -> tuple[str, str]:
        """
        Split the content of the pattern in the given content.

        Parameters:
            content (str): The content to be searched.

        Returns:
            (tuple[str, str]): The start and end content.
        """
        if located_content_index := self.locate_content(content):
            start_start, start_end, end_start, end_end = located_content_index
            return content[start_start:end_end], content[end_start:]
        return "", content


class MolOPPatternV2(MolOPPattern):
    @classmethod
    def from_pattern(cls, pattern: MolOPPattern) -> "MolOPPatternV2":
        return cls.model_validate(pattern.model_dump())

    def _locate_nth_start(self, content: str, start_pos: int) -> tuple[int, int] | None:
        if self.start_pattern_compiled:
            current_pos = start_pos
            start_index = start_end = current_pos
            for _ in range(self.start_offset + 1):
                match = self.start_pattern_compiled.search(content, pos=current_pos)
                if match is None:
                    return None
                start_index, start_end = match.start(), match.end()
                current_pos = match.end()
            return start_index, start_end
        if self.start_pattern:
            current_pos = start_pos
            start_index = start_end = current_pos
            for _ in range(self.start_offset + 1):
                match = find_no_regex(self.start_pattern, content, start_idx=current_pos)
                if match is None:
                    return None
                start_index, start_end = match
                current_pos = start_end
            return start_index, start_end
        return start_pos, start_pos

    def _locate_nth_end(self, content: str, start_index: int) -> tuple[int, int] | None:
        if self.end_pattern_compiled:
            current_pos = start_index
            end_index = end_end = current_pos
            for _ in range(self.end_offset + 1):
                match = self.end_pattern_compiled.search(content, pos=current_pos)
                if match is None:
                    return None
                end_index, end_end = match.start(), match.end()
                current_pos = match.end()
            return end_index, end_end
        if self.end_pattern:
            current_pos = start_index
            end_index = end_end = current_pos
            for _ in range(self.end_offset + 1):
                match = find_no_regex(self.end_pattern, content, start_idx=current_pos)
                if match is None:
                    return None
                end_index, end_end = match
                current_pos = end_end
            return end_index, end_end
        return start_index, len(content)

    def locate_content_from(
        self, content: str, start_pos: int = 0
    ) -> tuple[int, int, int, int] | None:
        if (start_match := self._locate_nth_start(content, start_pos)) is None:
            return None
        start_index, start_end = start_match
        if (end_match := self._locate_nth_end(content, start_index)) is None:
            return None
        end_index, end_end = end_match
        assert end_index >= start_index, (
            f"end_index should be greater than or equal to start_index, but got {end_index} < {start_index}"
        )
        assert end_end >= start_end, (
            f"end_pos should be greater than or equal to start_pos, but got {end_end} < {start_end}"
        )
        return start_index, start_end, end_index, end_end

    def match_content_from(
        self, content: str, start_pos: int = 0
    ) -> None | list[tuple[str | Any, ...]]:
        if located_content_index := self.locate_content_from(content, start_pos):
            start_start, _start_end, _end_start, end_end = located_content_index
            located_content = content[start_start:end_end]
            return self.get_matches(located_content)
        return None

    def split_content_from(self, content: str, start_pos: int = 0) -> tuple[str, int]:
        if located_content_index := self.locate_content_from(content, start_pos):
            start_start, _start_end, end_start, end_end = located_content_index
            return content[start_start:end_end], end_start
        return "", start_pos
