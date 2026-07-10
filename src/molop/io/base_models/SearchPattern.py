"""
Author: TMJ
Date: 2025-07-29 16:32:20
LastEditors: TMJ
LastEditTime: 2025-11-26 13:47:31
Description: 请填写简介
"""

from collections.abc import Generator, Iterator
from typing import TypeAlias

import regex
from pydantic import BaseModel, ConfigDict, Field, PrivateAttr, model_validator
from typing_extensions import Self


MolOPMatch: TypeAlias = regex.Match[str]


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
        "Non-zero values cap the number of matches by abs(value); regex flags such "
        "as (?r) control match direction.",
    )
    description: str = Field(default="", description="The description of the pattern.")

    _start_pattern_compiled: regex.Pattern[str] | None = PrivateAttr(None)
    _end_pattern_compiled: regex.Pattern[str] | None = PrivateAttr(None)
    _content_pattern_compiled: regex.Pattern[str] | None = PrivateAttr(None)

    @staticmethod
    def escape_literal(text: str) -> str:
        return regex.escape(text)

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

    @property
    def group_names(self) -> tuple[str, ...]:
        if not self.content_pattern_compiled:
            return ()
        return tuple(self.content_pattern_compiled.groupindex)

    def has_named_group(self, group_name: str) -> bool:
        return group_name in self.group_names

    def _ensure_named_group(self, group_name: str) -> None:
        if not self.has_named_group(group_name):
            raise IndexError(f"no such group: {group_name}")

    def get_named_group(self, matched: MolOPMatch, group_name: str) -> str | None:
        self._ensure_named_group(group_name)
        return matched.group(group_name)

    def require_named_group(self, matched: MolOPMatch, group_name: str) -> str:
        value = self.get_named_group(matched, group_name)
        if value is None:
            raise ValueError(f"group did not match: {group_name}")
        return value

    def named_group_dict(self, matched: MolOPMatch) -> dict[str, str | None]:
        return {name: matched.group(name) for name in self.group_names}

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
                    start_index, start_pos = match
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
                    end_index, end_pos = match
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

    def match_content(self, content: str) -> None | list[MolOPMatch]:
        """
        Match the content of the pattern in the given content.

        Parameters:
            content (str): The content to be searched.

        Returns:
            (None | list[MolOPMatch]): The matched content.
        """
        if located_content_index := self.locate_content(content):
            start_start, start_end, end_start, end_end = located_content_index
            located_content = content[start_start:end_end]
            return self.get_matches(located_content)
        return None

    def get_matches(self, located_content: str) -> None | list[MolOPMatch]:
        if not self.content_pattern_compiled:
            return None
        return self.find_matches(located_content)

    def find_matches(self, located_content: str) -> list[MolOPMatch]:
        return list(self.find_iter(located_content) or [])

    def find_named_group_values(self, located_content: str, group_name: str) -> list[str]:
        self._ensure_named_group(group_name)
        values: list[str] = []
        for matched in self.find_matches(located_content):
            value = matched.group(group_name)
            if value is not None:
                values.append(value)
        return values

    def find_first_named_group_value(self, located_content: str, group_name: str) -> str | None:
        self._ensure_named_group(group_name)
        for matched in self.find_matches(located_content):
            value = matched.group(group_name)
            if value is not None:
                return value
        return None

    def find_named_group_dicts(self, located_content: str) -> list[dict[str, str | None]]:
        return [self.named_group_dict(matched) for matched in self.find_matches(located_content)]

    def find_group_values(self, located_content: str, group_name: str) -> list[str]:
        return self.find_named_group_values(located_content, group_name)

    def find_first_group_value(self, located_content: str, group_name: str) -> str | None:
        return self.find_first_named_group_value(located_content, group_name)

    def get_group(self, matched: MolOPMatch, group_name: str) -> str | None:
        return self.get_named_group(matched, group_name)

    def require_group(self, matched: MolOPMatch, group_name: str) -> str:
        return self.require_named_group(matched, group_name)

    def group_dict(self, matched: MolOPMatch) -> dict[str, str | None]:
        return self.named_group_dict(matched)

    def find_groupdicts(self, located_content: str) -> list[dict[str, str | None]]:
        return self.find_named_group_dicts(located_content)

    def find_group_dicts(self, located_content: str) -> list[dict[str, str | None]]:
        return self.find_named_group_dicts(located_content)

    def search(self, content: str, pos: int = 0) -> MolOPMatch | None:
        if not self.content_pattern_compiled:
            return None
        return self.content_pattern_compiled.search(content, pos=pos)

    def match(self, content: str, pos: int = 0) -> MolOPMatch | None:
        if not self.content_pattern_compiled:
            return None
        return self.content_pattern_compiled.match(content, pos=pos)

    def find_iter(self, located_content: str) -> None | Iterator[MolOPMatch]:
        if not self.content_pattern_compiled:
            return None
        if self.content_repeat == 0:
            return self.content_pattern_compiled.finditer(located_content)
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

    def match_content_from(self, content: str, start_pos: int = 0) -> None | list[MolOPMatch]:
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

    def find_content_matches(self, content: str, start_pos: int = 0) -> None | list[MolOPMatch]:
        if located_content_index := self.locate_content_from(content, start_pos):
            start_start, _start_end, _end_start, end_end = located_content_index
            return self.find_matches(content[start_start:end_end])
        return None

    @classmethod
    def from_pattern(cls, pattern: "MolOPPattern") -> Self:
        return cls.model_validate(pattern.model_dump())


# Compatibility alias only. Keep MolOPPattern as the single implementation.
MolOPPatternV2: TypeAlias = MolOPPattern
