"""Common comment containers shared by file and frame models."""

from __future__ import annotations

from collections.abc import Iterable, Iterator
from typing import Any, Literal, overload

from pydantic import ConfigDict, Field

from ..Bases import BaseDataClassWithUnit


CommentKind = Literal["comment", "title", "annotation"]


class Comment(BaseDataClassWithUnit):
    """A format-neutral textual annotation.

    ``kind`` describes the semantic role when a format distinguishes a title
    card from an ordinary comment.  ``marker`` and ``raw`` retain enough
    format evidence for callers that need to inspect or reproduce the source
    syntax without making it part of the common access API.
    """

    model_config = ConfigDict(extra="forbid")

    text: str = Field(default="", description="Normalized comment text")
    kind: CommentKind = Field(default="comment", description="Semantic comment kind")
    marker: str | None = Field(
        default=None,
        description="Source comment marker such as '#' or '!'; absent for markerless formats",
    )
    source_format: str | None = Field(
        default=None,
        description="Format that produced this comment, when known",
    )
    source_line: int | None = Field(
        default=None,
        ge=0,
        description="Zero-based source line, when known",
    )
    raw: str | None = Field(
        default=None,
        description="Original source spelling, when it differs from ``text``",
    )
    metadata: dict[str, Any] = Field(
        default_factory=dict,
        description="Format-specific comment metadata that has no common field",
    )


class CommentContainer(BaseDataClassWithUnit):
    """Ordered, list-like storage for file-level or frame-level comments."""

    model_config = ConfigDict(extra="forbid")

    items: list[Comment] = Field(default_factory=list, description="Ordered comments")

    def __iter__(self) -> Iterator[Comment]:  # type: ignore[override]
        return iter(self.items)

    def __len__(self) -> int:
        return len(self.items)

    @overload
    def __getitem__(self, index: int) -> Comment: ...

    @overload
    def __getitem__(self, index: slice) -> list[Comment]: ...

    def __getitem__(self, index: int | slice) -> Comment | list[Comment]:
        return self.items[index]

    def add(
        self,
        comment: Comment | str,
        *,
        kind: CommentKind = "comment",
        marker: str | None = None,
        source_format: str | None = None,
        source_line: int | None = None,
        raw: str | None = None,
        metadata: dict[str, Any] | None = None,
    ) -> Comment:
        """Append one comment and return the stored model."""

        stored = (
            comment
            if isinstance(comment, Comment)
            else Comment(
                text=comment,
                kind=kind,
                marker=marker,
                source_format=source_format,
                source_line=source_line,
                raw=raw,
                metadata={} if metadata is None else metadata,
            )
        )
        self.items.append(stored)
        return stored

    def ensure(self, comment: Comment) -> Comment:
        """Append ``comment`` unless an identical item is already present."""

        if comment not in self.items:
            self.items.append(comment)
        return next(existing for existing in self.items if existing == comment)

    def extend(self, comments: Iterable[Comment | str]) -> None:
        """Append comments in order."""

        for comment in comments:
            self.add(comment)

    def first(self, *, kind: CommentKind | None = None) -> Comment | None:
        """Return the first comment, optionally restricted to one kind."""

        return next(
            (comment for comment in self.items if kind is None or comment.kind == kind),
            None,
        )

    def first_text(self, *, kind: CommentKind | None = None) -> str:
        """Return the first matching comment's text, or an empty string."""

        comment = self.first(kind=kind)
        return comment.text if comment is not None else ""

    @property
    def texts(self) -> list[str]:
        """Return comment text in source order."""

        return [comment.text for comment in self.items]


__all__ = ["Comment", "CommentContainer", "CommentKind"]
