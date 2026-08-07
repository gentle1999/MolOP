"""Low-dependency runtime protocols shared across the IO stack."""

from __future__ import annotations

from typing import Any, Protocol, runtime_checkable


@runtime_checkable
class DiskFileLike(Protocol):
    """Minimum file-model surface required by a disk batch."""

    file_path: str

    @property
    def filename(self) -> str: ...

    @property
    def file_format(self) -> str: ...

    def __len__(self) -> int: ...

    def __getitem__(self, index: int) -> object: ...

    def format_transform(self, *args: Any, **kwargs: Any) -> object: ...

    def to_summary_series(self, **kwargs: Any) -> object: ...

    def release_file_content(self) -> None: ...


__all__ = ["DiskFileLike"]
