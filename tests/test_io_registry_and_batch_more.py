from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, cast

import pytest

from molop.io.codec_exceptions import UnsupportedFormatError
from molop.io.codec_registry import Registry
from molop.io.codec_types import ParseResult, StructureLevel
from molop.io.FileBatchModelDisk import FileBatchModelDisk


@dataclass
class DummyReader:
    format_id: str
    extensions: frozenset[str]
    priority: int

    def read(self, path: str | Path, **_kwargs: Any) -> ParseResult[object]:
        return ParseResult(
            value={"path": str(path)},
            level=StructureLevel.COORDS,
            detected_format=self.format_id,
        )


class FakeDiskFile:
    def __init__(
        self, file_path: str, file_format: str, detected_format_id: str | None | object
    ) -> None:
        self.file_path = file_path
        self.filename = Path(file_path).name
        self.file_format = file_format
        if detected_format_id is not _MISSING:
            self.detected_format_id = detected_format_id

    def __len__(self) -> int:
        return 1

    def __getitem__(self, idx: int) -> object:
        return object()

    def format_transform(self, *_args: Any, **_kwargs: Any) -> str:
        return ""

    def to_summary_series(self, **_kwargs: Any) -> object:
        return {}

    def release_file_content(self) -> None:
        return None


_MISSING = object()


def test_registry_select_reader_raises_when_empty_and_autoload_disabled() -> None:
    reg = Registry(autoload_defaults=False)

    with pytest.raises(UnsupportedFormatError, match="No reader codecs registered"):
        reg.select_reader("input.xyz")


def test_registry_write_raises_when_empty_and_autoload_disabled() -> None:
    reg = Registry(autoload_defaults=False)

    with pytest.raises(UnsupportedFormatError, match="No writer codecs registered"):
        reg.write("xyz", value={})


def test_registry_normalizes_format_id_and_extensions_via_public_api() -> None:
    reg = Registry(autoload_defaults=False)
    reg.register_reader_factory(
        lambda: DummyReader(format_id="xyz", extensions=frozenset({".xyz"}), priority=2),
        format_id="  XyZ  ",
        extensions={"xyz", " .XYZ ", ""},
        priority=2,
    )

    selected = reg.select_reader("sample.XYZ", hint_format="  xyz  ")
    assert len(selected) == 1
    assert selected[0].format_id == "xyz"
    assert selected[0].extensions == frozenset({".xyz"})


def test_filebatch_filter_by_codec_id_missing_metadata_keep_drop_and_error() -> None:
    batch = cast(
        Any,
        FileBatchModelDisk(
            cast(
                Any,
                [
                    FakeDiskFile("/tmp/a.xyz", "xyz", "xyz"),
                    FakeDiskFile("/tmp/b.log", "log", "g16log"),
                    FakeDiskFile("/tmp/c.unknown", "unknown", _MISSING),
                ],
            )
        ),
    )

    dropped = batch.filter_by_codec_id("xyz", on_missing="drop")
    assert dropped.file_paths == ["/tmp/a.xyz"]

    kept = batch.filter_by_codec_id(" xyz ", on_missing="keep")
    assert kept.file_paths == ["/tmp/a.xyz", "/tmp/c.unknown"]

    with pytest.raises(ValueError, match="Missing detected_format_id"):
        batch.filter_by_codec_id("xyz", on_missing="error")


def test_filebatch_add_and_slice_semantics_with_fake_diskfiles() -> None:
    left = cast(
        Any,
        FileBatchModelDisk(
            cast(
                Any,
                [
                    FakeDiskFile("/tmp/b.xyz", "xyz", "xyz"),
                    FakeDiskFile("/tmp/a.xyz", "xyz", "xyz"),
                ],
            )
        ),
    )
    right = cast(
        Any,
        FileBatchModelDisk(
            cast(
                Any,
                [
                    FakeDiskFile("/tmp/b.xyz", "xyz", "xyz"),
                    FakeDiskFile("/tmp/c.xyz", "xyz", "xyz"),
                ],
            )
        ),
    )

    merged = left + right
    assert merged.file_paths == ["/tmp/a.xyz", "/tmp/b.xyz", "/tmp/c.xyz"]

    sliced = merged[1:]
    assert isinstance(sliced, FileBatchModelDisk)
    assert sliced.file_paths == ["/tmp/b.xyz", "/tmp/c.xyz"]
