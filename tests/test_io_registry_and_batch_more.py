from __future__ import annotations

import importlib
from dataclasses import dataclass
from pathlib import Path
from typing import Any, cast

import pytest

from molop.io.base_models._format_transform import FrameFormatTransformMixin
from molop.io.codec_exceptions import ConversionError
from molop.io.codec_exceptions import UnsupportedFormatError
from molop.io.codec_registry import Registry
from molop.io.codec_types import ParseResult, StructureLevel
from molop.io.FileBatchModelDisk import FileBatchModelDisk
from molop.io.FileBatchParserDisk import FileBatchParserDisk

filebatchparserdisk_module = importlib.import_module("molop.io.FileBatchParserDisk")


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


class DummyFrameTransform(FrameFormatTransformMixin):
    def __init__(self) -> None:
        self.frame_id = 0
        self.charge = 0
        self.multiplicity = 1
        self.rdmol = object()
        self.omol = type(
            "DummyOMol", (), {"write": lambda self, fmt: "<cml />" if fmt == "cml" else None}
        )()

    def model_dump(self, **_kwargs: Any) -> dict[str, Any]:
        return {
            "atoms": [1],
            "coords": [[0.0, 0.0, 0.0]],
            "charge": self.charge,
            "multiplicity": self.multiplicity,
        }

    def _render(self, **_kwargs: Any) -> str:
        return "dummy-frame"


class DummyStructuredFile:
    def __init__(self) -> None:
        self.frames = [object()]

    def model_dump(self, **_kwargs: Any) -> dict[str, Any]:
        return {"frames": 1}


@dataclass
class DummyWriter:
    format_id: str
    required_level: StructureLevel
    priority: int = 100

    def write(self, value: object, **_kwargs: Any) -> str:
        return f"{self.format_id}:{type(value).__name__}"


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


def test_registry_uses_writer_default_graph_policy_when_unspecified() -> None:
    reg = Registry(autoload_defaults=False)

    class GraphWriter:
        format_id = "dual"
        required_level = StructureLevel.GRAPH
        priority = 100

        def write(self, value: object, **_kwargs: Any) -> str:
            return f"graph:{type(value).__name__}"

    class CoordsWriter:
        format_id = "dual"
        required_level = StructureLevel.COORDS
        priority = 10

        def write(self, value: object, **_kwargs: Any) -> str:
            return f"coords:{type(value).__name__}"

    reg.register_writer_factory(
        lambda: GraphWriter(),
        format_id="dual",
        required_level=StructureLevel.GRAPH,
        domain="file",
        default_graph_policy="prefer",
        priority=100,
    )
    reg.register_writer_factory(
        lambda: CoordsWriter(),
        format_id="dual",
        required_level=StructureLevel.COORDS,
        domain="file",
        default_graph_policy="coords",
        priority=10,
    )

    rendered = reg.write("dual", DummyStructuredFile())

    assert rendered == "graph:DummyStructuredFile"


def test_registry_allows_manual_graph_policy_override() -> None:
    reg = Registry(autoload_defaults=False)

    class GraphWriter:
        format_id = "dual"
        required_level = StructureLevel.GRAPH
        priority = 100

        def write(self, value: object, **_kwargs: Any) -> str:
            return f"graph:{type(value).__name__}"

    class CoordsWriter:
        format_id = "dual"
        required_level = StructureLevel.COORDS
        priority = 10

        def write(self, value: object, **_kwargs: Any) -> str:
            return f"coords:{type(value).__name__}"

    reg.register_writer_factory(
        lambda: GraphWriter(),
        format_id="dual",
        required_level=StructureLevel.GRAPH,
        domain="frame",
        default_graph_policy="prefer",
        priority=100,
    )
    reg.register_writer_factory(
        lambda: CoordsWriter(),
        format_id="dual",
        required_level=StructureLevel.COORDS,
        domain="frame",
        default_graph_policy="coords",
        priority=10,
    )

    frame_value = DummyFrameTransform()

    assert reg.write_frame("dual", frame_value) == "graph:DummyFrameTransform"
    assert (
        reg.write_frame("dual", frame_value, graph_policy="coords") == "coords:DummyFrameTransform"
    )


def test_registry_raises_when_coords_override_has_no_coords_writer() -> None:
    reg = Registry(autoload_defaults=False)

    class GraphWriter:
        format_id = "sdf"
        required_level = StructureLevel.GRAPH
        priority = 100

        def write(self, value: object, **_kwargs: Any) -> str:
            return f"graph:{type(value).__name__}"

    reg.register_writer_factory(
        lambda: GraphWriter(),
        format_id="sdf",
        required_level=StructureLevel.GRAPH,
        domain="file",
        default_graph_policy="prefer",
        priority=100,
    )

    with pytest.raises(ConversionError, match="No coords-level writers registered"):
        reg.write("sdf", DummyStructuredFile(), graph_policy="coords")


def test_registry_builtin_default_graph_policies_match_format_semantics() -> None:
    reg = Registry(autoload_defaults=False)

    reg.register_writer_factory(
        lambda: DummyWriter("sdf", StructureLevel.GRAPH),
        format_id="sdf",
        required_level=StructureLevel.GRAPH,
        domain="file",
        default_graph_policy="strict",
        priority=100,
    )
    reg.register_writer_factory(
        lambda: DummyWriter("smi", StructureLevel.GRAPH),
        format_id="smi",
        required_level=StructureLevel.GRAPH,
        domain="file",
        default_graph_policy="strict",
        priority=100,
    )
    reg.register_writer_factory(
        lambda: DummyWriter("gjf", StructureLevel.COORDS),
        format_id="gjf",
        required_level=StructureLevel.COORDS,
        domain="file",
        default_graph_policy="prefer",
        priority=100,
    )

    assert reg._writers_by_format["sdf"][0].default_graph_policy == "strict"
    assert reg._writers_by_format["smi"][0].default_graph_policy == "strict"
    assert reg._writers_by_format["gjf"][0].default_graph_policy == "prefer"


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


def test_filebatch_iteration_is_reentrant_for_nested_loops() -> None:
    batch = cast(
        Any,
        FileBatchModelDisk(
            cast(
                Any,
                [
                    FakeDiskFile("/tmp/a.xyz", "xyz", "xyz"),
                    FakeDiskFile("/tmp/b.xyz", "xyz", "xyz"),
                ],
            )
        ),
    )

    pairs = [(left.file_path, right.file_path) for left in batch for right in batch]

    assert pairs == [
        ("/tmp/a.xyz", "/tmp/a.xyz"),
        ("/tmp/a.xyz", "/tmp/b.xyz"),
        ("/tmp/b.xyz", "/tmp/a.xyz"),
        ("/tmp/b.xyz", "/tmp/b.xyz"),
    ]


def test_filter_custom_keeps_alignment_with_generator_results(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    batch = cast(
        Any,
        FileBatchModelDisk(
            cast(
                Any,
                [
                    FakeDiskFile("/tmp/a.xyz", "xyz", "xyz"),
                    FakeDiskFile("/tmp/b.xyz", "xyz", "xyz"),
                    FakeDiskFile("/tmp/c.xyz", "xyz", "xyz"),
                ],
            )
        ),
    )

    def generator_parallel_execute(
        _func: Any, _desc: str = "", _n_jobs: int = 1, **_kwargs: Any
    ) -> Any:
        return (keep for keep in [True, False, True])

    monkeypatch.setattr(batch, "parallel_execute", generator_parallel_execute)

    filtered = batch.filter_custom(lambda _diskfile: True)

    assert filtered.file_paths == ["/tmp/a.xyz", "/tmp/c.xyz"]


def test_groupby_uses_generator_results_without_index_skew(monkeypatch: pytest.MonkeyPatch) -> None:
    batch = cast(
        Any,
        FileBatchModelDisk(
            cast(
                Any,
                [
                    FakeDiskFile("/tmp/a.xyz", "xyz", "xyz"),
                    FakeDiskFile("/tmp/b.log", "log", "g16log"),
                    FakeDiskFile("/tmp/c.xyz", "xyz", "xyz"),
                ],
            )
        ),
    )

    def generator_parallel_execute(
        _func: Any, _desc: str = "", _n_jobs: int = 1, **_kwargs: Any
    ) -> Any:
        return (key for key in ["xyz", "log", "xyz"])

    monkeypatch.setattr(batch, "parallel_execute", generator_parallel_execute)

    grouped = batch.groupby(lambda diskfile: diskfile.file_format)

    assert grouped["xyz"].file_paths == ["/tmp/a.xyz", "/tmp/c.xyz"]
    assert grouped["log"].file_paths == ["/tmp/b.log"]


def test_filebatchparser_parallel_generator_path_filters_results_safely(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    parser = FileBatchParserDisk(n_jobs=2)
    file_paths = ["/tmp/b.xyz", "/tmp/a.xyz", "/tmp/c.xyz"]

    monkeypatch.setattr(filebatchparserdisk_module.os.path, "isfile", lambda path: True)
    monkeypatch.setattr(filebatchparserdisk_module.os.path, "abspath", lambda path: cast(str, path))
    monkeypatch.setattr(
        filebatchparserdisk_module.os.path,
        "getsize",
        lambda path: {"/tmp/b.xyz": 30, "/tmp/a.xyz": 20, "/tmp/c.xyz": 10}[cast(str, path)],
    )
    monkeypatch.setattr(
        filebatchparserdisk_module.codec_registry,
        "select_reader",
        lambda _path, hint_format=None: (DummyReader("xyz", frozenset({".xyz"}), 1),),
    )

    parsed_by_path = {
        "/tmp/b.xyz": FakeDiskFile("/tmp/b.xyz", "xyz", "xyz"),
        "/tmp/a.xyz": None,
        "/tmp/c.xyz": FakeDiskFile("/tmp/c.xyz", "xyz", "xyz"),
    }

    monkeypatch.setattr(
        filebatchparserdisk_module,
        "single_file_parser",
        lambda **task: parsed_by_path[task["file_path"]],
    )

    monkeypatch.setattr(
        filebatchparserdisk_module, "delayed", lambda func: lambda **task: lambda: func(**task)
    )

    class StubParallel:
        def __init__(self, **_kwargs: Any) -> None:
            pass

        def __call__(self, iterable: Any) -> Any:
            return (callable_obj() for callable_obj in iterable)

    monkeypatch.setattr(filebatchparserdisk_module, "Parallel", StubParallel)

    batch = parser.parse(file_paths)

    assert batch.file_paths == ["/tmp/b.xyz", "/tmp/c.xyz"]


def test_frame_format_transform_routes_single_frame_through_codec_registry(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    frame = DummyFrameTransform()
    captured: dict[str, Any] = {}

    def fake_write(format_id: str, value: object, **kwargs: Any) -> str:
        captured["format_id"] = format_id
        captured["value"] = value
        captured["kwargs"] = kwargs
        return "rendered-block"

    monkeypatch.setattr(
        "molop.io.base_models._format_transform.codec_registry.write_frame", fake_write
    )

    rendered = frame.format_transform("xyz", graph_policy="coords")

    assert rendered == "rendered-block"
    assert captured["format_id"] == "xyz"
    assert captured["kwargs"]["graph_policy"] == "coords"
    assert captured["value"] is frame


def test_frame_format_transform_defers_default_graph_policy_to_registry(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    frame = DummyFrameTransform()
    captured: dict[str, Any] = {}

    def fake_write(format_id: str, value: object, **kwargs: Any) -> str:
        captured["format_id"] = format_id
        captured["value"] = value
        captured["kwargs"] = kwargs
        return "rendered-block"

    monkeypatch.setattr(
        "molop.io.base_models._format_transform.codec_registry.write_frame", fake_write
    )

    rendered = frame.format_transform("xyz")

    assert rendered == "rendered-block"
    assert captured["kwargs"]["graph_policy"] is None


def test_frame_format_transform_writes_requested_output_file(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    frame = DummyFrameTransform()
    monkeypatch.setattr(
        "molop.io.base_models._format_transform.codec_registry.write_frame",
        lambda *_args, **_kwargs: "xyz-block",
    )

    target = tmp_path / "frame.anything"
    rendered = frame.format_transform("xyz", file_path=target)

    assert rendered == "xyz-block"
    assert (tmp_path / "frame.xyz").read_text() == "xyz-block"


def test_frame_format_transform_cml_routes_through_codec_registry(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    frame = DummyFrameTransform()
    captured: dict[str, Any] = {}

    def fake_write(format_id: str, value: object, **kwargs: Any) -> str:
        captured["format_id"] = format_id
        captured["value"] = value
        captured["kwargs"] = kwargs
        return "<cml />"

    monkeypatch.setattr(
        "molop.io.base_models._format_transform.codec_registry.write_frame", fake_write
    )

    assert frame.format_transform("cml", engine="openbabel") == "<cml />"
    assert captured["format_id"] == "cml"
    assert captured["value"] is frame
