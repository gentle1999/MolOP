from __future__ import annotations

import importlib
from dataclasses import dataclass
from pathlib import Path
from typing import Any, cast

import pandas as pd
import pytest

from molop.io.base_models._format_transform import FrameFormatTransformMixin
from molop.io.codec_exceptions import ConversionError, FormatMismatchError, UnsupportedFormatError
from molop.io.codec_registry import Registry
from molop.io.codec_types import ParseResult, ParseWarning, StructureLevel
from molop.io.FileBatchModelDisk import FileBatchModelDisk
from molop.io.FileBatchParserDisk import FileBatchParserDisk
from molop.io.parse_outcomes import BatchParseResult, FileParseOutcome
from molop.unit import atom_ureg


io_module = importlib.import_module("molop.io")
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


@dataclass
class ProbingReader:
    format_id: str
    probe_result: bool

    def probe_file_format(self, path: str | Path) -> bool:
        _ = path
        return self.probe_result

    def read(self, path: str | Path, **_kwargs: Any) -> ParseResult[object]:
        return ParseResult(
            value=FakeDiskFile(str(path), Path(path).suffix.lstrip("."), self.format_id),
            level=StructureLevel.COORDS,
            detected_format=self.format_id,
        )


@dataclass
class FailingReader:
    format_id: str
    error: Exception

    def read(self, path: str | Path, **_kwargs: Any) -> ParseResult[object]:
        raise self.error


@dataclass
class FakeDiskFileReader:
    format_id: str

    def read(self, path: str | Path, **_kwargs: Any) -> ParseResult[object]:
        return ParseResult(
            value=FakeDiskFile(str(path), Path(path).suffix.lstrip("."), self.format_id),
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


class FakeFrame:
    def __init__(self, num_imaginary: int | None = None, num_frequencies: int = 1) -> None:
        self.vibrations = (
            None
            if num_imaginary is None
            else type(
                "FakeVibrations",
                (),
                {"num_imaginary": num_imaginary, "__len__": lambda _self: num_frequencies},
            )()
        )


class FakeStateDiskFile(FakeDiskFile):
    def __init__(self, file_path: str, frames: list[FakeFrame]) -> None:
        super().__init__(file_path, "log", "g16log")
        self.frames = frames

    def __len__(self) -> int:
        return len(self.frames)

    def __iter__(self) -> object:
        return iter(self.frames)

    def __getitem__(self, idx: int) -> FakeFrame:
        return self.frames[idx]


class FakeSummaryFrame:
    def __init__(self, frame_id: int) -> None:
        self.frame_id = frame_id

    def to_summary_series(self, brief: bool = True, **_kwargs: Any) -> pd.Series:
        return pd.Series(
            {
                ("General", "FrameID"): self.frame_id,
                ("General", "Brief"): brief,
            }
        )


class FakeSummaryDiskFile(FakeDiskFile):
    def __init__(self, file_path: str, frame_count: int) -> None:
        super().__init__(file_path, "log", "g16log")
        self.frames = [FakeSummaryFrame(frame_id) for frame_id in range(frame_count)]

    def __len__(self) -> int:
        return len(self.frames)

    def __getitem__(self, idx: int) -> FakeSummaryFrame:
        return self.frames[idx]

    def to_summary_series(self, brief: bool = True, **_kwargs: Any) -> pd.Series:
        return pd.Series(
            {
                ("File", "Path"): self.file_path,
                ("File", "Brief"): brief,
            }
        )


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


def test_filebatch_internal_sorted_constructor_preserves_sorted_order() -> None:
    batch = cast(
        Any,
        FileBatchModelDisk._new_batch_from_sorted_diskfiles(
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

    assert batch.file_paths == ["/tmp/a.xyz", "/tmp/b.xyz", "/tmp/c.xyz"]


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


def test_parallel_execute_suppresses_implicit_none_results_by_default() -> None:
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
    seen: list[str] = []

    def collect_path(diskfile: FakeDiskFile) -> None:
        seen.append(diskfile.file_path)

    result = batch.parallel_execute(collect_path, n_jobs=1)

    assert result is None
    assert seen == ["/tmp/a.xyz", "/tmp/b.xyz"]


def test_parallel_execute_returns_values_by_default() -> None:
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

    result = batch.parallel_execute(lambda diskfile: diskfile.file_path, n_jobs=1)

    assert result == ["/tmp/a.xyz", "/tmp/b.xyz"]


def test_parallel_execute_can_preserve_explicit_none_results() -> None:
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

    result = batch.parallel_execute(lambda _diskfile: None, n_jobs=1, return_results=True)

    assert result == [None, None]


def test_parallel_execute_uses_provided_snapshot(
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
                ],
            )
        ),
    )
    snapshot = batch._snapshot_diskfiles()

    def fail_snapshot() -> list[Any]:
        raise AssertionError("parallel_execute should reuse the provided snapshot")

    monkeypatch.setattr(batch, "_snapshot_diskfiles", fail_snapshot)

    result = batch.parallel_execute(
        lambda diskfile: diskfile.file_path,
        n_jobs=1,
        _diskfiles_snapshot=snapshot,
    )

    assert result == ["/tmp/a.xyz", "/tmp/b.xyz"]


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


def test_filter_custom_reuses_existing_snapshot(monkeypatch: pytest.MonkeyPatch) -> None:
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
    snapshot = batch._snapshot_diskfiles()
    snapshot_calls = 0

    def counted_snapshot() -> list[Any]:
        nonlocal snapshot_calls
        snapshot_calls += 1
        return snapshot

    monkeypatch.setattr(batch, "_snapshot_diskfiles", counted_snapshot)

    filtered = batch.filter_custom(lambda diskfile: diskfile.file_path != "/tmp/b.xyz")

    assert snapshot_calls == 1
    assert filtered.file_paths == ["/tmp/a.xyz", "/tmp/c.xyz"]


def test_to_summary_df_frame_all_returns_every_frame() -> None:
    batch = cast(
        Any,
        FileBatchModelDisk(
            cast(
                Any,
                [
                    FakeSummaryDiskFile("/tmp/a.log", 2),
                    FakeSummaryDiskFile("/tmp/b.log", 1),
                ],
            )
        ),
    )

    df = batch.to_summary_df(frame="all", n_jobs=1)

    assert df[("General", "FrameID", "")].tolist() == [0, 1, 0]


def test_to_summary_df_keeps_last_frame_default() -> None:
    batch = cast(
        Any,
        FileBatchModelDisk(
            cast(
                Any,
                [
                    FakeSummaryDiskFile("/tmp/a.log", 2),
                    FakeSummaryDiskFile("/tmp/b.log", 1),
                ],
            )
        ),
    )

    df = batch.to_summary_df(n_jobs=1)

    assert df[("General", "FrameID", "")].tolist() == [1, 0]


def test_to_summary_df_normalizes_negative_frame_sequence() -> None:
    batch = cast(
        Any,
        FileBatchModelDisk(
            cast(
                Any,
                [
                    FakeSummaryDiskFile("/tmp/a.log", 3),
                    FakeSummaryDiskFile("/tmp/b.log", 2),
                ],
            )
        ),
    )

    df = batch.to_summary_df(frame=[0, -1], n_jobs=1)

    assert df[("General", "FrameID", "")].tolist() == [0, 2, 0, 1]


def test_to_summary_df_brief_and_flatten_columns() -> None:
    batch = cast(
        Any,
        FileBatchModelDisk(cast(Any, [FakeSummaryDiskFile("/tmp/a.log", 1)])),
    )

    df = batch.to_summary_df(frame=0, n_jobs=1, brief=False, flatten_columns=True)

    assert df.columns.tolist() == ["General.FrameID", "General.Brief"]
    assert df["General.Brief"].tolist() == [False]


def test_to_summary_df_mode_file_passes_brief() -> None:
    batch = cast(
        Any,
        FileBatchModelDisk(cast(Any, [FakeSummaryDiskFile("/tmp/a.log", 1)])),
    )

    df = batch.to_summary_df(mode="file", n_jobs=1, brief=False, flatten_columns=True)

    assert df["File.Brief"].tolist() == [False]


def test_to_summary_df_missing_frame_error() -> None:
    batch = cast(
        Any,
        FileBatchModelDisk(cast(Any, [FakeSummaryDiskFile("/tmp/a.log", 1)])),
    )

    with pytest.raises(IndexError, match="Frame index 3 is out of range"):
        batch.to_summary_df(frame=3, n_jobs=1, on_missing_frame="error")


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


def test_single_file_parser_falls_back_only_on_format_mismatch() -> None:
    parsed = filebatchparserdisk_module.single_file_parser(
        file_path="/tmp/sample.out",
        possible_readers=(
            FailingReader("orcaout", FormatMismatchError("not ORCA")),
            FakeDiskFileReader("g16log"),
        ),
    )

    assert parsed is not None
    assert parsed.detected_format_id == "g16log"


def test_single_file_parser_does_not_fall_back_on_parse_error() -> None:
    parsed = filebatchparserdisk_module.single_file_parser(
        file_path="/tmp/sample.out",
        possible_readers=(
            FailingReader("orcaout", RuntimeError("parser bug")),
            FakeDiskFileReader("g16log"),
        ),
    )

    assert parsed is None


def test_single_file_parser_outcome_preserves_reader_warnings() -> None:
    warning = ParseWarning(code="TEST.WARNING", message="reader warning")

    class WarningReader(FakeDiskFileReader):
        def read(self, path: str | Path, **_kwargs: Any) -> ParseResult[object]:
            return ParseResult(
                value=FakeDiskFile(str(path), "xyz", self.format_id),
                level=StructureLevel.COORDS,
                warnings=(warning,),
                detected_format=self.format_id,
            )

    outcome = filebatchparserdisk_module.single_file_parser(
        file_path="/tmp/sample.xyz",
        possible_readers=(WarningReader("xyz"),),
        return_outcome=True,
    )

    assert isinstance(outcome, FileParseOutcome)
    assert outcome.status == "ok"
    assert outcome.warnings == (warning,)


def test_autoparser_report_distinguishes_success_and_missing_input(tmp_path: Path) -> None:
    valid_path = tmp_path / "valid.xyz"
    valid_path.write_text("1\nwater\nH 0.0 0.0 0.0\n", encoding="utf-8")
    missing_path = tmp_path / "missing.xyz"

    report = io_module.AutoParser(
        [valid_path, missing_path],
        n_jobs=1,
        parser_detection="xyz",
        return_report=True,
    )

    assert isinstance(report, BatchParseResult)
    assert report.batch.file_paths == [str(valid_path.resolve())]
    assert [outcome.status for outcome in report.outcomes] == ["missing", "ok"]
    assert report.failures[0].failure is not None
    assert report.failures[0].failure.kind == "missing_file"


def test_auto_detection_falls_back_from_orca_to_gaussian_for_shared_out_suffix() -> None:
    batch = FileBatchParserDisk(n_jobs=1).parse(
        [Path("tests/test_files/g16irc/irc.out")],
        parser_detection="auto",
        only_last_frame=True,
    )

    assert len(batch) == 1
    assert batch[0].detected_format_id == "g16log"
    assert batch[0].qm_software == "Gaussian"


def test_simple_coordinate_readers_defer_mismatch_to_parse_phase() -> None:
    from molop.io.logic.coords.parsers.SDFFileParser import SDFFileParserDisk
    from molop.io.logic.coords.parsers.SMIFileParser import SMIFileParserDisk
    from molop.io.logic.coords.parsers.XYZFileParser import XYZFileParserDisk

    for parser_cls in (XYZFileParserDisk, SMIFileParserDisk, SDFFileParserDisk):
        parser_cls._quick_check_file_format("/tmp/not-treated-as-a-path")


def test_filter_state_no_img_requires_frequency_frames() -> None:
    batch = cast(
        Any,
        FileBatchModelDisk(
            cast(
                Any,
                [
                    FakeStateDiskFile("/tmp/no-freq.log", [FakeFrame(None)]),
                    FakeStateDiskFile("/tmp/no-img.log", [FakeFrame(None), FakeFrame(0)]),
                    FakeStateDiskFile("/tmp/empty-freq.log", [FakeFrame(0, num_frequencies=0)]),
                    FakeStateDiskFile("/tmp/imaginary.log", [FakeFrame(1)]),
                ],
            )
        ),
    )

    filtered = batch.filter_state("no-img")

    assert filtered.file_paths == ["/tmp/no-img.log"]


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


def test_filebatchparser_progress_updates_after_parse_completion(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    parser = FileBatchParserDisk(n_jobs=1)
    file_paths = ["/tmp/a.xyz", "/tmp/b.xyz"]
    events: list[str] = []
    original_new_batch_from_sorted = (
        filebatchparserdisk_module.FileBatchModelDisk._new_batch_from_sorted_diskfiles
    )

    class FakeProgress:
        total: int | None

        def __init__(self, desc: str, total: int | None = None) -> None:
            self.desc = desc
            self.total = total
            events.append(f"start:{desc}")

        def update(self, value: int) -> None:
            events.append(f"progress:{self.desc}:{value}")

        def refresh(self) -> None:
            return None

        def close(self) -> None:
            events.append(f"closed:{self.desc}")

    def fake_progress(iterable: Any = None, **kwargs: Any) -> Any:
        if iterable is not None:
            raise AssertionError("FileBatchParserDisk should track completion manually")
        return FakeProgress(desc=kwargs["desc"], total=kwargs.get("total"))

    def fake_parse(**task: Any) -> FakeDiskFile:
        events.append(f"parse:{task['file_path']}")
        return FakeDiskFile(task["file_path"], "xyz", "xyz")

    def fake_new_batch_from_sorted(cls: Any, diskfiles: Any) -> FileBatchModelDisk[Any]:
        events.append("new_batch_from_sorted:start")
        batch = original_new_batch_from_sorted(list(diskfiles))
        events.append("new_batch_from_sorted:end")
        return batch

    monkeypatch.setattr(filebatchparserdisk_module, "AdaptiveProgress", fake_progress)
    monkeypatch.setattr(
        filebatchparserdisk_module.FileBatchModelDisk,
        "_new_batch_from_sorted_diskfiles",
        classmethod(fake_new_batch_from_sorted),
    )
    monkeypatch.setattr(filebatchparserdisk_module.os.path, "isfile", lambda path: True)
    monkeypatch.setattr(filebatchparserdisk_module.os.path, "abspath", lambda path: cast(str, path))
    monkeypatch.setattr(
        filebatchparserdisk_module.codec_registry,
        "select_reader",
        lambda _path, hint_format=None: (DummyReader("xyz", frozenset({".xyz"}), 1),),
    )
    monkeypatch.setattr(filebatchparserdisk_module, "single_file_parser", fake_parse)

    batch = parser.parse(file_paths)

    assert batch.file_paths == file_paths
    assert events == [
        "start:MolOP parsing with single process",
        "parse:/tmp/a.xyz",
        "progress:MolOP parsing with single process:1",
        "parse:/tmp/b.xyz",
        "progress:MolOP parsing with single process:1",
        "closed:MolOP parsing with single process",
        "new_batch_from_sorted:start",
        "new_batch_from_sorted:end",
    ]


def test_autoparser_glob_materializes_paths_for_progress_total(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    (tmp_path / "a.log").write_text("", encoding="utf-8")
    (tmp_path / "b.log").write_text("", encoding="utf-8")
    captured: dict[str, Any] = {}

    class StubFileBatchParserDisk:
        def __init__(self, n_jobs: int) -> None:
            captured["n_jobs"] = n_jobs

        def parse(self, file_paths: Any, **_kwargs: Any) -> FileBatchModelDisk[Any]:
            captured["is_list"] = isinstance(file_paths, list)
            captured["path_count"] = len(file_paths)
            return FileBatchModelDisk()

    monkeypatch.setattr(io_module, "FileBatchParserDisk", StubFileBatchParserDisk)

    io_module.AutoParser(str(tmp_path / "*.log"), n_jobs=3)

    assert captured == {
        "n_jobs": 3,
        "is_list": True,
        "path_count": 2,
    }


def _capture_autoparser_parse_call(monkeypatch: pytest.MonkeyPatch) -> dict[str, Any]:
    captured: dict[str, Any] = {}

    class StubFileBatchParserDisk:
        def __init__(self, n_jobs: int) -> None:
            captured["n_jobs"] = n_jobs

        def parse(self, file_paths: Any, **kwargs: Any) -> FileBatchModelDisk[Any]:
            captured["file_paths"] = file_paths
            captured["parse_kwargs"] = kwargs
            return FileBatchModelDisk()

    monkeypatch.setattr(io_module, "FileBatchParserDisk", StubFileBatchParserDisk)
    return captured


def test_autoparser_expands_mixed_path_and_glob_iterable_to_sorted_unique_absolute_paths(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    explicit_string = tmp_path / "c.xyz"
    explicit_path = tmp_path / "b.log"
    glob_only = tmp_path / "a.log"
    for path in (explicit_string, explicit_path, glob_only):
        path.write_text("", encoding="utf-8")
    captured = _capture_autoparser_parse_call(monkeypatch)

    result = io_module.AutoParser(
        [str(explicit_string), explicit_path, str(tmp_path / "*.log")],
        n_jobs=2,
    )

    file_paths = captured["file_paths"]
    assert isinstance(file_paths, list)
    assert [Path(path) for path in file_paths] == sorted(
        [glob_only.resolve(), explicit_path.resolve(), explicit_string.resolve()]
    )
    assert all(Path(path).is_absolute() for path in file_paths)
    assert captured["n_jobs"] == 2
    assert isinstance(result, FileBatchModelDisk)


def test_autoparser_materializes_generator_inputs_before_parsing(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    explicit = tmp_path / "b.xyz"
    globbed = tmp_path / "a.log"
    explicit.write_text("", encoding="utf-8")
    globbed.write_text("", encoding="utf-8")
    yielded: list[str | Path] = []
    inputs = [explicit, str(tmp_path / "*.log")]

    def generate_inputs():
        for value in inputs:
            yielded.append(value)
            yield value

    captured = _capture_autoparser_parse_call(monkeypatch)

    io_module.AutoParser(generate_inputs(), n_jobs=1)

    assert yielded == inputs
    assert isinstance(captured["file_paths"], list)
    assert [Path(path) for path in captured["file_paths"]] == sorted(
        [globbed.resolve(), explicit.resolve()]
    )


def test_autoparser_deduplicates_equivalent_relative_and_absolute_paths(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    path = tmp_path / "sample.log"
    path.write_text("", encoding="utf-8")
    monkeypatch.chdir(tmp_path)
    captured = _capture_autoparser_parse_call(monkeypatch)

    io_module.AutoParser(["sample.log", path.resolve(), Path("sample.log")], n_jobs=1)

    file_paths = captured["file_paths"]
    assert isinstance(file_paths, list)
    assert [Path(value) for value in file_paths] == [path.resolve()]


def test_autoparser_empty_and_unmatched_inputs_produce_empty_path_list(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured = _capture_autoparser_parse_call(monkeypatch)

    for value in ([], [str(tmp_path / "*.missing")], str(tmp_path / "*.missing")):
        io_module.AutoParser(value, n_jobs=1)
        assert captured["file_paths"] == []


def test_autoparser_passes_missing_literal_paths_to_batch_parser(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    missing = tmp_path / "missing.log"
    captured = _capture_autoparser_parse_call(monkeypatch)

    for value in (missing, [missing]):
        io_module.AutoParser(value, n_jobs=1)
        assert [Path(path) for path in captured["file_paths"]] == [missing.resolve()]


def test_autoparser_invalid_iterable_member_reports_its_index(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    valid = tmp_path / "valid.log"
    valid.write_text("", encoding="utf-8")
    captured = _capture_autoparser_parse_call(monkeypatch)

    with pytest.raises(TypeError) as exc_info:
        io_module.AutoParser([valid, object()], n_jobs=1)

    assert "[1]" in str(exc_info.value)
    assert "object" in str(exc_info.value)
    assert "file_paths" not in captured


def test_filebatchparser_parallel_orders_dispatch_buffer_by_file_size(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    parser = FileBatchParserDisk(n_jobs=2)
    file_paths = [
        "/tmp/tiny.xyz",
        "/tmp/large.xyz",
        "/tmp/small.xyz",
        "/tmp/medium.xyz",
        "/tmp/larger.xyz",
    ]
    file_sizes = {
        "/tmp/tiny.xyz": 1,
        "/tmp/large.xyz": 100,
        "/tmp/small.xyz": 10,
        "/tmp/medium.xyz": 50,
        "/tmp/larger.xyz": 80,
    }
    parsed_order: list[str] = []

    monkeypatch.setattr(filebatchparserdisk_module.os.path, "isfile", lambda path: True)
    monkeypatch.setattr(filebatchparserdisk_module.os.path, "abspath", lambda path: cast(str, path))
    monkeypatch.setattr(
        filebatchparserdisk_module.os.path,
        "getsize",
        lambda path: file_sizes[cast(str, path)],
    )
    monkeypatch.setattr(
        filebatchparserdisk_module.codec_registry,
        "select_reader",
        lambda _path, hint_format=None: (DummyReader("xyz", frozenset({".xyz"}), 1),),
    )

    def fake_parse(**task: Any) -> FakeDiskFile:
        parsed_order.append(task["file_path"])
        return FakeDiskFile(task["file_path"], "xyz", "xyz")

    monkeypatch.setattr(filebatchparserdisk_module, "single_file_parser", fake_parse)
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

    assert parsed_order == [
        "/tmp/large.xyz",
        "/tmp/larger.xyz",
        "/tmp/medium.xyz",
        "/tmp/small.xyz",
        "/tmp/tiny.xyz",
    ]
    assert batch.file_paths == parsed_order


def test_filebatchparser_tunes_parallel_jobs_when_task_count_matches_requested_jobs(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    parser = FileBatchParserDisk(n_jobs=3)
    file_paths = ["/tmp/a.xyz", "/tmp/b.xyz", "/tmp/c.xyz"]
    captured: dict[str, Any] = {}

    monkeypatch.setattr(filebatchparserdisk_module.os.path, "isfile", lambda path: True)
    monkeypatch.setattr(filebatchparserdisk_module.os.path, "abspath", lambda path: cast(str, path))
    monkeypatch.setattr(
        filebatchparserdisk_module.codec_registry,
        "select_reader",
        lambda _path, hint_format=None: (DummyReader("xyz", frozenset({".xyz"}), 1),),
    )
    monkeypatch.setattr(
        filebatchparserdisk_module,
        "single_file_parser",
        lambda **task: FakeDiskFile(task["file_path"], "xyz", "xyz"),
    )
    monkeypatch.setattr(
        filebatchparserdisk_module, "delayed", lambda func: lambda **task: lambda: func(**task)
    )

    class StubParallel:
        def __init__(self, **kwargs: Any) -> None:
            captured.update(kwargs)

        def __call__(self, iterable: Any) -> Any:
            return (callable_obj() for callable_obj in iterable)

    monkeypatch.setattr(filebatchparserdisk_module, "Parallel", StubParallel)

    batch = parser.parse(file_paths)

    assert captured["n_jobs"] == 2
    assert batch.file_paths == ["/tmp/a.xyz", "/tmp/b.xyz", "/tmp/c.xyz"]


def test_filebatchparser_does_not_stat_file_sizes_during_scheduling(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    parser = FileBatchParserDisk(n_jobs=1)
    file_paths = ["/tmp/good.xyz", "/tmp/no-size-scan.xyz"]

    monkeypatch.setattr(filebatchparserdisk_module.os.path, "isfile", lambda path: True)
    monkeypatch.setattr(filebatchparserdisk_module.os.path, "abspath", lambda path: cast(str, path))
    monkeypatch.setattr(
        filebatchparserdisk_module.os.path,
        "getsize",
        lambda _path: pytest.fail("FileBatchParserDisk should not stat file sizes"),
    )
    monkeypatch.setattr(
        filebatchparserdisk_module.codec_registry,
        "select_reader",
        lambda _path, hint_format=None: (DummyReader("xyz", frozenset({".xyz"}), 1),),
    )

    monkeypatch.setattr(
        filebatchparserdisk_module,
        "single_file_parser",
        lambda **task: FakeDiskFile(task["file_path"], "xyz", "xyz"),
    )

    batch = parser.parse(file_paths)

    assert batch.file_paths == ["/tmp/good.xyz", "/tmp/no-size-scan.xyz"]


def test_filebatchparser_auto_detection_filters_candidates_with_reader_probe(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    parser = FileBatchParserDisk(n_jobs=1)

    monkeypatch.setattr(filebatchparserdisk_module.os.path, "isfile", lambda path: True)
    monkeypatch.setattr(filebatchparserdisk_module.os.path, "abspath", lambda path: cast(str, path))
    monkeypatch.setattr(
        filebatchparserdisk_module.codec_registry,
        "select_reader",
        lambda _path, hint_format=None: (
            ProbingReader("orcaout", False),
            ProbingReader("g16log", True),
        ),
    )

    def fake_parse(**task: Any) -> FakeDiskFile:
        reader_ids = [reader.format_id for reader in task["possible_readers"]]
        assert reader_ids == ["g16log"]
        return FakeDiskFile(task["file_path"], "log", "g16log")

    monkeypatch.setattr(filebatchparserdisk_module, "single_file_parser", fake_parse)

    batch = parser.parse(["/tmp/shared.log"], parser_detection="auto")

    assert batch.file_paths == ["/tmp/shared.log"]
    assert batch[0].detected_format_id == "g16log"


def test_base_file_parser_disk_probe_uses_parser_quick_check() -> None:
    from molop.io.logic.gaussian.log.parsers.G16LogFileParser import G16LogFileParserDisk

    assert G16LogFileParserDisk.probe_file_format(
        "tests/test_files/g16log/000000000000_000016928457_00_conf_01_ts.107c60f3cfcb.log"
    )
    assert not G16LogFileParserDisk.probe_file_format(
        "tests/test_files/xyz/dsgdb9nsd_125600-5/0.xyz"
    )


def test_base_file_parser_finalizes_file_charge_and_multiplicity_from_first_frame() -> None:
    from molop.io.logic.coords.parsers.XYZFileParser import XYZFileParserMemory

    parsed = XYZFileParserMemory().parse(
        "\n".join(
            [
                "1",
                "charge 1 multiplicity 2",
                "H 0.0 0.0 0.0",
            ]
        )
    )

    assert parsed.charge == 1
    assert parsed.multiplicity == 2


def test_g16_file_parser_finalizes_temperature_from_later_frame() -> None:
    from molop.io.logic.gaussian.log.frame_models.G16LogFileFrame import G16LogFileFrameMemory
    from molop.io.logic.gaussian.log.parsers.G16LogFileParser import G16LogFileParserMemory

    parser = G16LogFileParserMemory()
    frame_data = {
        "atoms": [1],
        "coords": [[0.0, 0.0, 0.0]] * atom_ureg.angstrom,
        "charge": 0,
        "multiplicity": 1,
    }
    first = G16LogFileFrameMemory.model_validate(frame_data)
    second = G16LogFileFrameMemory.model_validate(frame_data | {"temperature": 350 * atom_ureg.K})
    parsed = parser._chem_file.model_validate({"qm_software": "Gaussian"})
    parsed.append(first)
    parsed.append(second)

    parser._update_file_metadata_from_frames(parsed, {})

    assert parsed.temperature is not None
    assert parsed.temperature.to("K").m == 350


def test_g16_file_parser_preserves_segment_status_over_frame_status() -> None:
    from molop.io.base_models.DataClasses import Status
    from molop.io.logic.gaussian.log.frame_models.G16LogFileFrame import G16LogFileFrameMemory
    from molop.io.logic.gaussian.log.parsers.G16LogFileParser import G16LogFileParserMemory

    parser = G16LogFileParserMemory()
    frame_data = {
        "atoms": [1],
        "coords": [[0.0, 0.0, 0.0]] * atom_ureg.angstrom,
        "charge": 0,
        "multiplicity": 1,
    }
    first = G16LogFileFrameMemory.model_validate(
        frame_data | {"status": Status(normal_terminated=True)}
    )
    second = G16LogFileFrameMemory.model_validate(
        frame_data | {"status": Status(normal_terminated=False)}
    )
    parsed = parser._chem_file.model_validate(
        {"qm_software": "Gaussian", "status": Status(normal_terminated=True)}
    )
    parsed.append(first)
    parsed.append(second)

    parser._update_file_metadata_from_frames(parsed, {"status": Status(normal_terminated=True)})

    assert parsed.status.normal_terminated is True


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
    rendered = frame.format_transform("xyz", file_path=target, write_to_disk=True)

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
