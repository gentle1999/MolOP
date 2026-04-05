from __future__ import annotations

from typing import Any

import pytest

from molop.io.codec_types import StructureLevel
from molop.io.codecs.cml_codec import CMLFrameWriter, CMLWriter


class DummyFrame:
    def __init__(self, frame_id: int, *, rdmol: object | None = None, omol: Any = None) -> None:
        self.frame_id = frame_id
        self.rdmol = rdmol
        self.omol = omol


class DummyFile:
    def __init__(self, frames: list[DummyFrame]) -> None:
        self.frames = frames


def test_cml_writer_renders_selected_frames_in_one_file(monkeypatch: pytest.MonkeyPatch) -> None:
    writer = CMLWriter(format_id="cml", required_level=StructureLevel.GRAPH, priority=100)
    monkeypatch.setattr("molop.io.codecs.cml_codec.Chem.MolToMrvBlock", lambda mol: f"cml:{mol}")

    rendered = writer.write(
        DummyFile([DummyFrame(0, rdmol="a"), DummyFrame(1, rdmol="b")]),
        frameID="all",
        embed_in_one_file=True,
    )

    assert rendered == "cml:a\ncml:b"


def test_cml_writer_supports_openbabel_engine() -> None:
    writer = CMLWriter(format_id="cml", required_level=StructureLevel.GRAPH, priority=100)
    omol = type("DummyOMol", (), {"write": lambda self, fmt: "<cml />" if fmt == "cml" else None})()

    rendered = writer.write(
        DummyFile([DummyFrame(0, omol=omol)]),
        frameID=[0],
        embed_in_one_file=False,
        engine="openbabel",
    )

    assert rendered == ["<cml />"]


def test_cml_writer_requires_frames() -> None:
    writer = CMLWriter(format_id="cml", required_level=StructureLevel.GRAPH, priority=100)

    with pytest.raises(TypeError, match="BaseChemFile-compatible"):
        writer.write(object())


def test_cml_frame_writer_renders_single_frame(monkeypatch: pytest.MonkeyPatch) -> None:
    writer = CMLFrameWriter(format_id="cml", required_level=StructureLevel.GRAPH, priority=100)
    monkeypatch.setattr("molop.io.codecs.cml_codec.Chem.MolToMrvBlock", lambda mol: f"cml:{mol}")

    rendered = writer.write(DummyFrame(0, rdmol="frame"), engine="rdkit")

    assert rendered == "cml:frame"
