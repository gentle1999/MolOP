from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any, Literal, cast

from rdkit import Chem

from molop.io.codec_registry import Registry
from molop.io.codec_types import StructureLevel, WriterCodec


def _normalize_frame_ids(value: object, frameID: Sequence[int] | int | Literal["all"]) -> list[int]:
    typed_value = cast(Any, value)
    frame_count = len(typed_value.frames)
    if isinstance(frameID, int):
        return [frameID if frameID >= 0 else frame_count + frameID]
    if frameID == "all":
        return list(range(frame_count))
    return list(frameID)


def _render_cml_frame(frame: object, *, engine: Literal["rdkit", "openbabel"] = "rdkit") -> str:
    typed_frame = cast(Any, frame)
    if engine == "rdkit":
        rdmol = getattr(typed_frame, "qm_embedded_rdmol", getattr(typed_frame, "rdmol", None))
        if rdmol is None:
            raise ValueError("CML building failed. No RDKit molecule recovered.")
        return Chem.MolToMrvBlock(rdmol)
    if engine == "openbabel":
        omol = getattr(typed_frame, "omol", None)
        if omol is None:
            raise ValueError("CML building failed. No Openbabel molecule recovered.")
        cml_text = omol.write("cml")
        if cml_text is None:
            raise ValueError("CML building failed. No CML text recovered.")
        return cml_text
    raise ValueError(f"Unsupported engine: {engine}")


@dataclass(frozen=True, slots=True)
class CMLWriter:
    format_id: str
    required_level: StructureLevel
    priority: int

    def write(
        self,
        value: object,
        *,
        frameID: Sequence[int] | int | Literal["all"] = -1,
        embed_in_one_file: bool = True,
        engine: Literal["rdkit", "openbabel"] = "rdkit",
        **kwargs: Any,
    ) -> object:
        if kwargs.get("file_path") is not None:
            kwargs.pop("file_path")
        if not hasattr(value, "frames"):
            raise TypeError("CML writer requires a BaseChemFile-compatible input.")
        typed_value = cast(Any, value)
        selected_frame_ids = _normalize_frame_ids(value, frameID)
        rendered_frames = [
            _render_cml_frame(frame, engine=engine)
            for frame in typed_value.frames
            if getattr(frame, "frame_id", None) in selected_frame_ids
        ]
        if embed_in_one_file:
            return "\n".join(rendered_frames)
        return rendered_frames


@dataclass(frozen=True, slots=True)
class CMLFrameWriter:
    format_id: str
    required_level: StructureLevel
    priority: int

    def write(
        self,
        value: object,
        *,
        engine: Literal["rdkit", "openbabel"] = "rdkit",
        **kwargs: Any,
    ) -> object:
        _ = kwargs.pop("file_path", None)
        kwargs.pop("frameID", None)
        kwargs.pop("embed_in_one_file", None)
        return _render_cml_frame(value, engine=engine)


def register(registry: Registry) -> None:
    priority = 100

    @registry.writer_factory(
        format_id="cml",
        required_level=StructureLevel.GRAPH,
        domain="file",
        priority=priority,
    )
    def _factory() -> WriterCodec:
        return cast(
            WriterCodec,
            CMLWriter(
                format_id="cml",
                required_level=StructureLevel.GRAPH,
                priority=priority,
            ),
        )

    @registry.writer_factory(
        format_id="cml",
        required_level=StructureLevel.GRAPH,
        domain="frame",
        priority=priority,
    )
    def _frame_factory() -> WriterCodec:
        return cast(
            WriterCodec,
            CMLFrameWriter(
                format_id="cml",
                required_level=StructureLevel.GRAPH,
                priority=priority,
            ),
        )


__all__ = ["register"]
