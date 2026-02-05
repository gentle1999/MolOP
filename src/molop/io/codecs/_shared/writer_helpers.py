from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any, Literal, cast

from molop.io.base_models.ChemFile import BaseChemFile, BaseCoordsFile
from molop.io.base_models.ChemFileFrame import BaseChemFileFrame, BaseCoordsFrame
from molop.io.base_models.Mixins import DiskStorageMixin, FileMixin, _HasRenderableFrames
from molop.io.codec_types import StructureLevel, WriterCodec


def clone_file_and_frames(
    source: object,
    file_cls: type[BaseChemFile],
    frame_cls: type[BaseChemFileFrame],
) -> Any:
    source_file = cast(BaseChemFile, source)
    cloned_file = file_cls.model_validate(source_file.model_dump())
    for frame in source_file.frames:
        cloned_frame = frame_cls.model_validate(
            cast(BaseChemFileFrame, frame).model_dump(exclude={"file_path"})
        )
        cloned_file.append(cast(Any, cloned_frame))
    return cloned_file


@dataclass
class FileRendererWriter:
    format_id: str
    required_level: StructureLevel
    file_cls: type[Any]
    frame_cls: type[BaseChemFileFrame]
    priority: int

    def write(
        self,
        value: object,
        *,
        frameID: Sequence[int] | int | Literal["all"] = -1,
        embed_in_one_file: bool = True,
        **kwargs: Any,
    ) -> object:
        if not hasattr(value, "model_dump") or not hasattr(value, "frames"):
            raise TypeError("Writer requires a BaseChemFile-compatible input.")
        cloned_file = clone_file_and_frames(value, self.file_cls, self.frame_cls)
        return FileMixin._render(
            cast(Any, cloned_file),
            frameID=frameID,
            embed_in_one_file=embed_in_one_file,
            **kwargs,
        )


class CMLFileFrameMixin:
    def _render(self, engine: Literal["rdkit", "openbabel"] = "rdkit", **kwargs: Any) -> str:
        typed_self = cast(BaseCoordsFrame, self)
        return typed_self.to_CML_block(engine=engine)


class CMLFileFrameDisk(
    DiskStorageMixin, CMLFileFrameMixin, BaseCoordsFrame["CMLFileFrameDisk"]
): ...


class CMLFileMixin(FileMixin):
    def _render_frames_in_one_file(self, frameID: Sequence[int], **kwargs: Any) -> str:
        typed_self = cast(_HasRenderableFrames, self)
        return "\n".join(
            frame._render(**kwargs) for frame in typed_self.frames if frame.frame_id in frameID
        )

    def _render_frames(self, frameID: Sequence[int], **kwargs: Any) -> list[str]:
        typed_self = cast(_HasRenderableFrames, self)
        return [frame._render(**kwargs) for frame in typed_self.frames if frame.frame_id in frameID]


class CMLFileDisk(DiskStorageMixin, CMLFileMixin, BaseCoordsFile[CMLFileFrameDisk]): ...


__all__ = [
    "CMLFileDisk",
    "CMLFileFrameDisk",
    "CMLFileFrameMixin",
    "CMLFileMixin",
    "FileRendererWriter",
    "StructureLevel",
    "WriterCodec",
    "clone_file_and_frames",
]
