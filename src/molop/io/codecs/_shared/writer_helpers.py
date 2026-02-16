from __future__ import annotations

import os
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
    *,
    file_path: str | None = None,
    frame_file_paths: dict[int, str] | None = None,
) -> Any:
    source_file = cast(BaseChemFile, source)
    file_payload = source_file.model_dump()
    if file_path is not None and issubclass(file_cls, DiskStorageMixin):
        file_payload["file_path"] = file_path
    cloned_file = file_cls.model_validate(file_payload)
    for frame in source_file.frames:
        typed_frame = cast(BaseChemFileFrame, frame)
        frame_payload = typed_frame.model_dump(exclude={"file_path"})
        if frame_file_paths is not None and issubclass(frame_cls, DiskStorageMixin):
            maybe_frame_file_path = frame_file_paths.get(typed_frame.frame_id)
            if maybe_frame_file_path is not None:
                frame_payload["file_path"] = maybe_frame_file_path
        cloned_frame = frame_cls.model_validate(frame_payload)

        source_private = getattr(typed_frame, "__pydantic_private__", None)
        if isinstance(source_private, dict):
            for key, value in source_private.items():
                if isinstance(key, str) and key.startswith("_orca_") and hasattr(cloned_frame, key):
                    setattr(cloned_frame, key, value)
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
        output_file_path = kwargs.pop("file_path", None)
        cloned_file = clone_file_and_frames(
            value,
            self.file_cls,
            self.frame_cls,
            file_path=self._derive_file_path(output_file_path),
            frame_file_paths=self._derive_frame_file_paths(
                value, output_file_path, embed_in_one_file
            ),
        )
        return FileMixin._render(
            cast(Any, cloned_file),
            frameID=frameID,
            embed_in_one_file=embed_in_one_file,
            **kwargs,
        )

    def _derive_file_path(self, file_path: os.PathLike[str] | str | None) -> str | None:
        if file_path is None:
            return None
        normalized_file_path = os.fspath(file_path)
        dir_path = os.path.dirname(normalized_file_path)
        base = os.path.basename(normalized_file_path).split(".")[0]
        return os.path.join(dir_path, f"{base}.{self.format_id}")

    def _derive_frame_file_paths(
        self,
        value: object,
        file_path: os.PathLike[str] | str | None,
        embed_in_one_file: bool,
    ) -> dict[int, str] | None:
        derived_file_path = self._derive_file_path(file_path)
        if derived_file_path is None:
            return None
        typed_value = cast(BaseChemFile, value)
        if embed_in_one_file:
            return {
                cast(BaseChemFileFrame, frame).frame_id: derived_file_path
                for frame in typed_value.frames
            }
        dir_path = os.path.dirname(derived_file_path)
        base = os.path.basename(derived_file_path).split(".")[0]
        return {
            frame.frame_id: os.path.join(
                dir_path,
                f"{base}{frame.frame_id:03d}.{self.format_id}",
            )
            for frame in cast(list[BaseChemFileFrame], typed_value.frames)
        }


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
