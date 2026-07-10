from __future__ import annotations

import os
from collections.abc import Sequence
from typing import Any, Protocol, cast

from molop.io import codec_registry
from molop.io.codec_types import GraphPolicy
from molop.io.frame_selection import FrameSelector, normalize_frame_selector


class _HasFrames(Protocol):
    @property
    def frames(self) -> Sequence[Any]: ...


class _HasFramePayload(Protocol):
    def model_dump(self, **kwargs: Any) -> dict[str, Any]: ...


def _source_file_path(value: object) -> str | None:
    file_path = getattr(value, "file_path", None)
    if isinstance(file_path, str) and file_path:
        return file_path
    return None


def _resolve_format_output_path(
    value: object,
    format: str,
    file_path: os.PathLike | str | None = None,
) -> str:
    source_path = os.fspath(file_path) if file_path is not None else _source_file_path(value)
    if source_path is None:
        raise ValueError("file_path is required when writing a memory-only object to disk")
    assert not os.path.isdir(source_path), "file_path should be a file path or None"
    dir_path = os.path.dirname(source_path)
    base = os.path.basename(source_path).split(".")[0]
    return os.path.join(dir_path, f"{base}.{format}")


class FrameFormatTransformMixin:
    def format_transform(
        self,
        format: str,
        file_path: os.PathLike | str | None = None,
        write_to_disk: bool = False,
        **kwargs: Any,
    ) -> str:
        normalized_format = format.strip().lower()
        graph_policy = kwargs.pop("graph_policy", None)
        writer_file_path = (
            _resolve_format_output_path(self, normalized_format, file_path)
            if write_to_disk
            else None
        )
        rendered = self._render_frame_format(
            normalized_format,
            file_path=writer_file_path,
            graph_policy=graph_policy,
            **kwargs,
        )
        if write_to_disk:
            output_path = cast(str, writer_file_path)
            with open(output_path, "w") as f:
                f.write(rendered)
        return rendered

    def _render_frame_format(
        self,
        format: str,
        *,
        file_path: os.PathLike | str | None = None,
        graph_policy: GraphPolicy | None = None,
        **kwargs: Any,
    ) -> str:
        rendered = codec_registry.write_frame(
            format,
            cast(_HasFramePayload, self),
            graph_policy=graph_policy,
            file_path=os.fspath(file_path) if file_path is not None else None,
            **kwargs,
        )
        if not isinstance(rendered, str):
            raise TypeError(
                f"Frame format transform expected str output for {format}, got {type(rendered)}"
            )
        return rendered


class FormatTransformMixin:
    def format_transform(
        self,
        format: str,
        frame: FrameSelector = -1,
        file_path: os.PathLike | str | None = None,
        embed_in_one_file: bool = True,
        write_to_disk: bool = False,
        **kwargs,
    ) -> str | list[str]:
        typed_self = cast(_HasFrames, self)
        frame_ids = normalize_frame_selector(
            frame,
            len(typed_self.frames),
            parameter_name="frame",
        )

        graph_policy = kwargs.pop("graph_policy", None)
        write_kwargs = dict(kwargs)
        writer_file_path = (
            _resolve_format_output_path(self, format, file_path) if write_to_disk else None
        )
        if writer_file_path is not None:
            write_kwargs["file_path"] = writer_file_path
        rendered = cast(
            str | list[str],
            codec_registry.write(
                format,
                self,
                frame=frame_ids,
                embed_in_one_file=embed_in_one_file,
                graph_policy=graph_policy,
                **write_kwargs,
            ),
        )
        if write_to_disk:
            output_file_path = cast(str, writer_file_path)
            dir_path = os.path.dirname(output_file_path)
            base = os.path.basename(output_file_path).split(".")[0]
            if isinstance(rendered, str):
                output_path = os.path.join(dir_path, f"{base}.{format}")
                with open(output_path, "w") as f:
                    f.write(rendered)
            elif isinstance(rendered, list):
                for idx, frame_content in zip(frame_ids, rendered, strict=True):
                    filename = base + f"{idx:03d}.{format}"
                    output_path = os.path.join(dir_path, filename)
                    with open(output_path, "w") as f:
                        f.write(frame_content)
        return rendered
