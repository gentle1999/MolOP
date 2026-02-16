from __future__ import annotations

import os
from collections.abc import Sequence
from typing import Any, Literal, Protocol, cast

from molop.io import codec_registry


class _HasFrames(Protocol):
    @property
    def frames(self) -> Sequence[Any]: ...


class FormatTransformMixin:
    def format_transform(
        self,
        format: str,
        frameID: Sequence[int] | int | Literal["all"] | slice = -1,
        file_path: os.PathLike | str | None = None,
        embed_in_one_file: bool = True,
        **kwargs,
    ) -> str | list[str]:
        typed_self = cast(_HasFrames, self)
        frame_key: object = frameID
        assert file_path is None or not os.path.isdir(file_path), (
            "file_path should be a file path or None"
        )
        if isinstance(frame_key, int):
            frame_ids = [frame_key if frame_key >= 0 else len(typed_self.frames) + frame_key]
        elif frame_key == "all":
            frame_ids = list(range(len(typed_self.frames)))
        elif isinstance(frame_key, slice):
            frame_ids = list(range(len(typed_self.frames)))[frame_key]
        elif isinstance(frame_key, Sequence):
            frame_ids = []
            for i in frame_key:
                assert isinstance(i, int), "frameID should be a sequence of integers"
                frame_ids.append(i)
        else:
            raise ValueError("frameID should be an integer, a sequence of integers, or 'all'")

        graph_policy = kwargs.pop("graph_policy", "prefer")
        write_kwargs = dict(kwargs)
        if file_path is not None:
            write_kwargs["file_path"] = os.fspath(file_path)
        rendered = cast(
            str | list[str],
            codec_registry.write(
                format,
                self,
                frameID=frame_ids,
                embed_in_one_file=embed_in_one_file,
                graph_policy=graph_policy,
                **write_kwargs,
            ),
        )
        if file_path:
            dir_path = os.path.dirname(os.fspath(file_path))
            base = os.path.basename(file_path).split(".")[0]
            if isinstance(rendered, str):
                filename = base + f".{format}"
                output_path = os.path.join(dir_path, filename)
                with open(output_path, "w") as f:
                    f.write(rendered)
            elif isinstance(rendered, list):
                for idx, frame_content in zip(frame_ids, rendered, strict=True):
                    filename = base + f"{idx:03d}.{format}"
                    output_path = os.path.join(dir_path, filename)
                    with open(output_path, "w") as f:
                        f.write(frame_content)
        return rendered
