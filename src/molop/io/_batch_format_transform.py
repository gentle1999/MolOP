"""
Author: TMJ
Date: 2026-02-15 20:33:01
LastEditors: TMJ
LastEditTime: 2026-02-15 22:14:01
Description: 请填写简介
"""

from __future__ import annotations

import os
from collections.abc import Sequence
from typing import Any, Literal, Protocol, cast

from molop.config import moloplogger


class _HasParallelExecute(Protocol):
    def parallel_execute(self, func: Any, desc: str = "", n_jobs: int = 1) -> list[Any]: ...


class BatchFormatTransformMixin:
    def format_transform(
        self,
        format: str,
        output_dir: str | None = None,
        frameID: int | Literal["all"] | Sequence[int] = -1,
        embed_in_one_file: bool = True,
        n_jobs: int = 1,
        **kwargs: Any,
    ) -> dict[str, str | list[str]]:

        if output_dir is not None:
            assert os.path.isdir(output_dir), f"{output_dir} is not a directory"

        typed_self = cast(_HasParallelExecute, self)

        def transform_func(diskfile: Any) -> tuple[str, str | list[str]]:
            frame_ids: object = range(len(diskfile)) if frameID == "all" else frameID
            try:
                res = diskfile.format_transform(
                    format,
                    frameID=frame_ids,
                    embed_in_one_file=embed_in_one_file,
                    file_path=os.path.join(
                        output_dir,
                        f"{os.path.splitext(diskfile.filename)[0]}.{format}",
                    )
                    if output_dir
                    else None,
                    **kwargs,
                )
                return diskfile.file_path, res
            except Exception as e:
                moloplogger.warning(
                    f"Format transform failed for {diskfile.filename}: {type(e).__name__}: {e}"
                )
                return diskfile.file_path, ("" if embed_in_one_file else [])

        desc = f"MolOP processing {format} format with {n_jobs} jobs"
        results = typed_self.parallel_execute(transform_func, desc, n_jobs)
        return dict(results)
