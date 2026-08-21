"""
Author: TMJ
Date: 2026-02-15 20:33:01
LastEditors: TMJ
LastEditTime: 2026-02-15 22:14:01
Description: 请填写简介
"""

from __future__ import annotations

import os
from typing import Any, Protocol, cast

from molop.config import molopconfig, moloplogger
from molop.io import codec_registry
from molop.io.frame_selection import FrameSelector
from molop.utils.progressbar import NativeReconstructionConcurrencyError


class _HasParallelExecute(Protocol):
    def parallel_execute(
        self,
        func: Any,
        desc: str = "",
        n_jobs: int = -1,
        *,
        return_results: bool | None = None,
    ) -> Any: ...


class BatchFormatTransformMixin:
    def format_transform(
        self,
        format: str,
        output_dir: str | None = None,
        frame: FrameSelector = -1,
        embed_in_one_file: bool = True,
        write_to_disk: bool = False,
        n_jobs: int = -1,
        **kwargs: Any,
    ) -> dict[str, str | list[str]]:
        if write_to_disk and output_dir is not None:
            assert os.path.isdir(output_dir), f"{output_dir} is not a directory"

        typed_self = cast(_HasParallelExecute, self)
        output_extension = codec_registry.get_writer_output_extension(format)
        graph_policy = kwargs.get("graph_policy")
        needs_graph = codec_registry.writer_requires_graph(
            format,
            graph_policy=graph_policy,
        ) or (format.strip().lower() == "gjf" and kwargs.get("add_gjf_connectivity", False))
        effective_jobs = (
            molopconfig.set_molgr_n_jobs(n_jobs) if needs_graph else molopconfig.set_n_jobs(n_jobs)
        )
        if needs_graph and molopconfig.prewarm_topologies:
            prewarm_topologies = getattr(typed_self, "_prewarm_topologies", None)
            if callable(prewarm_topologies):
                prewarm_topologies(
                    frame=frame,
                    max_workers=effective_jobs,
                )

        def transform_func(diskfile: Any) -> tuple[str, str | list[str]]:
            try:
                transform_kwargs = dict(kwargs)
                res = diskfile.format_transform(
                    format,
                    frame=frame,
                    embed_in_one_file=embed_in_one_file,
                    write_to_disk=write_to_disk,
                    file_path=os.path.join(
                        output_dir,
                        f"{os.path.splitext(diskfile.filename)[0]}.{output_extension}",
                    )
                    if write_to_disk and output_dir
                    else None,
                    **transform_kwargs,
                )
                return diskfile.file_path, res
            except NativeReconstructionConcurrencyError:
                # A missed parent prewarm is a scheduling violation, not a
                # per-file rendering failure. Do not turn it into empty output.
                raise
            except Exception as e:
                moloplogger.warning(
                    f"Format transform failed for {diskfile.filename}: {type(e).__name__}: {e}"
                )
                return diskfile.file_path, ("" if embed_in_one_file else [])

        desc = f"MolOP processing {format} format with {effective_jobs} jobs"
        results = typed_self.parallel_execute(
            transform_func,
            desc,
            effective_jobs,
            return_results=True,
        )
        return dict(results)
