"""
Author: TMJ
Date: 2025-08-20 22:55:18
LastEditors: TMJ
LastEditTime: 2026-02-11 22:20:34
Description: 请填写简介
"""

from __future__ import annotations

import os
import pathlib
from collections.abc import Iterable
from contextlib import suppress
from typing import Any, Protocol, cast

from joblib import Parallel, delayed

from molop.config import molopconfig, moloplogger
from molop.io.codec_types import ParseResult
from molop.io.FileBatchModelDisk import FileBatchModelDisk, FileDiskObj, _looks_like_disk_file
from molop.utils.progressbar import AdaptiveProgress

from . import codec_registry


class _FileReaderCodec(Protocol):
    format_id: str

    def read(self, path: str | pathlib.Path, **kwargs: Any) -> ParseResult[Any]: ...


def single_file_parser(
    file_path: str,
    possible_readers: tuple[_FileReaderCodec, ...],
    total_charge: int | None = None,
    total_multiplicity: int | None = None,
    only_extract_structure: bool = False,
    only_last_frame: bool = False,
    release_file_content: bool = False,
) -> FileDiskObj | None:
    for idx, reader in enumerate(possible_readers):
        try:
            result = reader.read(
                file_path,
                total_charge=total_charge,
                total_multiplicity=total_multiplicity,
                only_extract_structure=only_extract_structure,
                only_last_frame=only_last_frame,
                release_file_content=release_file_content,
            )
            value = result.value
            if hasattr(value, "detected_format_id"):
                detected = result.detected_format
                with suppress(Exception):
                    value.detected_format_id = detected.strip().lower() if detected else None
            if not _looks_like_disk_file(value):
                raise TypeError(
                    f"Reader {getattr(reader, 'format_id', reader.__class__.__name__)} returned "
                    f"unexpected value type: {type(value)}"
                )
            return cast(FileDiskObj, value)
        except Exception as e:
            reader_name = getattr(reader, "format_id", reader.__class__.__name__)
            if idx == len(possible_readers) - 1:
                moloplogger.error(f"Failed to parse file {file_path} with {reader_name}. {e}")
                return None
            moloplogger.debug(
                f"Failed to parse file {file_path} with {reader_name}, "
                f"trying {getattr(possible_readers[idx + 1], 'format_id', possible_readers[idx + 1].__class__.__name__)} "
                "instead"
            )
    return None


class FileBatchParserDisk:
    __n_jobs: int

    def __init__(self, n_jobs: int = -1):
        self.__n_jobs = molopconfig.set_n_jobs(n_jobs)

    def parse(
        self,
        file_paths: Iterable[str] | Iterable[pathlib.Path],
        total_charge: int | None = None,
        total_multiplicity: int | None = None,
        only_extract_structure: bool = False,
        only_last_frame: bool = False,
        release_file_content: bool = True,
        parser_detection: str = "auto",
    ) -> FileBatchModelDisk[FileDiskObj]:
        """
        Parses a list of input files and returns a FileBatchModelDisk object.

        Parameters:
            file_paths (Iterable[str]):
                A list of wildcard of input file paths.
            total_charge (int | None):
                forced charge of the molecule, if not given, will use the charge written in the file or 0.
            total_multiplicity (int | None):
                forced multiplicity of the molecule, if not given, will use the charge written in the file or 1.
            only_extract_structure (bool):
                if True, only extract the structure, else extract the whole file.
            only_last_frame (bool):
                if True, only extract the last frame, else extract all frames.
            parser_detection (str):
                if "auto", use the file extension to detect the parser, else use the given format id.
        """
        tasks: list[dict[str, Any]] = []
        hint_format = None if parser_detection == "auto" else parser_detection
        for file_path in file_paths:
            if isinstance(file_path, pathlib.Path):
                file_path = file_path.as_posix()
            if not os.path.isfile(file_path):
                moloplogger.warning(f"{file_path} is not a file.")
                continue
            if file_path.endswith("molop.log"):
                continue
            abs_path = os.path.abspath(file_path)
            try:
                possible_readers = cast(
                    tuple[_FileReaderCodec, ...],
                    codec_registry.select_reader(abs_path, hint_format=hint_format),
                )
            except codec_registry.UnsupportedFormatError:
                if parser_detection == "auto":
                    moloplogger.warning(f"Unsupported input file format: {abs_path}")
                    continue
                moloplogger.error(f"Unsupported input file format: {abs_path}")
                continue
            tasks.append(
                {
                    "file_path": abs_path,
                    "possible_readers": possible_readers,
                    "total_charge": total_charge,
                    "total_multiplicity": total_multiplicity,
                    "only_extract_structure": only_extract_structure,
                    "only_last_frame": only_last_frame,
                    "release_file_content": release_file_content,
                }
            )
        if len(tasks) == 0:
            moloplogger.error("No valid input files.")
            return FileBatchModelDisk()
        try:
            tasks.sort(key=lambda task: os.path.getsize(task["file_path"]), reverse=True)
        except OSError as e:
            moloplogger.warning(f"Could not get file size for sorting, proceeding without it: {e}")
        total_tasks = tasks
        # Determine if parallel processing should be used
        use_parallel = self.__n_jobs > 1 and len(total_tasks) > self.__n_jobs

        if use_parallel:
            moloplogger.info(f"Using {self.__n_jobs} processes for {len(total_tasks)} tasks.")
            results = Parallel(
                n_jobs=self.__n_jobs,
                maxtasks_per_child=50,
                return_as="generator",
                max_nbytes=molopconfig.parallel_max_size,
            )(
                delayed(single_file_parser)(**task)
                for task in AdaptiveProgress(
                    total_tasks,
                    total=len(total_tasks),
                    desc=f"MolOP parsing with {self.__n_jobs} processes",
                    disable=not molopconfig.show_progress_bar,
                )
            )
        else:
            results = (
                single_file_parser(**task)
                for task in AdaptiveProgress(
                    total_tasks,
                    total=len(total_tasks),
                    desc="MolOP parsing with single process",
                    disable=not molopconfig.show_progress_bar,
                )
            )
        typed_results = cast(Iterable[FileDiskObj | None], results)
        return FileBatchModelDisk.new_batch(
            result for result in typed_results if result is not None and len(result) > 0
        )

    @property
    def n_jobs(self) -> int:
        return self.__n_jobs

    @n_jobs.setter
    def n_jobs(self, n_jobs: int) -> None:
        self.__n_jobs = molopconfig.set_n_jobs(n_jobs)

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(n_jobs={self.__n_jobs})"
