"""
Author: TMJ
Date: 2025-08-20 22:55:18
LastEditors: TMJ
LastEditTime: 2026-05-02 18:24:16
Description: 请填写简介
"""

from __future__ import annotations

import os
import pathlib
from collections.abc import Iterable, Sequence, Sized
from contextlib import suppress
from typing import Any, Protocol, cast

from joblib import Parallel, delayed

from molop.config import molopconfig, moloplogger
from molop.io.codec_exceptions import FormatMismatchError
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
        except FormatMismatchError as e:
            reader_name = getattr(reader, "format_id", reader.__class__.__name__)
            if idx == len(possible_readers) - 1:
                moloplogger.error(f"Failed to parse file {file_path} with {reader_name}. {e}")
                return None
            moloplogger.debug(
                f"Failed to parse file {file_path} with {reader_name}, "
                f"trying {getattr(possible_readers[idx + 1], 'format_id', possible_readers[idx + 1].__class__.__name__)} "
                "instead"
            )
        except Exception as e:
            reader_name = getattr(reader, "format_id", reader.__class__.__name__)
            moloplogger.error(f"Failed to parse file {file_path} with {reader_name}. {e}")
            return None
    return None


def _length_or_none(value: object) -> int | None:
    if isinstance(value, Sized):
        return len(value)
    return None


def _safe_get_file_size(file_path: str | pathlib.Path) -> int | None:
    try:
        return os.path.getsize(file_path)
    except OSError:
        return None


def _known_paths_are_large_enough_for_processes(
    file_paths: Sequence[str] | Sequence[pathlib.Path],
) -> bool:
    total_size = 0
    for file_path in file_paths:
        size = _safe_get_file_size(file_path)
        if size is None:
            return True
        total_size += size
    return total_size >= molopconfig.parallel_max_size


def _tune_effective_jobs_for_known_paths(
    file_paths: Iterable[str] | Iterable[pathlib.Path],
    path_count: int | None,
    effective_jobs: int,
) -> int:
    if effective_jobs <= 1:
        return effective_jobs
    if path_count is None or path_count <= 1:
        return effective_jobs
    if not isinstance(file_paths, Sequence):
        return effective_jobs
    if path_count > effective_jobs * 2:
        return effective_jobs
    if not _known_paths_are_large_enough_for_processes(file_paths):
        return 1
    return min(effective_jobs, max(2, path_count // 2))


def _should_use_parallel(
    path_count: int | None,
    effective_jobs: int,
) -> bool:
    if effective_jobs <= 1:
        return False
    if path_count is None:
        return True
    return path_count > 1


def _reader_matches_file(reader: _FileReaderCodec, file_path: str) -> bool:
    probe = getattr(reader, "probe_file_format", None)
    if not callable(probe):
        return True
    try:
        return bool(probe(file_path))
    except Exception as exc:
        reader_name = getattr(reader, "format_id", reader.__class__.__name__)
        moloplogger.debug(
            f"Reader {reader_name} probe failed for {file_path}; keeping it as a candidate. {exc}"
        )
        return True


def _filter_readers_by_probe(
    file_path: str,
    readers: tuple[_FileReaderCodec, ...],
) -> tuple[_FileReaderCodec, ...]:
    probed = tuple(reader for reader in readers if _reader_matches_file(reader, file_path))
    return probed or readers


def _task_sort_size(task: dict[str, Any]) -> int:
    file_path = task.get("file_path")
    if not isinstance(file_path, str):
        return 0
    return _safe_get_file_size(file_path) or 0


def _iter_size_ordered_task_buffer(
    tasks: Iterable[dict[str, Any]],
    *,
    buffer_size: int,
) -> Iterable[dict[str, Any]]:
    buffer: list[dict[str, Any]] = []
    for task in tasks:
        buffer.append(task)
        if len(buffer) >= buffer_size:
            yield from sorted(buffer, key=_task_sort_size, reverse=True)
            buffer.clear()
    if buffer:
        yield from sorted(buffer, key=_task_sort_size, reverse=True)


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
        hint_format = None if parser_detection == "auto" else parser_detection
        path_count = _length_or_none(file_paths)

        def process_path(file_path: str | pathlib.Path):
            if isinstance(file_path, pathlib.Path):
                file_path = file_path.as_posix()
            if not os.path.isfile(file_path):
                moloplogger.warning(f"{file_path} is not a file.")
                return None
            if file_path.endswith("molop.log"):
                return None
            abs_path = os.path.abspath(file_path)
            try:
                possible_readers = cast(
                    tuple[_FileReaderCodec, ...],
                    codec_registry.select_reader(abs_path, hint_format=hint_format),
                )
            except codec_registry.UnsupportedFormatError:
                if parser_detection == "auto":
                    moloplogger.warning(f"Unsupported input file format: {abs_path}")
                    return None
                moloplogger.error(f"Unsupported input file format: {abs_path}")
                return None
            if parser_detection == "auto":
                possible_readers = _filter_readers_by_probe(abs_path, possible_readers)
            return {
                "file_path": abs_path,
                "possible_readers": possible_readers,
                "total_charge": total_charge,
                "total_multiplicity": total_multiplicity,
                "only_extract_structure": only_extract_structure,
                "only_last_frame": only_last_frame,
                "release_file_content": release_file_content,
            }

        effective_jobs = (
            self.__n_jobs if path_count is None else min(self.__n_jobs, max(path_count, 1))
        )
        effective_jobs = _tune_effective_jobs_for_known_paths(
            file_paths,
            path_count,
            effective_jobs,
        )
        use_parallel = _should_use_parallel(path_count, effective_jobs)
        desc = (
            f"MolOP parsing with {effective_jobs} processes"
            if use_parallel
            else "MolOP parsing with single process"
        )
        if use_parallel:
            task_label = path_count if path_count is not None else "streamed"
            moloplogger.info(f"Using {effective_jobs} processes for {task_label} input paths.")
        progress = AdaptiveProgress(
            None,
            disable=not molopconfig.show_progress_bar,
            desc=desc,
            total=path_count,
        )

        def register_input_path() -> None:
            if path_count is not None:
                return
            with suppress(Exception):
                progress.total = (progress.total or 0) + 1
                progress.refresh()

        def mark_input_path_done() -> None:
            with suppress(Exception):
                progress.update(1)

        def close_progress() -> None:
            with suppress(Exception):
                progress.close()

        def iter_raw_tasks():
            for file_path in file_paths:
                register_input_path()
                task = process_path(file_path)
                if task is not None:
                    yield task
                else:
                    mark_input_path_done()

        def iter_tasks():
            tasks = iter_raw_tasks()
            if not use_parallel:
                yield from tasks
                return
            yield from _iter_size_ordered_task_buffer(
                tasks,
                buffer_size=max(min(effective_jobs * 16, 256), 1),
            )

        try:
            if use_parallel:
                results = Parallel(
                    n_jobs=effective_jobs,
                    maxtasks_per_child=50,
                    return_as="generator_unordered",
                    max_nbytes=molopconfig.parallel_max_size,
                )(delayed(single_file_parser)(**task) for task in iter_tasks())
            else:
                results = (single_file_parser(**task) for task in iter_tasks())

            parsed_diskfiles: list[FileDiskObj] = []
            for result in cast(Iterable[FileDiskObj | None], results):
                if result is not None and len(result) > 0:
                    parsed_diskfiles.append(result)
                mark_input_path_done()

            close_progress()
            parsed_diskfiles.sort(key=lambda diskfile: diskfile.file_path)
            return FileBatchModelDisk._new_batch_from_sorted_diskfiles(parsed_diskfiles)
        except Exception:
            close_progress()
            raise

    @property
    def n_jobs(self) -> int:
        return self.__n_jobs

    @n_jobs.setter
    def n_jobs(self, n_jobs: int) -> None:
        self.__n_jobs = molopconfig.set_n_jobs(n_jobs)

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(n_jobs={self.__n_jobs})"
