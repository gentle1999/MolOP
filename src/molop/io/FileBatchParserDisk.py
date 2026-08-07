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
from typing import Any, Literal, Protocol, cast, overload

from joblib import Parallel, delayed

from molop.config import molopconfig, moloplogger
from molop.io.codec_exceptions import FormatMismatchError
from molop.io.codec_types import ParseOptions, ParseResult
from molop.io.FileBatchModelDisk import FileBatchModelDisk, FileDiskObj, _looks_like_disk_file
from molop.io.parse_outcomes import BatchParseResult, FileParseOutcome, ParseFailure
from molop.utils.progressbar import AdaptiveProgress

from . import codec_registry


class _FileReaderCodec(Protocol):
    format_id: str

    def read(self, path: str | pathlib.Path, **kwargs: Any) -> ParseResult[Any]: ...


@overload
def single_file_parser(
    file_path: str,
    possible_readers: tuple[_FileReaderCodec, ...],
    total_charge: int | None = None,
    total_multiplicity: int | None = None,
    only_extract_structure: bool = False,
    only_last_frame: bool = False,
    capture_source_evidence: bool = False,
    source_encoding: str = "utf-8",
    release_file_content: bool = False,
    *,
    input_index: int | None = None,
    parse_options: ParseOptions | None = None,
    return_outcome: Literal[False] = False,
) -> FileDiskObj | None: ...


@overload
def single_file_parser(
    file_path: str,
    possible_readers: tuple[_FileReaderCodec, ...],
    total_charge: int | None = None,
    total_multiplicity: int | None = None,
    only_extract_structure: bool = False,
    only_last_frame: bool = False,
    capture_source_evidence: bool = False,
    source_encoding: str = "utf-8",
    release_file_content: bool = False,
    *,
    input_index: int | None = None,
    parse_options: ParseOptions | None = None,
    return_outcome: Literal[True],
) -> FileParseOutcome[FileDiskObj]: ...


def single_file_parser(
    file_path: str,
    possible_readers: tuple[_FileReaderCodec, ...],
    total_charge: int | None = None,
    total_multiplicity: int | None = None,
    only_extract_structure: bool = False,
    only_last_frame: bool = False,
    capture_source_evidence: bool = False,
    source_encoding: str = "utf-8",
    release_file_content: bool = False,
    *,
    input_index: int | None = None,
    parse_options: ParseOptions | None = None,
    return_outcome: bool = False,
) -> FileDiskObj | FileParseOutcome[FileDiskObj] | None:
    outcome = _single_file_parse_outcome(
        file_path=file_path,
        possible_readers=possible_readers,
        total_charge=total_charge,
        total_multiplicity=total_multiplicity,
        only_extract_structure=only_extract_structure,
        only_last_frame=only_last_frame,
        capture_source_evidence=capture_source_evidence,
        source_encoding=source_encoding,
        release_file_content=release_file_content,
        input_index=input_index,
        parse_options=parse_options,
    )
    if return_outcome:
        return outcome
    return outcome.value if outcome.succeeded else None


def _single_file_parse_outcome(
    file_path: str,
    possible_readers: tuple[_FileReaderCodec, ...],
    total_charge: int | None = None,
    total_multiplicity: int | None = None,
    only_extract_structure: bool = False,
    only_last_frame: bool = False,
    capture_source_evidence: bool = False,
    source_encoding: str = "utf-8",
    release_file_content: bool = False,
    input_index: int | None = None,
    parse_options: ParseOptions | None = None,
) -> FileParseOutcome[FileDiskObj]:
    options = (
        parse_options
        or ParseOptions(
            total_charge=total_charge,
            total_multiplicity=total_multiplicity,
            only_extract_structure=only_extract_structure,
            only_last_frame=only_last_frame,
            capture_source_evidence=capture_source_evidence,
            source_encoding=source_encoding,
            release_file_content=release_file_content,
        )
    ).resolved()
    for idx, reader in enumerate(possible_readers):
        reader_name = getattr(reader, "format_id", reader.__class__.__name__)
        try:
            result = reader.read(
                file_path,
                total_charge=options.total_charge,
                total_multiplicity=options.total_multiplicity,
                only_extract_structure=options.only_extract_structure,
                only_last_frame=options.only_last_frame,
                capture_source_evidence=options.capture_source_evidence,
                source_encoding=options.source_encoding,
                release_file_content=options.release_file_content,
                parse_options=options,
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
            disk_file = cast(FileDiskObj, value)
            detected = result.detected_format
            normalized_detected = detected.strip().lower() if detected else None
            if len(disk_file) == 0:
                return FileParseOutcome(
                    file_path=file_path,
                    status="empty",
                    value=disk_file,
                    warnings=result.warnings,
                    detected_format=normalized_detected,
                    failure=ParseFailure(
                        kind="empty_result",
                        message="Reader returned a file model with no frames.",
                        reader_format=reader_name,
                    ),
                    input_index=input_index,
                )
            return FileParseOutcome(
                file_path=file_path,
                status="ok",
                value=disk_file,
                warnings=result.warnings,
                detected_format=normalized_detected,
                input_index=input_index,
            )
        except FormatMismatchError as e:
            if idx == len(possible_readers) - 1:
                moloplogger.error(f"Failed to parse file {file_path} with {reader_name}. {e}")
                return FileParseOutcome(
                    file_path=file_path,
                    status="mismatch",
                    failure=ParseFailure(
                        kind="format_mismatch",
                        message=str(e),
                        reader_format=reader_name,
                        exception_type=type(e).__name__,
                    ),
                    input_index=input_index,
                )
            moloplogger.debug(
                f"Failed to parse file {file_path} with {reader_name}, "
                f"trying {getattr(possible_readers[idx + 1], 'format_id', possible_readers[idx + 1].__class__.__name__)} "
                "instead"
            )
        except Exception as e:
            moloplogger.error(f"Failed to parse file {file_path} with {reader_name}. {e}")
            return FileParseOutcome(
                file_path=file_path,
                status="error",
                failure=ParseFailure(
                    kind="parse_error",
                    message=str(e),
                    reader_format=reader_name,
                    exception_type=type(e).__name__,
                ),
                input_index=input_index,
            )
    return FileParseOutcome(
        file_path=file_path,
        status="unsupported",
        failure=ParseFailure(
            kind="no_reader",
            message="No reader candidates were supplied.",
        ),
        input_index=input_index,
    )


def _run_single_file_task(**task: Any) -> FileParseOutcome[FileDiskObj]:
    """Normalize legacy or monkeypatched worker returns to the outcome contract."""

    result = single_file_parser(**task)
    if isinstance(result, FileParseOutcome):
        return result
    file_path = cast(str, task["file_path"])
    input_index = cast(int | None, task.get("input_index"))
    if result is not None and _looks_like_disk_file(result):
        disk_file = cast(FileDiskObj, result)
        return FileParseOutcome(
            file_path=file_path,
            status="ok" if len(disk_file) > 0 else "empty",
            value=disk_file,
            detected_format=getattr(disk_file, "detected_format_id", None),
            input_index=input_index,
        )
    return FileParseOutcome(
        file_path=file_path,
        status="error",
        failure=ParseFailure(
            kind="legacy_none_result",
            message="Parser task returned no structured outcome.",
        ),
        input_index=input_index,
    )


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
        capture_source_evidence: bool = False,
        source_encoding: str = "utf-8",
        release_file_content: bool = True,
        parser_detection: str = "auto",
        parse_options: ParseOptions | None = None,
    ) -> FileBatchModelDisk[FileDiskObj]:
        """Parse inputs and return the backward-compatible successful-file batch."""

        return self.parse_with_report(
            file_paths,
            total_charge=total_charge,
            total_multiplicity=total_multiplicity,
            only_extract_structure=only_extract_structure,
            only_last_frame=only_last_frame,
            capture_source_evidence=capture_source_evidence,
            source_encoding=source_encoding,
            release_file_content=release_file_content,
            parser_detection=parser_detection,
            parse_options=parse_options,
        ).batch

    def parse_with_report(
        self,
        file_paths: Iterable[str] | Iterable[pathlib.Path],
        total_charge: int | None = None,
        total_multiplicity: int | None = None,
        only_extract_structure: bool = False,
        only_last_frame: bool = False,
        capture_source_evidence: bool = False,
        source_encoding: str = "utf-8",
        release_file_content: bool = True,
        parser_detection: str = "auto",
        parse_options: ParseOptions | None = None,
    ) -> BatchParseResult[FileBatchModelDisk[FileDiskObj], FileDiskObj]:
        """
        Parse input files and retain one structured outcome for every input path.

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
            capture_source_evidence (bool):
                if True, capture format-supported source evidence while parsing.
            source_encoding (str):
                Text encoding used for strict source decoding and byte offsets.
            parser_detection (str):
                if "auto", use the file extension to detect the parser, else use the given format id.
        """
        hint_format = None if parser_detection == "auto" else parser_detection
        path_count = _length_or_none(file_paths)
        parse_options = (
            parse_options
            or ParseOptions(
                total_charge=total_charge,
                total_multiplicity=total_multiplicity,
                only_extract_structure=only_extract_structure,
                only_last_frame=only_last_frame,
                capture_source_evidence=capture_source_evidence,
                source_encoding=source_encoding,
                release_file_content=release_file_content,
            )
        ).resolved()

        preflight_outcomes: list[FileParseOutcome[FileDiskObj]] = []

        def process_path(file_path: str | pathlib.Path, input_index: int):
            if isinstance(file_path, pathlib.Path):
                file_path = file_path.as_posix()
            if not os.path.isfile(file_path):
                moloplogger.warning(f"{file_path} is not a file.")
                return FileParseOutcome(
                    file_path=os.path.abspath(file_path),
                    status="missing",
                    failure=ParseFailure(
                        kind="missing_file",
                        message="Input path is not an existing file.",
                    ),
                    input_index=input_index,
                )
            if file_path.endswith("molop.log"):
                return FileParseOutcome(
                    file_path=os.path.abspath(file_path),
                    status="skipped",
                    failure=ParseFailure(
                        kind="internal_log",
                        message="MolOP's own log file is excluded from parsing.",
                    ),
                    input_index=input_index,
                )
            abs_path = os.path.abspath(file_path)
            try:
                possible_readers = cast(
                    tuple[_FileReaderCodec, ...],
                    codec_registry.select_reader(abs_path, hint_format=hint_format),
                )
            except codec_registry.UnsupportedFormatError:
                if parser_detection == "auto":
                    moloplogger.warning(f"Unsupported input file format: {abs_path}")
                    return FileParseOutcome(
                        file_path=abs_path,
                        status="unsupported",
                        failure=ParseFailure(
                            kind="unsupported_format",
                            message="No reader codec is registered for this input path.",
                        ),
                        input_index=input_index,
                    )
                moloplogger.error(f"Unsupported input file format: {abs_path}")
                return FileParseOutcome(
                    file_path=abs_path,
                    status="unsupported",
                    failure=ParseFailure(
                        kind="unsupported_format",
                        message=f"No reader codec is registered for format {parser_detection!r}.",
                    ),
                    input_index=input_index,
                )
            if parser_detection == "auto":
                possible_readers = _filter_readers_by_probe(abs_path, possible_readers)
            return {
                "file_path": abs_path,
                "possible_readers": possible_readers,
                "total_charge": total_charge,
                "total_multiplicity": total_multiplicity,
                "only_extract_structure": only_extract_structure,
                "only_last_frame": only_last_frame,
                "capture_source_evidence": capture_source_evidence,
                "source_encoding": source_encoding,
                "release_file_content": release_file_content,
                "parse_options": parse_options,
                "input_index": input_index,
                "return_outcome": True,
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
            for input_index, file_path in enumerate(file_paths):
                register_input_path()
                task_or_outcome = process_path(
                    cast(str | pathlib.Path, file_path),
                    input_index,
                )
                if isinstance(task_or_outcome, FileParseOutcome):
                    preflight_outcomes.append(task_or_outcome)
                    mark_input_path_done()
                else:
                    yield task_or_outcome

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
                )(delayed(_run_single_file_task)(**task) for task in iter_tasks())
            else:
                results = (_run_single_file_task(**task) for task in iter_tasks())

            parsed_diskfiles: list[FileDiskObj] = []
            parse_outcomes = preflight_outcomes
            for outcome in cast(Iterable[FileParseOutcome[FileDiskObj]], results):
                parse_outcomes.append(outcome)
                if outcome.succeeded and outcome.value is not None:
                    parsed_diskfiles.append(outcome.value)
                mark_input_path_done()

            close_progress()
            parsed_diskfiles.sort(key=lambda diskfile: diskfile.file_path)
            parse_outcomes.sort(
                key=lambda outcome: (
                    outcome.input_index is None,
                    outcome.input_index if outcome.input_index is not None else 0,
                    outcome.file_path,
                )
            )
            batch = FileBatchModelDisk._new_batch_from_sorted_diskfiles(parsed_diskfiles)
            return BatchParseResult(batch=batch, outcomes=tuple(parse_outcomes))
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
