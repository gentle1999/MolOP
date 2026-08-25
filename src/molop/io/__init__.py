"""
Author: TMJ
Date: 2025-01-15 23:01:22
LastEditors: TMJ
LastEditTime: 2026-08-22 14:06:41
Description: 请填写简介
"""

import glob
import os
from collections.abc import Iterable
from pathlib import Path
from typing import TYPE_CHECKING, Literal, TypeAlias, cast, overload

from molop.io import codec_registry
from molop.io.codec_exceptions import FormatMismatchError, ParseError, UnsupportedFormatError
from molop.io.codec_types import ParseOptions
from molop.io.FileBatchModelDisk import FileBatchModelDisk
from molop.io.FileBatchParserDisk import (
    FileBatchParserDisk,
    _FileReaderCodec,
    _filter_readers_by_probe,
)
from molop.io.FileBatchParserDisk import (
    single_file_parser as _single_file_parser,
)
from molop.io.memory_parser import (
    AutoBytesParser,
    AutoMemoryParser,
    AutoParserMemory,
    AutoTextParser,
)
from molop.io.parse_outcomes import BatchParseResult


if TYPE_CHECKING:
    from molop.io.FileBatchModelDisk import FileDiskObj


PathSpec: TypeAlias = str | os.PathLike[str]
PathInput: TypeAlias = PathSpec | Iterable[PathSpec]


def AutoFileParser(
    file_path: PathSpec,
    *,
    total_charge: int | None = None,
    total_multiplicity: int | None = None,
    only_extract_structure: bool = False,
    only_last_frame: bool = False,
    capture_source_evidence: bool = False,
    source_encoding: str = "utf-8",
    release_file_content: bool = True,
    parser_detection: str = "auto",
    parse_options: ParseOptions | None = None,
) -> "FileDiskObj":
    """Parse one file with automatic format detection.

    Unlike :func:`AutoParser`, this entrypoint accepts exactly one file path and
    returns the parsed file model directly. It performs no glob expansion,
    batch preparation, or process scheduling, so callers can place the parse
    inside their own worker or process model.

    ``parser_detection="auto"`` first uses the path extension and reader
    probes to narrow the candidates, then lets the parser's format checks
    validate the complete source. Pass a format id to select a reader
    explicitly while retaining the same file-level return type.

    Raises:
        TypeError: If ``file_path`` is not a text path-like object.
        ValueError: If ``file_path`` is not an existing regular file.
        UnsupportedFormatError: If no reader is registered for the input.
        FormatMismatchError: If all candidate readers reject the input format.
        ParseError: If a selected reader fails during parsing.
    """
    try:
        path_value = os.fspath(file_path)
    except TypeError as exc:
        raise TypeError("file_path must be a path-like object") from exc
    if not isinstance(path_value, str):
        raise TypeError("file_path must resolve to a text path, not bytes")

    path = Path(os.path.abspath(path_value))
    if not path.is_file():
        raise ValueError(f"File {path} does not exist or is not a regular file.")

    hint_format = None if parser_detection == "auto" else parser_detection
    possible_readers = cast(
        tuple[_FileReaderCodec, ...],
        codec_registry.select_reader(path, hint_format=hint_format),
    )
    if parser_detection == "auto":
        possible_readers = _filter_readers_by_probe(str(path), possible_readers)

    outcome = _single_file_parser(
        file_path=str(path),
        possible_readers=possible_readers,
        total_charge=total_charge,
        total_multiplicity=total_multiplicity,
        only_extract_structure=only_extract_structure,
        only_last_frame=only_last_frame,
        capture_source_evidence=capture_source_evidence,
        source_encoding=source_encoding,
        release_file_content=release_file_content,
        parse_options=parse_options,
        return_outcome=True,
    )
    if outcome.succeeded and outcome.value is not None:
        return outcome.value

    failure = outcome.failure
    message = failure.message if failure is not None else f"Failed to parse file {path}."
    if outcome.status == "unsupported":
        raise UnsupportedFormatError(message)
    if outcome.status == "mismatch":
        raise FormatMismatchError(message)
    raise ParseError(message)


def split_path_pattern(path_str: str) -> tuple[Path, str]:
    p = Path(path_str)
    parts = p.parts

    split_index = len(parts)
    for i, part in enumerate(parts):
        if glob.has_magic(part):
            split_index = i
            break
    base_path = Path(*parts[:split_index])
    pattern = str(Path(*parts[split_index:])) if split_index < len(parts) else ""

    return base_path, pattern


def _expand_path_spec(path_spec: PathSpec) -> list[Path]:
    path_str = os.fspath(path_spec)
    if not isinstance(path_str, str):
        raise TypeError("path-like objects must resolve to str, not bytes")
    if os.path.isfile(path_str):
        return [Path(path_str)]
    if not glob.has_magic(path_str):
        return [Path(path_str)]
    base_path, pattern = split_path_pattern(path_str)
    return sorted(base_path.glob(pattern))


def _normalize_file_paths(file_path: PathInput) -> list[Path]:
    if isinstance(file_path, (str, os.PathLike)):
        path_specs: Iterable[PathSpec] = (file_path,)
    else:
        try:
            path_specs = iter(file_path)
        except TypeError as exc:
            raise TypeError(
                "file_path must be a path-like object or an iterable of path-like objects"
            ) from exc

    paths_by_key: dict[str, Path] = {}
    for index, path_spec in enumerate(path_specs):
        if not isinstance(path_spec, (str, os.PathLike)):
            raise TypeError(
                f"file_path[{index}] must be a path-like object, got {type(path_spec).__name__}"
            )
        try:
            expanded_paths = _expand_path_spec(path_spec)
        except TypeError as exc:
            raise TypeError(
                f"file_path[{index}] must resolve to a text path, got {type(path_spec).__name__}"
            ) from exc
        for path in expanded_paths:
            absolute_path = Path(os.path.abspath(path))
            key = os.path.normcase(os.fspath(absolute_path))
            paths_by_key.setdefault(key, absolute_path)

    return [paths_by_key[key] for key in sorted(paths_by_key)]


@overload
def AutoParser(
    file_path: PathInput,
    *,
    total_charge: int | None = None,
    total_multiplicity: int | None = None,
    n_jobs: int = -1,
    only_extract_structure: bool = False,
    only_last_frame: bool = False,
    capture_source_evidence: bool = False,
    source_encoding: str = "utf-8",
    release_file_content: bool = True,
    parser_detection: str = "auto",
    parse_options: ParseOptions | None = None,
    return_report: Literal[False] = False,
) -> FileBatchModelDisk["FileDiskObj"]: ...


@overload
def AutoParser(
    file_path: PathInput,
    *,
    total_charge: int | None = None,
    total_multiplicity: int | None = None,
    n_jobs: int = -1,
    only_extract_structure: bool = False,
    only_last_frame: bool = False,
    capture_source_evidence: bool = False,
    source_encoding: str = "utf-8",
    release_file_content: bool = True,
    parser_detection: str = "auto",
    parse_options: ParseOptions | None = None,
    return_report: Literal[True],
) -> BatchParseResult[FileBatchModelDisk["FileDiskObj"], "FileDiskObj"]: ...


def AutoParser(
    file_path: PathInput,
    *,
    total_charge: int | None = None,
    total_multiplicity: int | None = None,
    n_jobs: int = -1,
    only_extract_structure=False,
    only_last_frame=False,
    capture_source_evidence: bool = False,
    source_encoding: str = "utf-8",
    release_file_content: bool = True,
    parser_detection: str = "auto",
    parse_options: ParseOptions | None = None,
    return_report: bool = False,
) -> (
    FileBatchModelDisk["FileDiskObj"]
    | BatchParseResult[FileBatchModelDisk["FileDiskObj"], "FileDiskObj"]
):
    """
    The Entrypoint of MolOP

    Parameters:
        file_path (PathInput):
            A path, glob pattern, or iterable containing paths and glob patterns.
        total_charge (int | None):
            forced charge of the molecule, if not given, will use the charge written in the file or 0.
        total_multiplicity (int | None):
            forced multiplicity of the molecule, if not given, will use the charge written in the file or 1.
        n_jobs (int):
            number of jobs to use, if -1, use cpu with default max number.
        only_extract_structure (bool):
            if True, only extract the structure, else extract the whole file.
        only_last_frame (bool):
            if True, only extract the last frame, else extract all frames.
        capture_source_evidence (bool):
            if True, capture format-supported source evidence while parsing.
        source_encoding (str):
            Text encoding used for strict source decoding and byte offsets.
        release_file_content (bool):
            if True, release the file content after parsing, else keep the file content in memory.
        parser_detection (str):
            if "auto", use the file extension to detect the parser, else use the given format id.
        parse_options (ParseOptions | None):
            immutable parsing options; when provided, these take precedence over individual options.
        return_report (bool):
            if True, return the successful batch together with one structured outcome per input.

    Returns:
        FileBatchModelDisk: Parsed files sorted by absolute file path.

    How to use
    ----------
    ```python
    from molop import AutoParser

    parser = AutoParser("/path/to/file")  # parse one file
    parser = AutoParser("/path/to/files/*.log")  # parse one glob pattern
    parser = AutoParser(["a.log", "group_b/*.log"])  # parse multiple inputs
    ```
    """
    files = _normalize_file_paths(file_path)
    parser = FileBatchParserDisk(n_jobs=n_jobs)
    if return_report:
        return parser.parse_with_report(
            files,
            total_charge=total_charge,
            total_multiplicity=total_multiplicity,
            only_extract_structure=only_extract_structure,
            only_last_frame=only_last_frame,
            capture_source_evidence=capture_source_evidence,
            source_encoding=source_encoding,
            release_file_content=release_file_content,
            parser_detection=parser_detection,
            parse_options=parse_options,
        )
    return parser.parse(
        files,
        total_charge=total_charge,
        total_multiplicity=total_multiplicity,
        only_extract_structure=only_extract_structure,
        only_last_frame=only_last_frame,
        capture_source_evidence=capture_source_evidence,
        source_encoding=source_encoding,
        release_file_content=release_file_content,
        parser_detection=parser_detection,
        parse_options=parse_options,
    )


__all__ = [
    "AutoBytesParser",
    "AutoFileParser",
    "AutoMemoryParser",
    "AutoParser",
    "AutoParserMemory",
    "AutoTextParser",
    "FileBatchModelDisk",
    "FileBatchParserDisk",
    "split_path_pattern",
]
