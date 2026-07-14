"""
Author: TMJ
Date: 2025-01-15 23:01:22
LastEditors: TMJ
LastEditTime: 2026-04-01 14:02:46
Description: 请填写简介
"""

import glob
import os
from collections.abc import Iterable
from pathlib import Path
from typing import TYPE_CHECKING, TypeAlias

from molop.io.FileBatchModelDisk import FileBatchModelDisk
from molop.io.FileBatchParserDisk import FileBatchParserDisk


if TYPE_CHECKING:
    from molop.io.FileBatchModelDisk import FileDiskObj


PathSpec: TypeAlias = str | os.PathLike[str]
PathInput: TypeAlias = PathSpec | Iterable[PathSpec]


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
) -> FileBatchModelDisk["FileDiskObj"]:
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
    return FileBatchParserDisk(n_jobs=n_jobs).parse(
        files,
        total_charge=total_charge,
        total_multiplicity=total_multiplicity,
        only_extract_structure=only_extract_structure,
        only_last_frame=only_last_frame,
        capture_source_evidence=capture_source_evidence,
        source_encoding=source_encoding,
        release_file_content=release_file_content,
        parser_detection=parser_detection,
    )
