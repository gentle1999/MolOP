"""
Author: TMJ
Date: 2025-08-20 22:55:18
LastEditors: TMJ
LastEditTime: 2025-12-16 01:31:00
Description: 请填写简介
"""

import os
import pathlib
from collections.abc import Iterable
from enum import Enum
from typing import Any, Literal

from joblib import Parallel, delayed

from molop.config import molopconfig, moloplogger
from molop.io.FileBatchModelDisk import FileBatchModelDisk
from molop.io.types import (
    FILEDISK,
    PARSERDISK,
    PARSERS_DICT,
    G16LogFileParserDisk,
    GJFFileParserDisk,
    SDFFileParserDisk,
    XYZFileParserDisk,
)
from molop.utils.progressbar import AdaptiveProgress


class FileParser(Enum):
    gjf = GJFFileParserDisk
    xyz = XYZFileParserDisk
    sdf = SDFFileParserDisk
    g16log = G16LogFileParserDisk

    def init(self, **kwargs: Any) -> PARSERDISK:
        return self.value(**kwargs)

    def execute(
        self,
        file_path: str,
        total_charge: int | None = None,
        total_multiplicity: int | None = None,
        **kwargs: Any,
    ) -> FILEDISK:
        return self.init(**kwargs).parse(
            file_path,
            total_charge=total_charge,
            total_multiplicity=total_multiplicity,
        )


def single_file_parser(
    file_path: str,
    possible_parsers: tuple[PARSERDISK, ...],
    release_file_content: bool = False,
) -> FILEDISK | None:
    for idx, parser in enumerate(possible_parsers):
        try:
            diskfile = parser.parse(file_path, release_file_content=release_file_content)
            return diskfile
        except Exception as e:
            if idx == len(possible_parsers) - 1:
                moloplogger.error(
                    f"Failed to parse file {file_path} with {parser.__class__.__name__}. {e}"
                )
                return None
            moloplogger.debug(
                f"Failed to parse file {file_path} with {parser.__class__.__name__}, "
                f"trying {possible_parsers[idx + 1].__name__} instead"
            )
    return None


def worker_wrapper(args_dict: dict[str, Any]) -> FILEDISK | None:
    """
    Wrapper to call single_file_parser with keyword arguments from a dictionary.
    multiprocessing.Pool.imap needs a function that takes a single argument.
    """
    return single_file_parser(**args_dict)


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
        parser_detection: Literal["auto", "gjf", "xyz", "sdf", "g16log"] = "auto",
    ) -> FileBatchModelDisk:
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
            parser_detection (Literal["auto", "gjf", "xyz", "sdf", "g16log"]):
                if "auto", use the file extension to detect the parser, else use the given parser.
        """
        valid_file_paths = []
        for file_path in file_paths:
            if isinstance(file_path, pathlib.Path):
                file_path = file_path.as_posix()
            if not os.path.isfile(file_path):
                moloplogger.warning(f"{file_path} is not a file.")
                continue
            if parser_detection == "auto" and os.path.splitext(file_path)[1] not in PARSERS_DICT:
                moloplogger.warning(f"Unsupported input file format: {file_path}")
                continue
            if file_path.endswith("molop.log"):
                continue
            valid_file_paths.append(os.path.abspath(file_path))
        if len(valid_file_paths) == 0:
            moloplogger.error("No valid input files.")
            return FileBatchModelDisk()
        if parser_detection == "auto":
            parsers = {
                key: [
                    parser(
                        forced_charge=total_charge,
                        forced_multiplicity=total_multiplicity,
                        only_extract_structure=only_extract_structure,
                        only_last_frame=only_last_frame,
                    )
                    for parser in val
                ]
                for key, val in PARSERS_DICT.items()
            }
        else:
            parsers = {
                key: [
                    FileParser.__getitem__(parser_detection).init(
                        forced_charge=total_charge,
                        forced_multiplicity=total_multiplicity,
                        only_extract_structure=only_extract_structure,
                        only_last_frame=only_last_frame,
                    )
                    for parser in val
                ]
                for key, val in PARSERS_DICT.items()
            }
        try:
            valid_file_paths.sort(key=os.path.getsize, reverse=True)
        except OSError as e:
            moloplogger.warning(f"Could not get file size for sorting, proceeding without it: {e}")
        total_tasks = [
            {
                "file_path": fp,
                "possible_parsers": parsers.get(os.path.splitext(fp)[1], ()),
                "release_file_content": release_file_content,
            }
            for fp in valid_file_paths
        ]
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
                delayed(worker_wrapper)(task)
                for task in AdaptiveProgress(
                    total_tasks,
                    total=len(total_tasks),
                    desc=f"MolOP parsing with {self.__n_jobs} processes",
                    disable=not molopconfig.show_progress_bar,
                )
            )
        else:
            results = (
                worker_wrapper(task)
                for task in AdaptiveProgress(
                    total_tasks,
                    total=len(total_tasks),
                    desc="MolOP parsing with single process",
                    disable=not molopconfig.show_progress_bar,
                )
            )
        return FileBatchModelDisk.new_batch(
            result for result in results if result is not None and len(result) > 0
        )

    @property
    def n_jobs(self) -> int:
        return self.__n_jobs

    @n_jobs.setter
    def n_jobs(self, n_jobs: int) -> None:
        self.__n_jobs = molopconfig.set_n_jobs(n_jobs)

    def __repr__(self) -> str:
        return f"FileBatchParserDisk(n_jobs={self.__n_jobs})"
