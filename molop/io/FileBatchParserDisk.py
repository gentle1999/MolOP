"""
Author: TMJ
Date: 2025-08-20 22:55:18
LastEditors: TMJ
LastEditTime: 2025-08-21 00:15:32
Description: 请填写简介
"""

import os
from enum import Enum
from typing import Any, Dict, Iterable, Literal, Optional, Tuple

from joblib import Parallel, delayed
from tqdm import tqdm

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


class FileParser(Enum):
    gjf = GJFFileParserDisk
    xyz = XYZFileParserDisk
    sdf = SDFFileParserDisk
    g16log = G16LogFileParserDisk

    def init(self, **kwargs):
        return self.value(**kwargs)

    def execute(
        self,
        file_path: str,
        total_charge: Optional[int] = None,
        total_multiplicity: Optional[int] = None,
        **kwargs,
    ):
        return self.init(**kwargs).parse(
            file_path,
            total_charge=total_charge,
            total_multiplicity=total_multiplicity,
        )


def single_file_parser(
    file_path: str,
    possible_parsers: Tuple[PARSERDISK, ...],
) -> Optional[FILEDISK]:
    for idx, parser in enumerate(possible_parsers):
        try:
            diskfile = parser.parse(file_path)
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


def worker_wrapper(args_dict: Dict[str, Any]) -> Optional[FILEDISK]:
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
        file_paths: Iterable[str],
        total_charge: Optional[int] = None,
        total_multiplicity: Optional[int] = None,
        only_extract_structure=False,
        only_last_frame=False,
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
            if not os.path.isfile(file_path):
                moloplogger.warning(f"{file_path} is not a file.")
                continue
            if parser_detection == "auto":
                if os.path.splitext(file_path)[1] not in PARSERS_DICT:
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
            moloplogger.warning(
                f"Could not get file size for sorting, proceeding without it: {e}"
            )
        total_tasks = [
            {"file_path": fp, "possible_parsers": parsers.get(os.path.splitext(fp)[1], ())}
            for fp in valid_file_paths
        ]
        # Determine if parallel processing should be used
        use_parallel = self.__n_jobs > 1 and len(total_tasks) > self.__n_jobs

        if use_parallel:
            moloplogger.info(
                f"Using {self.__n_jobs} processes for {len(total_tasks)} tasks."
            )
            results = Parallel(
                n_jobs=self.__n_jobs,
                maxtasks_per_child=50,
                return_as="generator",
                max_nbytes=molopconfig.parallel_max_size,
            )(
                delayed(worker_wrapper)(task)
                for task in tqdm(
                    total_tasks,
                    desc=f"MolOP parsing with {self.__n_jobs} processes",
                    disable=not molopconfig.show_progress_bar,
                )
            )
        else:
            desc = "MolOP parsing with single process"
            results = (
                worker_wrapper(task)
                for task in tqdm(
                    total_tasks,
                    desc=desc,
                    disable=not molopconfig.show_progress_bar,
                )
            )
        return FileBatchModelDisk.new_batch(
            result for result in results if result is not None and len(result) > 0
        )

    @property
    def n_jobs(self):
        return self.__n_jobs

    @n_jobs.setter
    def n_jobs(self, n_jobs: int):
        self.__n_jobs = molopconfig.set_n_jobs(n_jobs)

    def __repr__(self) -> str:
        return f"FileBatchParserDisk(n_jobs={self.__n_jobs})"
