"""
Author: TMJ
Date: 2025-08-20 22:55:18
LastEditors: TMJ
LastEditTime: 2025-08-21 00:15:32
Description: 请填写简介
"""

import os
import time
from typing import Dict, List, Optional, Tuple, Type

from joblib import Parallel, cpu_count, delayed
from tqdm import tqdm

from molop.config import molopconfig
from molop.io.coords_models import GJFFileDisk, SDFFileDisk, XYZFileDisk
from molop.io.coords_parsers import (
    GJFFileParserDisk,
    SDFFileParserDisk,
    XYZFileParserDisk,
)
from molop.io.FileBatchModelDisk import FileBatchModelDisk
from molop.io.QM_models import G16LogFileDisk
from molop.io.QM_parsers import G16LogFileParserDisk
from molop.logger.logger import moloplogger

FILEDISK = G16LogFileDisk | GJFFileDisk | XYZFileDisk | SDFFileDisk
PARSER = (
    GJFFileParserDisk | XYZFileParserDisk | SDFFileParserDisk | G16LogFileParserDisk
)
PARSERTYPE = Type[PARSER]
PARSERS_DICT: Dict[str, Tuple[PARSERTYPE, ...]] = {
    ".gjf": (GJFFileParserDisk,),
    ".gau": (G16LogFileParserDisk,),
    ".com": (GJFFileParserDisk,),
    ".gjc": (GJFFileParserDisk,),
    ".log": (G16LogFileParserDisk,),
    ".g16": (G16LogFileParserDisk,),
    ".gal": (G16LogFileParserDisk,),
    ".xyz": (XYZFileParserDisk,),
    ".sdf": (SDFFileParserDisk,),
    ".mol": (SDFFileParserDisk,),
    ".out": (G16LogFileParserDisk,),
    ".irc": (G16LogFileParserDisk,),
}


def single_file_parser(
    file_path: str,
    possible_parsers: Tuple[PARSER, ...],
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
                f"trying {possible_parsers[idx+1].__name__} instead"
            )


class FileBatchParserDisk:
    __n_jobs: int

    def __init__(self, n_jobs: int = -1):
        self.__n_jobs = self._set_n_jobs(n_jobs)

    def parse(
        self,
        file_paths: List[str],
        total_charge: Optional[int] = None,
        total_multiplicity: Optional[int] = None,
        only_extract_structure=False,
        only_last_frame=False,
    ) -> FileBatchModelDisk:
        time1 = time.time()
        valid_file_paths = []
        for file_path in file_paths:
            if not os.path.isfile(file_path):
                moloplogger.warning(f"{file_path} is not a file.")
                continue
            if os.path.splitext(file_path)[1] not in PARSERS_DICT:
                moloplogger.warning(f"Unsupported input file format: {file_path}")
                continue
            if file_path.endswith("molop.log"):
                continue
            valid_file_paths.append(os.path.abspath(file_path))
        if len(valid_file_paths) == 0:
            moloplogger.error("No valid input files.")
            return FileBatchModelDisk()
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
        task_list_small_files = [
            {
                "file_path": file_path,
                "possible_parsers": parsers[os.path.splitext(file_path)[1]],
            }
            for file_path in valid_file_paths
            if os.path.getsize(file_path) < molopconfig.parallel_max_size
        ]
        task_list_large_files = [
            {
                "file_path": file_path,
                "possible_parsers": parsers[os.path.splitext(file_path)[1]],
            }
            for file_path in valid_file_paths
            if os.path.getsize(file_path) >= molopconfig.parallel_max_size
        ]
        total_tasks = task_list_small_files + task_list_large_files
        time2 = time.time()
        moloplogger.warning(f"Time to create parsers: {time2 - time1:.2f}s")
        diskfiles: List[FILEDISK] = []
        if self.__n_jobs > 1 and len(valid_file_paths) > self.__n_jobs:
            if only_last_frame:
                if molopconfig.show_progress_bar:
                    diskfiles.extend(
                        diskfile
                        for diskfile in tqdm(
                            Parallel(
                                return_as="generator",
                                n_jobs=self.__n_jobs,
                                pre_dispatch="1.5*n_jobs",
                            )(
                                delayed(single_file_parser)(**arguments)
                                for arguments in total_tasks
                            ),
                            desc=f"MolOP parsing with {self.__n_jobs} jobs",
                            total=len(total_tasks),
                        )
                        if diskfile is not None and len(diskfile) > 0
                    )
                else:
                    diskfiles.extend(
                        diskfile
                        for diskfile in Parallel(
                            return_as="generator",
                            n_jobs=self.__n_jobs,
                            pre_dispatch="1.5*n_jobs",
                        )(
                            delayed(single_file_parser)(**arguments)
                            for arguments in total_tasks
                        )
                        if diskfile is not None and len(diskfile) > 0
                    )
            else:
                if len(task_list_small_files) > 0:
                    if molopconfig.show_progress_bar:
                        diskfiles.extend(
                            diskfile
                            for diskfile in tqdm(
                                Parallel(
                                    return_as="generator",
                                    n_jobs=self.__n_jobs,
                                    pre_dispatch="1.5*n_jobs",
                                )(
                                    delayed(single_file_parser)(**arguments)
                                    for arguments in task_list_small_files
                                ),
                                desc=f"MolOP parsing with {self.__n_jobs} jobs",
                                total=len(task_list_small_files),
                            )
                            if diskfile is not None and len(diskfile) > 0
                        )
                    else:
                        diskfiles.extend(
                            diskfile
                            for diskfile in Parallel(
                                return_as="generator",
                                n_jobs=self.__n_jobs,
                                pre_dispatch="1.5*n_jobs",
                            )(
                                delayed(single_file_parser)(**arguments)
                                for arguments in task_list_small_files
                            )
                            if diskfile is not None and len(diskfile) > 0
                        )
                if len(task_list_large_files) > 0:
                    if molopconfig.show_progress_bar:
                        diskfiles.extend(
                            diskfile
                            for diskfile in tqdm(
                                (
                                    single_file_parser(**arguments)
                                    for arguments in (task_list_large_files)
                                ),
                                desc="MolOP parsing with single process",
                                total=len(task_list_large_files),
                            )
                            if diskfile is not None and len(diskfile) > 0
                        )
                    else:
                        diskfiles.extend(
                            diskfile
                            for diskfile in (
                                single_file_parser(**arguments)
                                for arguments in (task_list_large_files)
                            )
                            if diskfile is not None and len(diskfile) > 0
                        )
        else:
            if molopconfig.show_progress_bar:
                diskfiles.extend(
                    diskfile
                    for diskfile in tqdm(
                        (
                            single_file_parser(**arguments)
                            for arguments in total_tasks
                        ),
                        total=len(total_tasks),
                        desc="MolOP parsing with single process",
                    )
                    if diskfile is not None and len(diskfile) > 0
                )
            else:
                diskfiles.extend(
                    diskfile
                    for diskfile in (
                        single_file_parser(**arguments)
                        for arguments in total_tasks
                    )
                    if diskfile is not None and len(diskfile) > 0
                )
        time3 = time.time()
        moloplogger.warning(f"Time to parse: {time3 - time2:.2f}s")
        return FileBatchModelDisk.new_batch(diskfiles)

    @property
    def n_jobs(self):
        return self.__n_jobs

    @n_jobs.setter
    def n_jobs(self, n_jobs: int):
        self.__n_jobs = self._set_n_jobs(n_jobs)

    def _set_n_jobs(self, n_jobs: int):
        return (
            min(n_jobs, cpu_count(), molopconfig.max_jobs)
            if n_jobs > 0
            else min(cpu_count(), molopconfig.max_jobs)
        )
