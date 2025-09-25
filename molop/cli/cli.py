"""
Author: TMJ
Date: 2024-02-04 11:04:35
LastEditors: TMJ
LastEditTime: 2024-02-12 16:49:38
Description: Command-line interface for MolOP, providing file parsing and batch processing utilities.
"""

import os
from typing import Literal

import fire

from molop import AutoParser
from molop.config import molopconfig
from molop.io import FileBatchModelDisk


class MolOPCLI:
    """
    CLI for MolOP.
    """

    def __init__(self) -> None:
        molopconfig.enable_file_logging()
        self._file_batch: FileBatchModelDisk | None = None
        self._temp_batch: FileBatchModelDisk | None = None

    @property
    def temp_batch(self):
        """
        Get the file batch object.
        """
        # If _temp_batch is None, fall back to _file_batch; this allows switching between temporary and main batch objects.
        if self._temp_batch is None:
            return self._file_batch
        else:
            return self._temp_batch

    def auto(self, file_pattern: str = "*.log", hong=True):
        """
        Auto process the current directory.

        Parameters:
        self._file_batch = AutoParser(r"*.log", only_last_frame=True)
        """
        self._file_batch = AutoParser(file_pattern, only_last_frame=True)
        # self.temp_batch.to_summary_csv(use_hong_style=hong)

    def read(
        self,
        file_path: str,
        *,
        charge: int = 0,
        multiplicity: int = 1,
        n_jobs: int = -1,
        only_extract_structure: bool = False,
        only_last_frame: bool = False,
    ):
        """
        Read the files given and set the file batch object.
        Parameters:
            file_path (str):
                use regax to match files.
            charge (int):
                forced charge of the molecule, if not given, will use the charge written in the file or 0.
            multiplicity (int):
                forced multiplicity of the molecule, if not given, will use the multiplicity written in the file or 1.
            n_jobs (int):
                number of jobs to run in parallel, if -1, use all available cores.
            only_extract_structure (bool):
                if True, only extract the structure, else extract the whole file.
            only_last_frame (bool):
                if True, only extract the last frame, else extract all frames.
        """
        self._file_batch = AutoParser(
            file_path,
            total_charge=charge,
            total_multiplicity=multiplicity,
            n_jobs=n_jobs,
            only_extract_structure=only_extract_structure,
            only_last_frame=only_last_frame,
        )
        return self

    @property
    def __checked_file_batch(self):
        if self.temp_batch is None:
            raise ValueError("No file batch found, please use `read` command first.")
        return self.temp_batch

    # def summary(
    #     self,
    #     file_dir: str = None,
    #     full: bool = True,
    #     with_units: bool = True,
    #     use_hong_style=False,
    # ):
    #     """
    #     Save the summary csv file.

    #     Parameters:
    #         file_dir str:
    #             the directory to save the summary csv file.
    #         full bool:
    #             if True, save all the information, else save only the essential information.
    #         with_units bool:
    #             if True, save the information with units, else save the information without units.
    #     """
    #     self.temp_batch.to_summary_csv(
    #         file_dir, full, with_units, use_hong_style=use_hong_style
    #     )
    #     return self

    def transform(
        self,
        format: Literal["xyz", "sdf", "cml", "gjf", "smi"],
        output_dir: str = os.getcwd(),
        frame_id: int | Literal["all"] = -1,
        embed_in_one_file: bool = True,
    ):
        """
        Transform the file format of the batch of files.

        Parameters:
            format ("xyz" | "sdf" | "cml" | "gjf" | "smi"): The target format to transform.
            output_dir (str): The directory to store the transformed files.
            frame_id (int | "all"): The frame ID to transform. If "all", all frames will be transformed.
            embed_in_one_file (bool): Whether to embed the transformed data in one file.
        """
        self.__checked_file_batch.format_transform(
            format=format,
            output_dir=output_dir,
            frameID=frame_id,
            embed_in_one_file=embed_in_one_file,
        )
        return self

    def filter_state(
        self, state: Literal["ts", "error", "opt", "normal"], negate: bool = False
    ):
        """
        Filter the file batch by the state of the files.

        Parameters:
            state ("ts" | "error" | "opt" | "normal"): The state to filter.
            negate (bool): Whether to negate the filter.
        """
        self._temp_batch = self.__checked_file_batch.filter_state(state, negate)
        return self

    def paths(self):
        """
        Print the file paths of each file.
        """
        for _file in self.__checked_file_batch:
            print(_file.file_path)
        return self

    def end(self):
        """
        End the command chain, stop printing help comments.
        """
        self._file_batch = None

    def quiet(self):
        """
        Quiet the command chain, stop printing logs.
        """
        molopconfig.quiet()
        return self

    def verbose(self):
        """
        Verbose the command chain, print logs.
        """
        molopconfig.verbose()
        return self

    def log_off(self):
        """
        Turn off the log.
        """
        molopconfig.disable_file_logging()
        return self

    def log_on(self):
        """
        Turn on the log.
        """
        molopconfig.enable_file_logging()
        return self


def app():
    fire.Fire(MolOPCLI)
