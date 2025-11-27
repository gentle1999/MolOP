"""
Author: TMJ
Date: 2024-02-04 11:04:35
LastEditors: TMJ
LastEditTime: 2025-11-23 19:59:06
Description: Command-line interface for MolOP, providing file parsing and batch processing utilities.
"""

import os
from typing import Literal, Sequence

import fire

from molop import AutoParser
from molop.config import molopconfig, moloplogger
from molop.io import FileBatchModelDisk


class MolOPCLI:
    """
    CLI for MolOP.

    This class provides a command-line interface for the MolOP library,
    allowing for file parsing, batch processing, and analysis of chemical files.

    The CLI supports a **chainable interface**, where commands can be strung together
    to perform a sequence of operations. Most commands return an instance of the CLI,
    allowing further commands to be appended. Use `-` as the separator to build a chain,
    and end with `end` to execute the chain.

    Examples:
        - Read all .log files, filter for transition states (TS), and transform them to
        .sdf format in the current directory:

        `molop read "*.log" - filter_state --state ts - transform --format sdf --output_dir "." - end`

        - Read GJF files, filter for molecules with a charge of 0, and print their file paths:

        `molop read "*.gjf" - filter_value --target charge --value 0 - paths - end`

        - Read files and generate a summary report, saving it to 'my_summary.csv':

        `molop read "*.log" - summary --output_path my_summary.csv - end`

    Returns:
        FileBatchModelDisk: The processed file batch object.
    """

    def __init__(self) -> None:
        molopconfig.enable_file_logging()
        self._file_batch: FileBatchModelDisk | None = None
        self._temp_batch: FileBatchModelDisk | None = None

    @property
    def temp_batch(self) -> FileBatchModelDisk | None:
        """
        Get the file batch object.
        """
        # If _temp_batch is None, fall back to _file_batch; this allows switching between temporary and main batch objects.
        if self._temp_batch is None:
            return self._file_batch
        else:
            return self._temp_batch

    def auto(self, file_pattern: str = "*.log", only_last_frame: bool = True):
        """
        Auto process the current directory.
        A shortcut for read(file_pattern, only_last_frame=True).

        Parameters:
            file_pattern (str): The file pattern to match. Defaults to "*.log".
            only_last_frame (bool): If True, only extract the last frame. Defaults to True.
        """
        self._file_batch = AutoParser(file_pattern, only_last_frame=only_last_frame)
        return self

    def read(
        self,
        file_path: str,
        *,
        total_charge: int | None = None,
        total_multiplicity: int | None = None,
        n_jobs: int = -1,
        only_extract_structure: bool = False,
        only_last_frame: bool = False,
        release_file_content: bool = True,
        parser_detection: Literal["auto", "gjf", "xyz", "sdf", "g16log"] = "auto",
    ):
        """
        Read the files given and set the file batch object.
        Parameters:
            file_path (str):
                use glob pattern to match files.
            total_charge (int | None):
                forced charge of the molecule, if not given, will use the charge written in the file or 0.
            total_multiplicity (int | None):
                forced multiplicity of the molecule, if not given, will use the multiplicity written in the file or 1.
            n_jobs (int):
                number of jobs to run in parallel, if -1, use all available cores.
            only_extract_structure (bool):
                if True, only extract the structure, else extract the whole file.
            only_last_frame (bool):
                if True, only extract the last frame, else extract all frames.
            release_file_content (bool):
                if True, release file content after parsing.
            parser_detection ("auto" | "gjf" | "xyz" | "sdf" | "g16log"):
                The parser to use. Defaults to "auto".
        """
        self._file_batch = AutoParser(
            file_path,
            total_charge=total_charge,
            total_multiplicity=total_multiplicity,
            n_jobs=n_jobs,
            only_extract_structure=only_extract_structure,
            only_last_frame=only_last_frame,
            release_file_content=release_file_content,
            parser_detection=parser_detection,
        )
        self._temp_batch = None  # reset temp batch
        return self

    def __checked_file_batch(self) -> FileBatchModelDisk:
        if self.temp_batch is None:
            raise ValueError(
                "No file batch found, please use `read` or `auto` command first."
            )
        return self.temp_batch

    def summary(
        self,
        output_path: str = "summary.csv",
        mode: Literal["file", "frame"] = "frame",
        frame_ids: int | Sequence[int] = -1,
    ):
        """
        Save the summary to a csv file.

        Parameters:
            output_path (str): The path to save the summary csv file. Defaults to "summary.csv".
            mode ("file" | "frame"): The mode to convert. If "file", each file will be a row. If "frame", each frame will be a row. Defaults to "frame".
            frame_ids (int | Sequence[int]): The frame IDs to convert. Defaults to -1 (last frame).
        """
        df = self.__checked_file_batch().to_summary_df(mode=mode, frameIDs=frame_ids)
        df.to_csv(output_path)
        moloplogger.info(f"Summary saved to {output_path}")
        return self

    def transform(
        self,
        format: Literal["xyz", "sdf", "cml", "gjf", "smi"],
        output_dir: str = os.getcwd(),
        frame_id: int | Literal["all"] = -1,
        embed_in_one_file: bool = True,
        **kwargs,
    ):
        """
        Transform the file format of the batch of files.

        Parameters:
            format ("xyz" | "sdf" | "cml" | "gjf" | "smi"): The target format to transform.
            output_dir (str): The directory to store the transformed files.
            frame_id (int | "all"): The frame ID to transform. If "all", all frames will be transformed. Defaults to -1.
            embed_in_one_file (bool): Whether to embed the transformed data in one file. Defaults to True.
            **kwargs: Additional keyword arguments to pass to the format_transform method.
        """
        self.__checked_file_batch().format_transform(
            format=format,
            output_dir=output_dir,
            frameID=frame_id,
            embed_in_one_file=embed_in_one_file,
            **kwargs,
        )
        return self

    def filter_state(
        self, state: Literal["ts", "error", "opt", "normal"], negate: bool = False
    ):
        """
        Filter the file batch by the state of the files.

        Parameters:
            state ("ts" | "error" | "opt" | "normal"): The state to filter.
            negate (bool): Whether to negate the filter. Defaults to False.
        """
        self._temp_batch = self.__checked_file_batch().filter_state(state, negate)
        return self

    def filter_value(
        self,
        target: Literal["charge", "multiplicity", "format"],
        value: str | float | int,
        compare: Literal["==", "!=", ">", "<", ">=", "<="] = "==",
    ):
        """
        Filter the files based on their properties.

        Parameters:
            target ("charge" | "multiplicity" | "format"): The property to filter.
            value (str | float | int): The value to filter.
            compare ("==" | "!=" | ">" | "<" | ">=" | "<="): The comparison operator. Defaults to "==".
        """
        self._temp_batch = self.__checked_file_batch().filter_value(
            target=target, value=value, compare=compare
        )
        return self

    def paths(self):
        """
        Print the file paths of each file.
        """
        for _file in self.__checked_file_batch():
            print(_file.file_path)
        return self

    def end(self):
        """
        End the command chain, stop printing help comments.
        """
        self._file_batch = None
        self._temp_batch = None

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
