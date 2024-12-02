"""
Author: TMJ
Date: 2024-01-25 22:46:35
LastEditors: TMJ
LastEditTime: 2024-01-25 23:10:59
Description: 请填写简介
"""

import os
from collections import OrderedDict
from collections.abc import MutableMapping
from typing import Dict, List, Union, overload

import pandas as pd
from joblib import Parallel, cpu_count, delayed
from tqdm import tqdm

from molop.config import molopconfig
from molop.io.bases.BaseMolFileParser import BaseMolFileParser
from molop.io.bases.BaseMolFrameParser import MolFrameType
from molop.io.coords_file.GJFFileParser import GJFFileParser
from molop.io.coords_file.SDFFileParser import SDFFileParser
from molop.io.coords_file.XYZFileParser import XYZFileParser
from molop.io.qm_file.G16FchkFileParser import G16FchkFileParser
from molop.io.qm_file.G16LogFileParser import G16LogFileParser
from molop.io.qm_file.XTBFileParser import XTBFileParser
from molop.io.types import PARSERTYPES
from molop.logger.logger import moloplogger

parsers: Dict[str, PARSERTYPES] = {
    ".gjf": (GJFFileParser,),
    ".gau": (GJFFileParser,),
    ".com": (GJFFileParser,),
    ".gjc": (GJFFileParser,),
    ".log": (G16LogFileParser,),
    ".g16": (G16LogFileParser,),
    ".gal": (G16LogFileParser,),
    ".xyz": (XYZFileParser,),
    ".sdf": (SDFFileParser,),
    ".mol": (SDFFileParser,),
    ".out": (G16LogFileParser, XTBFileParser),
    ".irc": (G16LogFileParser,),
    ".fchk": (G16FchkFileParser,),
    ".fck": (G16FchkFileParser,),
    ".fch": (G16FchkFileParser,),
}

qm_parsers = (G16LogFileParser, G16FchkFileParser, XTBFileParser)


def singlefile_parser(
    file_path: str,
    charge=0,
    multiplicity=1,
    only_extract_structure=False,
    only_last_frame=False,
) -> PARSERTYPES:
    """
    A function to parse a single file based on its format.

    Core Logic:
        - Checks if the provided `file_path` is an actual file.
        - Determines the file format by getting the file extension.
        - If the format is not supported, logs an error and returns None.
        - Iterates through the available parsers for that format.
        - Tries to instantiate and return the parser object with the given parameters.
        - If an exception occurs during parsing, it logs the error and tries the next parser in the list.
        - If all parsers fail, logs the final error and returns None.
        - If the `file_path` points to a directory, logs an error and returns None.

    Parameters:
        file_path (str):
            The path of the file.
        charge (int):
            The molecular charge. Defaults to None.
        multiplicity (int):
            The molecular spin multiplicity. Defaults to None.
        only_extract_structure (bool):
            Whether to extract structure only. Defaults to False.
        only_last_frame (bool):
            Whether to extract the last frame only. Defaults to False.

    Returns:
        PARSERTYPES: An instance of the appropriate parser class.
    """
    if os.path.isfile(file_path):
        _, file_format = os.path.splitext(file_path)
        if file_format not in parsers:
            moloplogger.error(
                "Unknown file format: {}, {}".format(file_format, file_path)
            )
            return None
        for idx, parser in enumerate(parsers[file_format]):
            try:
                # Instantiate and return specific parser classes with their respective arguments
                if parser in qm_parsers:
                    return parser(
                        file_path=file_path,
                        charge=charge,
                        multiplicity=multiplicity,
                        only_extract_structure=only_extract_structure,
                        only_last_frame=only_last_frame,
                    )
                elif parser in (GJFFileParser,):
                    return parser(
                        file_path=file_path, charge=charge, multiplicity=multiplicity
                    )
                elif parser in (XYZFileParser,):
                    return parser(
                        file_path=file_path,
                        charge=charge,
                        multiplicity=multiplicity,
                        only_last_frame=only_last_frame,
                    )
                elif parser in (SDFFileParser,):
                    return parser(file_path=file_path, only_last_frame=only_last_frame)
                else:
                    return parser(file_path=file_path)
            except Exception as e:
                if idx == len(parsers[file_format]) - 1:
                    moloplogger.error(
                        f"Failed to parse file {file_path} with {parser.__name__}. {e}"
                    )
                    return None
                moloplogger.debug(
                    f"Failed to parse file {file_path} with {parser.__name__}, trying {parsers[file_format][idx+1].__name__} instead"
                )
    elif os.path.isdir(file_path):
        moloplogger.error(f"{file_path} is not a file.")
        return None


class FileParserBatch(MutableMapping):
    __n_jobs: int
    __parsers: Dict[str, PARSERTYPES]

    def __init__(
        self,
        n_jobs: int = -1,
    ) -> None:
        self.__n_jobs = self._set_n_jobs(n_jobs)
        self.__parsers: Dict[str, PARSERTYPES] = OrderedDict()

    def __contains__(self, key: Union[str, PARSERTYPES]) -> bool:
        if isinstance(key, str):
            return key in self.__parsers.keys()
        else:
            return key in self.__parsers.values()

    @overload
    def __getitem__(self, key: int) -> PARSERTYPES: ...

    @overload
    def __getitem__(self, key: slice) -> "FileParserBatch": ...

    @overload
    def __getitem__(self, key: str) -> PARSERTYPES: ...

    def __getitem__(self, key: Union[int, str, slice]) -> PARSERTYPES:
        if isinstance(key, int):
            return list(self.__parsers.values())[key]
        if isinstance(key, slice):
            slicedkeys = list(self.__parsers.keys())[key]
            return self.__new_batch([self.__parsers[k] for k in slicedkeys])
        else:
            return self.__parsers[key]

    def __setitem__(self, key: str, value: PARSERTYPES):
        if not isinstance(key, str):
            raise TypeError(f"key should be a string, got {type(key)}")
        self.__parsers[key] = value

    def __delitem__(self, key: str):
        if not isinstance(key, str):
            raise TypeError(f"key should be a string, got {type(key)}")
        del self.__parsers[key]

    def __len__(self) -> int:
        return len(self.__parsers)

    def __iter__(self):
        self.__index = 0
        return self

    def __next__(
        self,
    ) -> PARSERTYPES:
        if self.__index >= len(self):
            raise StopIteration
        else:
            self.__index += 1
            return self[self.__index - 1]

    def add_files(
        self,
        files: List[str],
        charge=0,
        multiplicity=1,
        only_extract_structure=False,
        only_last_frame=False,
    ):
        file_paths = []
        for file_path in files:
            if not os.path.isfile(file_path):
                moloplogger.warning(f"{file_path} is not a file.")
                continue
            if os.path.splitext(file_path)[1] not in parsers:
                moloplogger.warning(f"Unsupported input file format: {file_path}")
                continue
            if file_path.endswith("molop.log"):
                continue
            file_paths.append(os.path.abspath(file_path))

        arguments_list_small_files = [
            {
                "file_path": file_path,
                "charge": charge,
                "multiplicity": multiplicity,
                "only_extract_structure": only_extract_structure,
                "only_last_frame": only_last_frame,
            }
            for file_path in file_paths
            if os.path.getsize(file_path) < molopconfig.parallel_max_size
        ]
        arguments_list_large_files = [
            {
                "file_path": file_path,
                "charge": charge,
                "multiplicity": multiplicity,
                "only_extract_structure": only_extract_structure,
                "only_last_frame": only_last_frame,
            }
            for file_path in file_paths
            if os.path.getsize(file_path) >= molopconfig.parallel_max_size
        ]
        if self.__n_jobs > 1 and len(file_paths) > self.__n_jobs:
            if molopconfig.show_progress_bar:
                self.__parsers.update(
                    {
                        **{
                            parser.file_path: parser
                            for parser in tqdm(
                                Parallel(
                                    return_as="generator",
                                    n_jobs=self.__n_jobs,
                                    pre_dispatch="1.5*n_jobs",
                                )(
                                    delayed(singlefile_parser)(**arguments)
                                    for arguments in arguments_list_small_files
                                ),
                                total=len(arguments_list_small_files),
                                desc=f"MolOP parsing with {self.__n_jobs} jobs",
                            )
                            if parser is not None and len(parser) > 0
                        },
                        **{
                            parser.file_path: parser
                            for parser in tqdm(
                                (
                                    singlefile_parser(**arguments)
                                    for arguments in (arguments_list_large_files)
                                ),
                                total=len(arguments_list_large_files),
                                desc=f"MolOP parsing with single thread for large files",
                            )
                            if parser is not None and len(parser) > 0
                        },
                    }
                )
            else:
                self.__parsers.update(
                    {
                        **{
                            parser.file_path: parser
                            for parser in Parallel(
                                return_as="generator",
                                n_jobs=self.__n_jobs,
                                pre_dispatch="1.5*n_jobs",
                            )(
                                delayed(singlefile_parser)(**arguments)
                                for arguments in arguments_list_small_files
                            )
                            if parser is not None and len(parser) > 0
                        },
                        **{
                            parser.file_path: parser
                            for parser in (
                                singlefile_parser(**arguments)
                                for arguments in (arguments_list_large_files)
                            )
                            if parser is not None and len(parser) > 0
                        },
                    }
                )
        else:
            if molopconfig.show_progress_bar:
                self.__parsers.update(
                    {
                        parser.file_path: parser
                        for parser in tqdm(
                            (
                                singlefile_parser(**arguments)
                                for arguments in (
                                    arguments_list_small_files
                                    + arguments_list_large_files
                                )
                            ),
                            total=len(
                                arguments_list_small_files + arguments_list_large_files
                            ),
                            desc=f"MolOP parsing with single thread",
                        )
                        if parser is not None and len(parser) > 0
                    }
                )
            else:
                self.__parsers.update(
                    {
                        parser.file_path: parser
                        for parser in (
                            singlefile_parser(**arguments)
                            for arguments in (
                                arguments_list_small_files + arguments_list_large_files
                            )
                        )
                        if parser is not None and len(parser) > 0
                    }
                )
        self.__parsers = OrderedDict(sorted(self.__parsers.items()))

        moloplogger.info(
            f"{len(file_paths) - len(self.__parsers)} files failed to parse, {len(self.__parsers)} successfully parsed"
        )

    def add_file_parsers(self, file_parsers: List[PARSERTYPES]) -> None:
        for parser in file_parsers:
            if not isinstance(parser, BaseMolFileParser):
                raise TypeError(
                    f"file_parsers should be a list of BaseFileParser, got {type(parser)}"
                )
            if parser.file_path not in self.__parsers:
                self[parser.file_path] = parser
            else:
                moloplogger.warning(
                    f"File {parser.file_path} already exists in the batch, skipped"
                )

    @property
    def file_paths(self) -> List[str]:
        """return a list of file paths"""
        return [parser.file_path for parser in self]

    @property
    def file_names(self) -> List[str]:
        """return a list of file names"""
        return [parser.file_name for parser in self]

    def to_GJF_file(
        self,
        file_path: str = None,
        charge: int = None,
        multiplicity: int = None,
        prefix: str = "",
        suffix="",
        template: str = None,
        chk: bool = True,
        oldchk: bool = False,
        frameID=-1,
    ) -> None:
        """
        Write the GJF file.

        Parameters:
            file_path (str):
                The path to write the GJF file. If not specified, will be generated in situ.
            charge (int):
                The forced charge. If specified, will be used to overwrite the charge in the gjf file.
            multiplicity (int):
                The forced multiplicity. If specified, will be used to overwrite the multiplicity in the gjf file.
            template (str):
                path to read a gjf file as a template.
            prefix (str):
                prefix to add to the beginning of the gjf file, priority is higher than template.
            suffix (str):
                suffix to add to the end of the gjf file, priority is higher than template.
            chk (bool):
                If true, add the chk keyword to the link0 section. Will use the file name as the chk file name.
            oldchk (bool):
                If true, add the oldchk keyword to the link0 section. Will use the file name as the chk file name.
            frameID (int):
                The frame ID to write.
        Returns:
            str: The path to the GJF file.
        """
        if file_path is None:
            file_path = os.path.curdir
        if not os.path.isdir(file_path):
            raise NotADirectoryError(f"{file_path} is not a directory")
        for parser in self:
            parser.to_GJF_file(
                os.path.join(file_path, parser.pure_filename + ".gjf"),
                charge=charge,
                multiplicity=multiplicity,
                prefix=prefix,
                suffix=suffix,
                template=template,
                chk=chk,
                oldchk=oldchk,
                frameID=frameID,
            )
        moloplogger.info(f"gjf files saved to {os.path.abspath(file_path)}")

    def to_XYZ_file(self, file_path: str = None) -> None:
        if file_path is None:
            file_path = os.path.curdir
        if not os.path.isdir(file_path):
            raise NotADirectoryError(f"{file_path} is not a directory")
        for parser in self:
            parser.to_XYZ_file(os.path.join(file_path, parser.pure_filename + ".xyz"))
        moloplogger.info(f"xyz files saved to {os.path.abspath(file_path)}")

    def to_SDF_file(self, file_path: str = None) -> None:
        if file_path is None:
            file_path = os.path.curdir
        if not os.path.isdir(file_path):
            raise NotADirectoryError(f"{file_path} is not a directory")
        for parser in self:
            parser.to_SDF_file(os.path.join(file_path, parser.pure_filename + ".sdf"))
        moloplogger.info(f"sdf files saved to {os.path.abspath(file_path)}")

    def to_chemdraw(self, file_path: str = None, frameID=-1, keep3D=True) -> None:
        if file_path is None:
            file_path = os.path.curdir
        if not os.path.isdir(file_path):
            raise NotADirectoryError(f"{file_path} is not a directory")
        for parser in self:
            parser.to_chemdraw(
                os.path.join(file_path, parser.pure_filename + ".cdxml"),
                frameID=frameID,
                keep3D=keep3D,
            )
        moloplogger.info(f"chemdraw files saved to {os.path.abspath(file_path)}")

    def replace_substituent(
        self,
        query_smi: str,
        replacement_smi: str,
        bind_idx: int = None,
        replace_all: bool = False,
        attempt_num: int = 10,
        crowding_threshold: float = 0.75,
        angle_split: int = 10,
        randomSeed: int = 114514,
        *,
        replacement_relative_idx: int = 0,
        replacement_absolute_idx: Union[int, None] = None,
        frameID=-1,
    ) -> "FileParserBatch":
        """
        Replace the substituent with the given SMARTS. The substituent is defined by the query_smi, and the new substituent is defined by the replacement_smi.

        Parameters:
            query_smi (str):
                The SMARTS to query the substituent in the original molecule.
            replacement_smi (str):
                The SMARTS of new substituent. The bind atom is the first atom of the replacement_smi.
            bind_idx (int):
                The index of the atom to bind the new substituent. The default is None, which means to replace the first legal atom in original molecule.
                If specified, try to replace the atom. User should meke sure the atom is legal.
                Detail example in (Repalce Substituent)[Repalce Substituent]
            replace_all (bool):
                If True, replace all the substituent queried in the original molecule.
            attempt_num (int):
                Max attempt times to replace the substituent. Each time a new substituent conformation will be used for substitution.
            crowding_threshold (float):
                The threshold of crowding. If the new substituent is too crowded, the substitution will be rejected.
            angle_split (int):
                The number of angles to rotate the new substituent to find a legal conformation.
            randomSeed (int):
                The random seed for the random number generator.
                key between `start_idx` and `end_idx` and replace the base group where `end_idx` is located
            replacement_relative_idx (int):
                The relative index of the radical atom in the replacement molecule to be
                transformed to the first atom.
            replacement_absolute_idx (Union[int, None]):
                Priority is higher than replacement_relative_idx.
                The absolute index of the radical atom in the replacement molecule to be
                transformed to the first atom.
                If None, the function will try to find the first atom in the replacement
                molecule that is a radical atom.
            frameID (int):
                The frameID to replace.
        Returns:
            FileParserBatch: The new `FileParserBatch`.
        """
        new_parsers = []
        if molopconfig.show_progress_bar:
            for parser in tqdm(self):
                try:
                    temp_parser = parser[frameID].replace_substituent(
                        query=query_smi,
                        replacement=replacement_smi,
                        bind_idx=bind_idx,
                        replace_all=replace_all,
                        attempt_num=attempt_num,
                        crowding_threshold=crowding_threshold,
                        angle_split=angle_split,
                        randomSeed=randomSeed,
                        replacement_relative_idx=replacement_relative_idx,
                        replacement_absolute_idx=replacement_absolute_idx,
                    )
                    temp_file_path = temp_parser.to_SDF_file()
                    new_parsers.append(
                        SDFFileParser(file_path=temp_file_path, only_last_frame=True)
                    )
                    os.remove(temp_file_path)
                except:
                    moloplogger.error(
                        f"Failed to replace substituent from {query_smi} to {replacement_smi} in {parser.file_path}, {parser.file_name}"
                    )
        else:
            for parser in self:
                try:
                    temp_parser = parser[frameID].replace_substituent(
                        query=query_smi,
                        replacement=replacement_smi,
                        bind_idx=bind_idx,
                        replace_all=replace_all,
                        attempt_num=attempt_num,
                        crowding_threshold=crowding_threshold,
                        angle_split=angle_split,
                        randomSeed=randomSeed,
                        replacement_relative_idx=replacement_relative_idx,
                        replacement_absolute_idx=replacement_absolute_idx,
                    )
                    temp_file_path = temp_parser.to_SDF_file()
                    new_parsers.append(
                        SDFFileParser(file_path=temp_file_path, only_last_frame=True)
                    )
                    os.remove(temp_file_path)
                except:
                    moloplogger.error(
                        f"Failed to replace substituent from {query_smi} to {replacement_smi} in {parser.file_path}, {parser.file_name}"
                    )

        new_batch = FileParserBatch()
        new_batch.add_file_parsers(new_parsers)
        moloplogger.info(
            f"{len(new_parsers)} files successfully replaced, {len(self.__parsers) - len(new_parsers)} files failed to replace"
        )
        return new_batch

    def reset_atom_index(
        self, mapping_smarts: str, frameID: int = -1
    ) -> "FileParserBatch":
        """only consider the last frame of each file"""
        new_parsers = []
        if molopconfig.show_progress_bar:
            for parser in tqdm(self):
                try:
                    temp_parser = parser[frameID].reset_atom_index(
                        mapping_smarts=mapping_smarts
                    )
                    temp_file_path = temp_parser.to_SDF_file()
                    new_parsers.append(
                        SDFFileParser(file_path=temp_file_path, only_last_frame=True)
                    )
                    os.remove(temp_file_path)
                except Exception as e:
                    moloplogger.error(
                        f"{e}: Failed to reset atom index by {mapping_smarts} in {parser.file_path}, {parser.file_name}"
                    )
        else:
            for parser in self:
                try:
                    temp_parser = parser[frameID].reset_atom_index(
                        mapping_smarts=mapping_smarts
                    )
                    temp_file_path = temp_parser.to_SDF_file()
                    new_parsers.append(
                        SDFFileParser(file_path=temp_file_path, only_last_frame=True)
                    )
                    os.remove(temp_file_path)
                except Exception as e:
                    moloplogger.error(
                        f"{e}: Failed to reset atom index by {mapping_smarts} in {parser.file_path}, {parser.file_name}"
                    )
        new_batch = FileParserBatch()
        new_batch.add_file_parsers(new_parsers)
        moloplogger.info(
            f"{len(new_parsers)} files successfully replaced, {len(self.__parsers) - len(new_parsers)} files failed to reset_index"
        )
        return new_batch

    def standard_orient(
        self,
        anchor_list: List[int],
        frameID: int = -1,
    ) -> "FileParserBatch":
        """
        Depending on the input `idx_list`, `translate_anchor`, `rotate_anchor_to_axis`, and `rotate_anchor_to_plane` are executed in order to obtain the normalized oriented molecule.

        Sub-functions:
            - `translate_anchor`: Translate the entire molecule so that the specified atom reaches the origin.
            - `rotate_anchor_to_axis`: Rotate the specified second atom along the axis passing through the origin so that it reaches the positive half-axis of the X-axis.
            - `rotate_anchor_to_plane`: Rotate along the axis passing through the origin so that the specified third atom reaches quadrant 1 or 2 of the XY plane.

        Parameters:
            anchor_list (List[int]):
                A list of indices of the atoms to be translated to origin, rotated to X axis, and rotated again to XY face:

                - If length is 1, execute `translate_anchor`
                - If length is 2, execute `translate_anchor` and `rotate_anchor_to_axis`
                - If length is 3, execute `translate_anchor`, `rotate_anchor_to_axis` and `rotate_anchor_to_plane`
                - If the length of the input `idx_list` is greater than 3, subsequent atomic numbers are not considered.


        Returns:
            FileParserBatch: The new `FileParserBatch`.
        """
        new_parsers = []
        if molopconfig.show_progress_bar:
            for parser in tqdm(self):
                try:
                    temp_parser = parser[frameID].standard_orient(anchor_list)
                    temp_file_path = temp_parser.to_SDF_file()
                    new_parsers.append(
                        SDFFileParser(file_path=temp_file_path, only_last_frame=True)
                    )
                    os.remove(temp_file_path)
                except Exception as e:
                    moloplogger.error(
                        f"{e}: Failed to standard_orient by {anchor_list} in {parser.file_path}, {parser.file_name}"
                    )
        else:
            for parser in self:
                try:
                    temp_parser = parser[frameID].standard_orient(anchor_list)
                    temp_file_path = temp_parser.to_SDF_file()
                    new_parsers.append(
                        SDFFileParser(file_path=temp_file_path, only_last_frame=True)
                    )
                    os.remove(temp_file_path)
                except Exception as e:
                    moloplogger.error(
                        f"{e}: Failed to standard_orient by {anchor_list} in {parser.file_path}, {parser.file_name}"
                    )
        new_batch = FileParserBatch()
        new_batch.add_file_parsers(new_parsers)
        moloplogger.info(
            f"{len(new_parsers)} files successfully replaced, {len(self.__parsers) - len(new_parsers)} files failed to standard_orient"
        )
        return new_batch

    @classmethod
    def __new_batch(cls, parsers: List[PARSERTYPES]):
        new_batch = cls()
        new_batch.add_file_parsers(parsers)
        return new_batch

    def filter_TS(self) -> "FileParserBatch":
        return self.__new_batch(
            [
                parser
                for parser in self
                if parser.__class__ in qm_parsers and parser[-1].is_TS
            ]
        )

    def filter_error(self) -> "FileParserBatch":
        """
        Return a new `FileParserBatch` with all the QM parsers that are flagged as errors.

        Returns:
            FileParserBatch: The new `FileParserBatch`.
        """
        return self.__new_batch(
            [
                parser
                for parser in self
                if parser.__class__ in qm_parsers and parser[-1].is_error
            ]
        )

    def filter_normal(self) -> "FileParserBatch":
        """
        Return a new `FileParserBatch` with all the QM parsers that are flagged as normal.

        Returns:
            FileParserBatch: The new `FileParserBatch`.
        """
        return self.__new_batch(
            [
                parser
                for parser in self
                if parser.__class__ in qm_parsers and parser[-1].is_normal
            ]
        )

    def filter_by_charge(self, charge: int) -> "FileParserBatch":
        return self.__new_batch(
            [parser for parser in self if parser[-1].charge == charge]
        )

    def filter_by_multi(self, multi: int) -> "FileParserBatch":
        return self.__new_batch(
            [parser for parser in self if parser[-1].multiplicity == multi]
        )

    def filter_by_format(self, format: str) -> "FileParserBatch":
        if not format.startswith("."):
            format = "." + format
        return self.__new_batch(
            [parser for parser in self if parser.file_format == format]
        )

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({len(self)})"

    def hong_style_summary(self, frameID=-1) -> pd.DataFrame:
        empty_df = pd.DataFrame(
            [],
            columns=[
                "file_name",
                "charge",
                "multiplicity",
                "Canonical SMILES",
                "opt_software",
                "opt_version",
                "opt_method",
                "opt_functional",
                "opt_basis",
                "opt_solvent_model",
                "opt_solvent",
                "opt_normal terminated",
                "sp_software",
                "sp_version",
                "sp_method",
                "sp_functional",
                "sp_basis",
                "sp_solvent_model",
                "sp_solvent",
                "sp_normal terminated",
                "SP(hartree)",
                "ZPE(kcal/mol)",
                "TCH(kcal/mol)",
                "TCG(kcal/mol)",
                "H-Opt(kcal/mol)",
                "G-Opt(kcal/mol)",
                "Freq1(cm-1)",
                "Freq2(cm-1)",
                "H-Sum(kcal/mol)",
                "G-Sum(kcal/mol)",
            ],
        )
        total_df = pd.concat(
            [parser[frameID].hong_style_summary_series() for parser in self], axis=1
        ).T
        if "sp_normal terminated" in total_df.columns:
            sp_sub_df_index = total_df.apply(
                lambda x: x["file_name"].lower().endswith("sp"), axis=1
            )
            sp_sub_df = total_df[
                sp_sub_df_index & (total_df["sp_normal terminated"] == True)
            ]
            sp_error_terminated_df = total_df[
                sp_sub_df_index & (total_df["sp_normal terminated"] == False)
            ]
            if sp_error_terminated_df.shape[0] > 0:
                moloplogger.warning(
                    f"The following files have terminated with error in SP calculation:\n{sp_error_terminated_df['file_name'].tolist()}"
                )
        else:
            sp_sub_df = empty_df.copy()
        if "opt_normal terminated" in total_df.columns:
            opt_sub_df_index = total_df.apply(
                lambda x: x["file_name"].lower().endswith("opt"), axis=1
            )
            opt_sub_df = total_df[
                opt_sub_df_index & (total_df["opt_normal terminated"] == True)
            ]
            opt_error_terminated_df = total_df[
                opt_sub_df_index & (total_df["opt_normal terminated"] == False)
            ]
            if opt_error_terminated_df.shape[0] > 0:
                moloplogger.warning(
                    f"The following files have terminated with error in OPT calculation:\n{opt_error_terminated_df['file_name'].tolist()}"
                )
        else:
            opt_sub_df = empty_df.copy()
        if opt_sub_df.shape[0] != 0:
            opt_sub_df.loc[:, "file_name"] = opt_sub_df["file_name"].map(
                lambda x: x[:-4]
            )
        if sp_sub_df.shape[0] != 0:
            sp_sub_df.loc[:, "file_name"] = sp_sub_df["file_name"].map(lambda x: x[:-3])

        def clac_sum(row: pd.Series):
            if (
                pd.isna(row["SP(hartree)"])
                or pd.isna(row["TCH(kcal/mol)"])
                or pd.isna(row["TCG(kcal/mol)"])
            ):
                return row
            else:
                row["H-Sum(kcal/mol)"] = (
                    row["SP(hartree)"] * 627.5095 + row["TCH(kcal/mol)"]
                )
                row["G-Sum(kcal/mol)"] = (
                    row["SP(hartree)"] * 627.5095 + row["TCG(kcal/mol)"]
                )
                return row

        if opt_sub_df.shape[0] != 0 and sp_sub_df.shape[0] != 0:
            df = pd.merge(
                pd.concat(
                    [
                        sp_sub_df.iloc[:, :4],
                        sp_sub_df.iloc[:, 4:].dropna(how="all", axis=1),
                    ],
                    axis=1,
                ),
                pd.concat(
                    [
                        opt_sub_df.iloc[:, :4],
                        opt_sub_df.iloc[:, 4:]
                        .drop("SP(hartree)", axis=1)
                        .dropna(how="all", axis=1),
                    ],
                    axis=1,
                ),
                on=["file_name", "charge", "multiplicity", "Canonical SMILES"],
                how="outer",
            )
        elif opt_sub_df.shape[0] != 0 and sp_sub_df.shape[0] == 0:
            df = pd.concat(
                [
                    opt_sub_df,
                    empty_df[
                        [
                            "sp_software",
                            "sp_version",
                            "sp_method",
                            "sp_functional",
                            "sp_basis",
                            "sp_solvent_model",
                            "sp_solvent",
                            "sp_normal terminated",
                        ]
                    ],
                ],
                axis=1,
            )
        elif opt_sub_df.shape[0] == 0 and sp_sub_df.shape[0] != 0:
            df = pd.concat(
                [
                    sp_sub_df,
                    empty_df[
                        [
                            "opt_software",
                            "opt_version",
                            "opt_method",
                            "opt_functional",
                            "opt_basis",
                            "opt_solvent_model",
                            "opt_solvent",
                            "opt_normal terminated",
                            "ZPE(kcal/mol)",
                            "TCH(kcal/mol)",
                            "TCG(kcal/mol)",
                            "H-Opt(kcal/mol)",
                            "G-Opt(kcal/mol)",
                            "Freq1(cm-1)",
                            "Freq2(cm-1)",
                            "H-Sum(kcal/mol)",
                            "G-Sum(kcal/mol)",
                        ]
                    ],
                ],
                axis=1,
            )
        else:
            df = empty_df.copy()
        return (
            df.apply(clac_sum, axis=1)
            .reset_index(drop=True)[
                [
                    "file_name",
                    "charge",
                    "multiplicity",
                    "Canonical SMILES",
                    "opt_software",
                    "opt_version",
                    "opt_method",
                    "opt_functional",
                    "opt_basis",
                    "opt_solvent_model",
                    "opt_solvent",
                    "opt_normal terminated",
                    "sp_software",
                    "sp_version",
                    "sp_method",
                    "sp_functional",
                    "sp_basis",
                    "sp_solvent_model",
                    "sp_solvent",
                    "sp_normal terminated",
                    "SP(hartree)",
                    "ZPE(kcal/mol)",
                    "TCH(kcal/mol)",
                    "TCG(kcal/mol)",
                    "H-Opt(kcal/mol)",
                    "G-Opt(kcal/mol)",
                    "Freq1(cm-1)",
                    "Freq2(cm-1)",
                    "H-Sum(kcal/mol)",
                    "G-Sum(kcal/mol)",
                ]
            ]
            .sort_values(["file_name"])
            .reset_index(drop=True)
        )

    def to_summary_df(
        self,
        full: bool = False,
        with_units: bool = True,
        *,
        frameID=-1,
        use_hong_style: bool = False,
    ) -> pd.DataFrame:
        if use_hong_style:
            return self.hong_style_summary(frameID)
        return pd.concat(
            [parser[frameID].to_summary_series(full, with_units) for parser in self],
            axis=1,
        ).T

    def to_summary_csv(
        self,
        file_path: str = None,
        full: bool = False,
        with_units: bool = True,
        *,
        frameID=-1,
        use_hong_style: bool = False,
    ):
        if not file_path:
            file_path = os.path.join(os.path.curdir, "summary.csv")
        if os.path.isdir(file_path):
            file_path = os.path.join(file_path, "summary.csv")
        self.to_summary_df(
            full, with_units, frameID=frameID, use_hong_style=use_hong_style
        ).to_csv(file_path)
        moloplogger.info(f"summary csv saved to {os.path.abspath(file_path)}")

    def to_summary_excel(
        self,
        file_path: str = None,
        full: bool = False,
        with_units: bool = True,
        *,
        frameID=-1,
        ues_hong_style: bool = False,
    ):
        if not file_path:
            file_path = os.path.join(os.path.curdir, "summary.xlsx")
        if os.path.isdir(file_path):
            file_path = os.path.join(file_path, "summary.xlsx")
        self.to_summary_df(
            full, with_units, frameID=frameID, use_hong_style=ues_hong_style
        ).to_excel(file_path)
        moloplogger.info(f"summary xlsx saved to {os.path.abspath(file_path)}")

    def geometry_analysis(
        self,
        key_atoms: List[List[int]],
        file_path: str = None,
        precision: int = 1,
        one_start=False,
    ) -> List[str]:
        """
        Get the geometry infos among the atoms with all frames in each file, and save them to seperated csv files.

        Parameters:
            key_atoms L(ist[List[int]]):
                A list of list of index of the atoms, starts from 0
                    If the length of atom_idxs is 2, the bond length with unit Angstrom between the two atoms will be returned.

                    If the length of atom_idxs is 3, the angle with unit degree between  the three atoms will be returned.

                    If the length of atom_idxs is 4, the dihedral angle with unit degree between the four atoms will be returned.
            file_path (str):
                The path of the csv file to be saved. If None, the file will be saved in the same directory of the file_path.
            precision (int):
                The precision of the geometry analysis. Default is 1. e.g. 1 means 1.0001 ==> 1.0
            one_start (bool):
                If true, consider atom index starts from 1, so let index value subtracts 1 for all the atoms
        Returns:
            List[str]: file_paths
        """
        if file_path is None:
            file_path = os.path.curdir
        return [
            parser.geometry_analysis(
                key_atoms,
                os.path.join(file_path, parser.file_name),
                precision,
                one_start,
            )
            for parser in self
        ]

    def recover_structures(self, n_jobs: int = -1) -> List[str]:
        """
        Recover the structures of the QM parsers.
        """
        _n_jobs = self._set_n_jobs(n_jobs)
        total_smiles_list = []
        if _n_jobs > 1 and len(self) > _n_jobs:
            if molopconfig.show_progress_bar:
                for smiles_list in tqdm(
                    Parallel(
                        return_as="generator",
                        n_jobs=_n_jobs,
                        pre_dispatch="1.5*n_jobs",
                    )(delayed(parser.recover_structures)() for parser in self),
                    total=len(self),
                    desc=f"MolOP structure recovery with {_n_jobs} jobs",
                ):
                    total_smiles_list.extend(smiles_list)
            else:
                for smiles_list in Parallel(
                    return_as="generator",
                    n_jobs=self.__n_jobs,
                    pre_dispatch="1.5*n_jobs",
                )(delayed(parser.recover_structures)() for parser in self):
                    total_smiles_list.extend(smiles_list)
        else:
            if molopconfig.show_progress_bar:
                for parser in tqdm(
                    self,
                    total=len(self),
                    desc=f"MolOP structure recovery with single jobs",
                ):
                    total_smiles_list.extend(parser.recover_structures())
            else:
                for parser in self:
                    total_smiles_list.extend(parser.recover_structures())
        return total_smiles_list

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

    @property
    def closest_optimized_frames(self):
        """
        Find the closest optimized frames for each file.

        Returns:
            A list of optimized frame parsers closest to the optimized state.
        """
        return [parser.closest_optimized_frame for parser in self]
