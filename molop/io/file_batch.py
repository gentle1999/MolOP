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
from typing import Dict, Generator, List, Union, Sequence

import pandas as pd
from joblib import Parallel, cpu_count, delayed
from tqdm import tqdm

from molop.config import molopconfig
from molop.io.bases.file_base import BaseBlockParser, BaseFileParser
from molop.io.coords_file.gjf_parser import GJFParser
from molop.io.coords_file.sdf_parser import SDFBlockParser, SDFParser
from molop.io.coords_file.xyz_parser import XYZParser
from molop.io.qm_file.g16fchk_parser import G16FCHKParser
from molop.io.qm_file.g16log_parser import G16LOGParser
from molop.io.qm_file.xtbout_parser import XTBOUTParser
from molop.io.types import PARSERTYPES
from molop.logger.logger import logger

parsers: Dict[str, PARSERTYPES] = {
    ".gjf": (GJFParser,),
    ".gau": (GJFParser,),
    ".com": (GJFParser,),
    ".gjc": (GJFParser,),
    ".log": (G16LOGParser,),
    ".g16": (G16LOGParser,),
    ".gal": (G16LOGParser,),
    ".xyz": (XYZParser,),
    ".sdf": (SDFParser,),
    ".mol": (SDFParser,),
    ".out": (G16LOGParser, XTBOUTParser),
    ".irc": (G16LOGParser,),
    ".fchk": (G16FCHKParser,),
    ".fck": (G16FCHKParser,),
    ".fch": (G16FCHKParser,),
}

qm_parsers = (G16LOGParser, G16FCHKParser, XTBOUTParser)


def singlefile_parser(
    file_path: str,
    charge=None,
    multiplicity=None,
    only_extract_structure=False,
    only_last_frame=False,
) -> PARSERTYPES:
    """
    A function to parse a single file based on its format.

    Args:
        file_path (str): The path of the file.
        charge (int, optional): The molecular charge. Defaults to None.
        multiplicity (int, optional): The molecular spin multiplicity. Defaults to None.
        only_extract_structure (bool, optional): Whether to extract structure only. Defaults to False.
        only_last_frame (bool, optional): Whether to extract the last frame only. Defaults to False.

    Returns:
        PARSERTYPES: An instance of the appropriate parser class.

    Core Logic:
        - Checks if the provided `file_path` is an actual file.
        - Determines the file format by getting the file extension.
        - If the format is not supported, logs an error and returns None.
        - Iterates through the available parsers for that format.
        - Tries to instantiate and return the parser object with the given parameters.
        - If an exception occurs during parsing, it logs the error and tries the next parser in the list.
        - If all parsers fail, logs the final error and returns None.
        - If the `file_path` points to a directory, logs an error and returns None.

    Raises:
        None

    """
    if os.path.isfile(file_path):
        _, file_format = os.path.splitext(file_path)
        if file_format not in parsers:
            logger.error("Unknown file format: {}, {}".format(file_format, file_path))
            return None
        for idx, parser in enumerate(parsers[file_format]):
            try:
                # Instantiate and return specific parser classes with their respective arguments
                if parser in qm_parsers:
                    p = parser(
                        file_path,
                        charge=charge,
                        multiplicity=multiplicity,
                        only_extract_structure=only_extract_structure,
                        only_last_frame=only_last_frame,
                    )
                    return p
                elif parser in (GJFParser,):
                    return parser(file_path, charge=charge, multiplicity=multiplicity)
                elif parser in (XYZParser,):
                    return parser(
                        file_path,
                        charge=charge,
                        multiplicity=multiplicity,
                        only_last_frame=only_last_frame,
                    )
                elif parser in (SDFParser,):
                    return parser(file_path, only_last_frame=only_last_frame)
                else:
                    return parser(file_path)
            except Exception as e:
                if idx == len(parsers[file_format]) - 1:
                    logger.error(
                        f"Failed to parse file {file_path} with {parser.__name__}. {e}"
                    )
                    return None
                logger.debug(
                    f"Failed to parse file {file_path} with {parser.__name__}, trying {parsers[file_format][idx+1].__name__} instead"
                )
    elif os.path.isdir(file_path):
        logger.error(f"{file_path} is not a file.")
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

    def __getitem__(self, key: Union[int, str, slice]) -> PARSERTYPES:
        if isinstance(key, int):
            return list(self.__parsers.values())[key]
        if isinstance(key, slice):
            slicedkeys = list(self.__parsers.keys())[key]
            return self.__new_batch([self.__parsers[k] for k in slicedkeys])
        else:
            return self.__parsers[key]

    def __setitem__(self, key: str, value: PARSERTYPES):
        self.__parsers[key] = value

    def __delitem__(self, key: str):
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
        charge=None,
        multiplicity=None,
        only_extract_structure=False,
        only_last_frame=False,
    ):
        file_paths = []
        for file_path in files:
            if not os.path.isfile(file_path):
                logger.warning(f"{file_path} is not a file.")
                continue
            if os.path.splitext(file_path)[1] not in parsers:
                logger.warning(f"Unsupported input file format: {file_path}")
                continue
            if file_path.endswith("molop.log"):
                continue
            file_paths.append(os.path.abspath(file_path))

        arguments_list = [
            {
                "file_path": file_path,
                "charge": charge,
                "multiplicity": multiplicity,
                "only_extract_structure": only_extract_structure,
                "only_last_frame": only_last_frame,
            }
            for file_path in file_paths
        ]
        if self.__n_jobs > 1 and len(file_paths) > self.__n_jobs:
            if molopconfig.show_progress_bar:
                self.__parsers.update(
                    {
                        parser._file_path: parser
                        for parser in tqdm(
                            Parallel(
                                return_as="generator",
                                n_jobs=self.__n_jobs,
                                pre_dispatch="1.5*n_jobs",
                            )(
                                delayed(singlefile_parser)(**arguments)
                                for arguments in arguments_list
                            ),
                            total=len(arguments_list),
                            desc=f"MolOP parsing with {self.__n_jobs} jobs",
                        )
                        if parser is not None and len(parser) > 0
                    }
                )
            else:
                self.__parsers.update(
                    {
                        parser._file_path: parser
                        for parser in Parallel(
                            return_as="generator",
                            n_jobs=self.__n_jobs,
                            pre_dispatch="1.5*n_jobs",
                        )(
                            delayed(singlefile_parser)(**arguments)
                            for arguments in arguments_list
                        )
                        if parser is not None and len(parser) > 0
                    }
                )
        else:
            if molopconfig.show_progress_bar:
                self.__parsers.update(
                    {
                        parser._file_path: parser
                        for parser in tqdm(
                            (
                                singlefile_parser(**arguments)
                                for arguments in arguments_list
                            ),
                            total=len(arguments_list),
                            desc=f"MolOP parsing with single thread",
                        )
                        if parser is not None and len(parser) > 0
                    }
                )
            else:
                self.__parsers.update(
                    {
                        parser._file_path: parser
                        for parser in (
                            singlefile_parser(**arguments)
                            for arguments in arguments_list
                        )
                        if parser is not None and len(parser) > 0
                    }
                )

        logger.info(
            f"{len(file_paths) - len(self.__parsers)} files failed to parse, {len(self.__parsers)} successfully parsed"
        )

    def add_file_parsers(self, file_parsers: List[PARSERTYPES]) -> None:
        for parser in file_parsers:
            if not isinstance(parser, BaseFileParser):
                raise TypeError(
                    f"file_parsers should be a list of BaseFileParser, got {type(parser)}"
                )
            if parser.file_path not in self.__parsers:
                self[parser.file_path] = parser
            else:
                logger.warning(
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
        prefix: str = "#p opt b3lyp def2svp freq EmpiricalDispersion=GD3BJ NoSymm\n",
        suffix="\n\n",
        template: str = None,
        chk: bool = True,
        oldchk: bool = False,
        frameID=-1,
    ) -> None:
        """file_path should be a directory, prefix and suffix are optional"""
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
        logger.info(f"gjf files saved to {os.path.abspath(file_path)}")

    def to_XYZ_file(self, file_path: str = None) -> None:
        if file_path is None:
            file_path = os.path.curdir
        if not os.path.isdir(file_path):
            raise NotADirectoryError(f"{file_path} is not a directory")
        for parser in self:
            parser.to_XYZ_file(os.path.join(file_path, parser.pure_filename + ".xyz"))
        logger.info(f"xyz files saved to {os.path.abspath(file_path)}")

    def to_SDF_file(self, file_path: str = None) -> None:
        if file_path is None:
            file_path = os.path.curdir
        if not os.path.isdir(file_path):
            raise NotADirectoryError(f"{file_path} is not a directory")
        for parser in self:
            parser.to_SDF_file(os.path.join(file_path, parser.pure_filename + ".sdf"))
        logger.info(f"sdf files saved to {os.path.abspath(file_path)}")

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
        logger.info(f"chemdraw files saved to {os.path.abspath(file_path)}")

    def replace_substituent(
        self,
        query_smi: str,
        replacement_smi: str,
        bind_idx: int = None,
        replace_all: bool = False,
        attempt_num: int = 10,
        frameID=-1,
    ) -> "FileParserBatch":
        """
        Replace the substituent with the given SMARTS. The substituent is defined by the query_smi, and the new substituent is defined by the replacement_smi.

        Parameters:
            query_smi str:
                The SMARTS to query the substituent in the original molecule.
            replacement_smi str:
                The SMARTS of new substituent. The bind atom is the first atom of the replacement_smi.
            bind_idx int:
                The index of the atom to bind the new substituent. The default is None, which means to replace the first legal atom in original molecule.
                If specified, try to replace the atom. User should meke sure the atom is legal.
                Detail example in (Repalce Substituent)[Repalce Substituent]
            replace_all bool:
                If True, replace all the substituent queried in the original molecule.
            attempt_num int:
                Max attempt times to replace the substituent. Each time a new substituent conformation will be used for substitution.
            frameID int:
                The frameID to replace.
        Returns:
            The new `FileParserBatch`.
        """
        new_parsers = []
        if molopconfig.show_progress_bar:
            for parser in tqdm(self):
                try:
                    temp_parser = parser[frameID].replace_substituent(
                        query_smi=query_smi,
                        replacement_smi=replacement_smi,
                        bind_idx=bind_idx,
                        replace_all=replace_all,
                        attempt_num=attempt_num,
                    )
                    temp_file_path = temp_parser.to_SDF_file()
                    new_parsers.append(SDFParser(temp_file_path, only_last_frame=True))
                    os.remove(temp_file_path)
                except:
                    logger.error(
                        f"Failed to replace substituent from {query_smi} to {replacement_smi} in {parser.file_path}, {parser.file_name}"
                    )
        else:
            for parser in self:
                try:
                    temp_parser = parser[frameID].replace_substituent(
                        query_smi=query_smi,
                        replacement_smi=replacement_smi,
                        bind_idx=bind_idx,
                        replace_all=replace_all,
                        attempt_num=attempt_num,
                    )
                    temp_file_path = temp_parser.to_SDF_file()
                    new_parsers.append(SDFParser(temp_file_path, only_last_frame=True))
                    os.remove(temp_file_path)
                except:
                    logger.error(
                        f"Failed to replace substituent from {query_smi} to {replacement_smi} in {parser.file_path}, {parser.file_name}"
                    )

        new_batch = FileParserBatch()
        new_batch.add_file_parsers(new_parsers)
        logger.info(
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
                    new_parsers.append(SDFParser(temp_file_path, only_last_frame=True))
                    os.remove(temp_file_path)
                except Exception as e:
                    logger.error(
                        f"{e}: Failed to reset atom index by {mapping_smarts} in {parser.file_path}, {parser.file_name}"
                    )
        else:
            for parser in self:
                try:
                    temp_parser = parser[frameID].reset_atom_index(
                        mapping_smarts=mapping_smarts
                    )
                    temp_file_path = temp_parser.to_SDF_file()
                    new_parsers.append(SDFParser(temp_file_path, only_last_frame=True))
                    os.remove(temp_file_path)
                except Exception as e:
                    logger.error(
                        f"{e}: Failed to reset atom index by {mapping_smarts} in {parser.file_path}, {parser.file_name}"
                    )
        new_batch = FileParserBatch()
        new_batch.add_file_parsers(new_parsers)
        logger.info(
            f"{len(new_parsers)} files successfully replaced, {len(self.__parsers) - len(new_parsers)} files failed to reset_index"
        )
        return new_batch

    def standard_orient(
        self,
        anchor_list: List[int],
        frameID: int = -1,
    ) -> "FileParserBatch":
        """
        Depending on the input `idx_list`, `translate_anchor`, `rotate_anchor_to_X`, and `rotate_anchor_to_XY` are executed in order to obtain the normalized oriented molecule.

        Sub-functions:
            - `translate_anchor`: Translate the entire molecule so that the specified atom reaches the origin.
            - `rotate_anchor_to_X`: Rotate the specified second atom along the axis passing through the origin so that it reaches the positive half-axis of the X-axis.
            - `rotate_anchor_to_XY`: Rotate along the axis passing through the origin so that the specified third atom reaches quadrant 1 or 2 of the XY plane.

        Parameters:
            anchor_list List[int]:
                A list of indices of the atoms to be translated to origin, rotated to X axis, and rotated again to XY face:

                - If length is 1, execute `translate_anchor`
                - If length is 2, execute `translate_anchor` and `rotate_anchor_to_X`
                - If length is 3, execute `translate_anchor`, `rotate_anchor_to_X` and `rotate_anchor_to_XY`
                - If the length of the input `idx_list` is greater than 3, subsequent atomic numbers are not considered.


        Returns:
            The new `FileParserBatch`.
        """
        new_parsers = []
        if molopconfig.show_progress_bar:
            for parser in tqdm(self):
                try:
                    temp_parser = parser[frameID].standard_orient(anchor_list)
                    temp_file_path = temp_parser.to_SDF_file()
                    new_parsers.append(SDFParser(temp_file_path, only_last_frame=True))
                    os.remove(temp_file_path)
                except Exception as e:
                    logger.error(
                        f"{e}: Failed to standard_orient by {anchor_list} in {parser.file_path}, {parser.file_name}"
                    )
        else:
            for parser in self:
                try:
                    temp_parser = parser[frameID].standard_orient(anchor_list)
                    temp_file_path = temp_parser.to_SDF_file()
                    new_parsers.append(SDFParser(temp_file_path, only_last_frame=True))
                    os.remove(temp_file_path)
                except Exception as e:
                    logger.error(
                        f"{e}: Failed to standard_orient by {anchor_list} in {parser.file_path}, {parser.file_name}"
                    )
        new_batch = FileParserBatch()
        new_batch.add_file_parsers(new_parsers)
        logger.info(
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
                if parser.__class__ in qm_parsers and parser[-1].is_TS()
            ]
        )

    def filter_error(self) -> "FileParserBatch":
        """
        Return a new `FileParserBatch` with all the QM parsers that are flagged as errors.

        Returns:
            The new `FileParserBatch`.
        """
        return self.__new_batch(
            [
                parser
                for parser in self
                if parser.__class__ in qm_parsers and parser[-1].is_error()
            ]
        )

    def filter_normal(self) -> "FileParserBatch":
        """
        Return a new `FileParserBatch` with all the QM parsers that are flagged as normal.

        Returns:
            The new `FileParserBatch`.
        """
        return self.__new_batch(
            [
                parser
                for parser in self
                if parser.__class__ in qm_parsers and not parser[-1].is_error()
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
            [parser for parser in self if parser._file_format == format]
        )

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({len(self)})"

    def to_summary_df(self):
        return pd.concat([parser[-1].to_summary_series() for parser in self], axis=1).T

    def to_summary_csv(self, file_path: str = None):
        if not file_path:
            file_path = os.path.join(os.path.curdir, "summary.csv")
        if os.path.isdir(file_path):
            file_path = os.path.join(file_path, "summary.csv")
        self.to_summary_df().to_csv(file_path)
        logger.info(f"summary csv saved to {os.path.abspath(file_path)}")

    def to_summary_excel(self, file_path: str = None):
        if not file_path:
            file_path = os.path.join(os.path.curdir, "summary.xlsx")
        if os.path.isdir(file_path):
            file_path = os.path.join(file_path, "summary.xlsx")
        self.to_summary_df().to_excel(file_path)
        logger.info(f"summary xlsx saved to {os.path.abspath(file_path)}")

    def geometry_analysis(
        self,
        key_atoms: List[List[int]],
        file_path: str = None,
        precision: int = 1,
        one_start=False,
    ):
        """
        Get the geometry infos among the atoms with all frames in each file, and save them to seperated csv files.

        Parameters:
            key_atoms List[List[int]]:
                A list of list of index of the atoms, starts from 0
                    If the length of atom_idxs is 2, the bond length with unit Angstrom between the two atoms will be returned.

                    If the length of atom_idxs is 3, the angle with unit degree between  the three atoms will be returned.

                    If the length of atom_idxs is 4, the dihedral angle with unit degree between the four atoms will be returned.
            file_path str:
                The path of the csv file to be saved. If None, the file will be saved in the same directory of the file_path.
            precision int:
                The precision of the geometry analysis. Default is 1. e.g. 1 means 1.0001 ==> 1.0
            one_start bool:
                If true, consider atom index starts from 1, so let index value subtracts 1 for all the atoms
        Returns:
            file_paths
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
