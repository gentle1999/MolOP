"""
Author: TMJ
Date: 2024-01-25 22:46:35
LastEditors: TMJ
LastEditTime: 2024-01-25 23:10:59
Description: 请填写简介
"""
import os
from typing import List, Union
from molop.config import molopconfig

import pandas as pd
from joblib import Parallel, cpu_count, delayed
from tqdm import tqdm

from molop.config import molopconfig
from molop.io.bases.file_base import BaseFileParser, BlockType
from molop.io.coords_file.gjf_parser import GJFParser
from molop.io.coords_file.sdf_parser import SDFBlockParser, SDFParser
from molop.io.coords_file.xyz_parser import XYZParser
from molop.io.qm_file.g16fchk_parser import G16FCHKParser
from molop.io.qm_file.g16irc_parser import G16IRCParser
from molop.io.qm_file.g16log_parser import G16LOGParser
from molop.io.qm_file.xtbout_parser import XTBOUTParser
from molop.logger.logger import logger

PARSERTYPES = Union[
    GJFParser,
    SDFParser,
    XYZParser,
    G16FCHKParser,
    G16IRCParser,
    G16LOGParser,
    XTBOUTParser,
]

parsers = {
    ".gjf": (GJFParser,),
    ".log": (G16LOGParser,),
    ".xyz": (XYZParser,),
    ".sdf": (SDFParser,),
    ".out": (XTBOUTParser, G16IRCParser),
    ".irc": (G16IRCParser,),
    ".fchk": (G16FCHKParser,),
    ".fck": (G16FCHKParser,),
}

qm_parsers = (G16LOGParser, G16IRCParser, G16FCHKParser, XTBOUTParser)


def singlefile_parser(
    file_path,
    charge=None,
    multiplicity=None,
    only_extract_structure=False,
    only_last_frame=False,
) -> Union[PARSERTYPES, None]:
    if os.path.isfile(file_path):
        _, file_format = os.path.splitext(file_path)
        if file_format not in parsers:
            raise NotImplementedError("Unknown file format: {}".format(file_format))
        for idx, parser in enumerate(parsers[file_format]):
            try:
                if parser in (G16LOGParser, XTBOUTParser, G16IRCParser, G16FCHKParser):
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
            except:
                if idx == len(parsers[file_format]) - 1:
                    logger.error(
                        f"Failed to parse file {file_path} with {parser.__name__}"
                    )
                    raise Exception(f"Failed to parse file {file_path}")
                logger.debug(
                    f"Failed to parse file {file_path} with {parser.__name__}, try {parsers[file_format][idx+1].__name__} instead"
                )
    elif os.path.isdir(file_path):
        raise IsADirectoryError(f"{file_path} is not a file.")


class FileParserBatch:
    __n_jobs: int
    __file_paths: List[str]
    __parsers = List[PARSERTYPES]

    def __init__(
        self,
        files: List[str],
        charge=None,
        multiplicity=None,
        n_jobs: int = -1,
        only_extract_structure=False,
        only_last_frame=False,
    ) -> None:
        self.__n_jobs = n_jobs if n_jobs > 0 else cpu_count()
        self.__file_paths = []
        for file_path in files:
            if not os.path.isfile(file_path):
                raise FileNotFoundError(f"{file_path} is not a file.")
            if os.path.splitext(file_path)[1] not in parsers:
                logger.warning(f"Unsupported input file format: {file_path}")
                continue
            self.__file_paths.append(os.path.abspath(file_path))
        self.__parsers: List[PARSERTYPES] = []
        self._charge = charge
        self._multiplicity = multiplicity
        self._only_extract_structure = only_extract_structure
        self._only_last_frame = only_last_frame

        arguments_list = [
            {
                "file_path": file_path,
                "charge": self._charge,
                "multiplicity": self._multiplicity,
                "only_extract_structure": self._only_extract_structure,
                "only_last_frame": self._only_last_frame,
            }
            for file_path in self.__file_paths
        ]
        if self.__n_jobs > 1 and len(self.__file_paths) > 30:
            if molopconfig.show_progress_bar:
                self.__parsers = [
                    parser
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
                        desc=f"MolOP parsing with {cpu_count()} jobs",
                    )
                    if len(parser) > 0
                ]
            else:
                self.__parsers = [
                    parser
                    for parser in Parallel(
                        return_as="generator",
                        n_jobs=self.__n_jobs,
                        pre_dispatch="1.5*n_jobs",
                    )(
                        delayed(singlefile_parser)(**arguments)
                        for arguments in arguments_list
                    )
                    if len(parser) > 0
                ]
        else:
            if molopconfig.show_progress_bar:
                self.__parsers = [
                    parser
                    for parser in tqdm(
                        (
                            singlefile_parser(**arguments)
                            for arguments in arguments_list
                        ),
                        total=len(arguments_list),
                        desc=f"MolOP parsing with single thread",
                    )
                ]
            else:
                self.__parsers = [
                    parser
                    for parser in (
                        singlefile_parser(**arguments) for arguments in arguments_list
                    )
                ]

        logger.info(
            f"{len(self.__file_paths) - len(self.__parsers)} files failed to parse, {len(self.__parsers)} successfully parsed"
        )

    def to_GJF_file(
        self,
        file_path: str = None,
        charge: int = None,
        multiplicity: int = None,
        prefix: str = f"# g16 gjf \n",
        suffix="\n\n",
    ) -> None:
        """file_path should be a directory, prefix and suffix are optional"""
        for parser in self.__parsers:
            parser.to_GJF_file(
                os.path.join(file_path, parser.file_name),
                charge=charge,
                multiplicity=multiplicity,
                prefix=prefix,
                suffix=suffix,
            )

    def to_XYZ_file(self, file_path: str = None) -> None:
        if file_path is None:
            file_path = os.path.curdir
        for parser in self.__parsers:
            parser.to_XYZ_file(os.path.join(file_path, parser.file_name))

    def to_SDF_file(self, file_path: str = None) -> None:
        if file_path is None:
            file_path = os.path.curdir
        for parser in self.__parsers:
            parser.to_SDF_file(os.path.join(file_path, parser.file_name))

    def to_chemdraw(self, file_path: str = None, frameID=-1, keep3D=False) -> None:
        if file_path is None:
            file_path = os.path.curdir
        for parser in self.__parsers:
            parser.to_chemdraw(
                os.path.join(file_path, parser.file_name),
                frameID=frameID,
                keep3D=keep3D,
            )

    def replace_substituent(
        self,
        query_smi: str,
        replacement_smi: str,
        bind_idx: int = None,
        replace_all: bool = False,
        attempt_num: int = 10,
    ) -> List[BlockType]:
        """only consider the last frame of each file"""
        new_parsers = []
        if molopconfig.show_progress_bar:
            for parser in tqdm(self.__parsers):
                try:
                    temp_parser = parser[-1].replace_substituent(
                        query_smi=query_smi,
                        replacement_smi=replacement_smi,
                        bind_idx=bind_idx,
                        replace_all=replace_all,
                        attempt_num=attempt_num,
                    )
                    new_parsers.append(
                        SDFBlockParser(
                            temp_parser.to_SDF_block(),
                            os.path.splitext(temp_parser._file_path)[0] + ".sdf",
                        )
                    )
                except:
                    logger.error(
                        f"Failed to replace substituent from {query_smi} to {replacement_smi} in {parser.file_path}, {parser.file_name}"
                    )
        else:
            for parser in self.__parsers:
                try:
                    temp_parser = parser[-1].replace_substituent(
                        query_smi=query_smi,
                        replacement_smi=replacement_smi,
                        bind_idx=bind_idx,
                        replace_all=replace_all,
                        attempt_num=attempt_num,
                    )
                    new_parsers.append(
                        SDFBlockParser(
                            temp_parser.to_SDF_block(),
                            os.path.splitext(temp_parser._file_path)[0] + ".sdf",
                        )
                    )
                except:
                    logger.error(
                        f"Failed to replace substituent from {query_smi} to {replacement_smi} in {parser.file_path}, {parser.file_name}"
                    )

        logger.info(
            f"{len(new_parsers)} files successfully replaced, {len(self.__parsers) - len(new_parsers)} files failed to replace"
        )
        return new_parsers

    def reset_atom_index(self, mapping_smarts: str) -> List[BlockType]:
        """only consider the last frame of each file"""
        new_parsers = []
        if molopconfig.show_progress_bar:
            for parser in tqdm(self.__parsers):
                try:
                    temp_parser = parser[-1].reset_atom_index(
                        mapping_smarts=mapping_smarts
                    )
                    new_parsers.append(
                        SDFBlockParser(
                            temp_parser.to_SDF_block(),
                            os.path.splitext(temp_parser._file_path)[0] + ".sdf",
                        )
                    )
                except Exception as e:
                    logger.error(
                        f"{e}: Failed to reset atom index by {mapping_smarts} in {parser.file_path}, {parser.file_name}"
                    )
        else:
            for parser in self.__parsers:
                try:
                    temp_parser = parser[-1].reset_atom_index(
                        mapping_smarts=mapping_smarts
                    )
                    new_parsers.append(
                        SDFBlockParser(
                            temp_parser.to_SDF_block(),
                            os.path.splitext(temp_parser._file_path)[0] + ".sdf",
                        )
                    )
                except Exception as e:
                    logger.error(
                        f"{e}: Failed to reset atom index by {mapping_smarts} in {parser.file_path}, {parser.file_name}"
                    )
        logger.info(
            f"{len(new_parsers)} files successfully replaced, {len(self.__parsers) - len(new_parsers)} files failed to reset_index"
        )
        return new_parsers

    def __getitem__(self, parserID: int) -> PARSERTYPES:
        return self.__parsers[parserID]

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
            return self.__parsers[self.__index - 1]

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({len(self)})"

    def to_summary_df(self):
        def get_generator_value(g, idx):
            i = 0
            temp_g = g
            try:
                while i <= idx:
                    temp_g = next(g)
                    i += 1
                else:
                    return temp_g
            except StopIteration:
                return {"is imaginary": None, "freq": None}

        if self._only_extract_structure:
            return pd.DataFrame(
                {
                    "parser": [parser.__class__.__name__ for parser in self.__parsers],
                    "file_name": [parser.file_name for parser in self.__parsers],
                    "file_path": [parser.file_path for parser in self.__parsers],
                    "file_format": [parser._file_format for parser in self.__parsers],
                    "charge": [parser[-1].charge for parser in self.__parsers],
                    "multiplicity": [
                        parser[-1].multiplicity for parser in self.__parsers
                    ],
                    "SMILES": [
                        parser[-1].to_standard_SMILES() for parser in self.__parsers
                    ],
                }
            )
        else:
            return pd.DataFrame(
                {
                    "parser": [parser.__class__.__name__ for parser in self.__parsers],
                    "file_name": [parser.file_name for parser in self.__parsers],
                    "file_path": [parser.file_path for parser in self.__parsers],
                    "file_format": [parser._file_format for parser in self.__parsers],
                    "charge": [parser[-1].charge for parser in self.__parsers],
                    "multiplicity": [
                        parser[-1].multiplicity for parser in self.__parsers
                    ],
                    "SMILES": [
                        parser[-1].to_standard_SMILES() for parser in self.__parsers
                    ],
                    "status": [
                        parser[-1].state if parser.__class__ in qm_parsers else None
                        for parser in self.__parsers
                    ],
                    "ZPE": [
                        parser[-1].dimensionless_sum_energy["zero-point correction"]
                        if parser.__class__ in qm_parsers
                        else None
                        for parser in self.__parsers
                    ],
                    "TCE": [
                        parser[-1].dimensionless_sum_energy["TCE"]
                        if parser.__class__ in qm_parsers
                        else None
                        for parser in self.__parsers
                    ],
                    "TCH": [
                        parser[-1].dimensionless_sum_energy["TCH"]
                        if parser.__class__ in qm_parsers
                        else None
                        for parser in self.__parsers
                    ],
                    "TCG": [
                        parser[-1].dimensionless_sum_energy["TCG"]
                        if parser.__class__ in qm_parsers
                        else None
                        for parser in self.__parsers
                    ],
                    "ZPE-Gas": [
                        parser[-1].dimensionless_sum_energy["zero-point gas"]
                        if parser.__class__ in qm_parsers
                        else None
                        for parser in self.__parsers
                    ],
                    "E-Gas": [
                        parser[-1].dimensionless_sum_energy["E gas"]
                        if parser.__class__ in qm_parsers
                        else None
                        for parser in self.__parsers
                    ],
                    "H-Gas": [
                        parser[-1].dimensionless_sum_energy["H gas"]
                        if parser.__class__ in qm_parsers
                        else None
                        for parser in self.__parsers
                    ],
                    "G-Gas": [
                        parser[-1].dimensionless_sum_energy["G gas"]
                        if parser.__class__ in qm_parsers
                        else None
                        for parser in self.__parsers
                    ],
                    "sp": [
                        parser[-1].dimensionless_energy
                        if parser.__class__ in qm_parsers
                        else None
                        for parser in self.__parsers
                    ],
                    "HOMO": [
                        parser[-1].dimensionless_alpha_energy["homo"]
                        if parser.__class__ in qm_parsers
                        else None
                        for parser in self.__parsers
                    ],
                    "LUMO": [
                        parser[-1].dimensionless_alpha_energy["lumo"]
                        if parser.__class__ in qm_parsers
                        else None
                        for parser in self.__parsers
                    ],
                    "GAP": [
                        parser[-1].dimensionless_alpha_energy["gap"]
                        if parser.__class__ in qm_parsers
                        else None
                        for parser in self.__parsers
                    ],
                    "first freq": [
                        get_generator_value(parser[-1].dimensionless_frequencies, 0)[
                            "freq"
                        ]
                        if parser.__class__ in qm_parsers
                        else None
                        for parser in self.__parsers
                    ],
                    "first freq tag": [
                        get_generator_value(parser[-1].dimensionless_frequencies, 0)[
                            "is imaginary"
                        ]
                        if parser.__class__ in qm_parsers
                        else None
                        for parser in self.__parsers
                    ],
                    "second freq": [
                        get_generator_value(parser[-1].dimensionless_frequencies, 1)[
                            "freq"
                        ]
                        if parser.__class__ in qm_parsers
                        else None
                        for parser in self.__parsers
                    ],
                    "second freq tag": [
                        get_generator_value(parser[-1].dimensionless_frequencies, 1)[
                            "is imaginary"
                        ]
                        if parser.__class__ in qm_parsers
                        else None
                        for parser in self.__parsers
                    ],
                    "S**2": [
                        parser[-1].spin_eigenvalue
                        if parser.__class__ in qm_parsers
                        else None
                        for parser in self.__parsers
                    ],
                    "S": [
                        parser[-1].spin_multiplicity
                        if parser.__class__ in qm_parsers
                        else None
                        for parser in self.__parsers
                    ],
                }
            )

    def to_summary_csv(self, file_path: str = None):
        if not file_path:
            file_path = os.path.join(os.path.curdir, "summary.csv")
        self.to_summary_df().to_csv(file_path)

    def to_summary_excel(self, file_path: str = None):
        if not file_path:
            file_path = os.path.join(os.path.curdir, "summary.xlsx")
        self.to_summary_df().to_excel(file_path)

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
            for parser in self.__parsers
        ]