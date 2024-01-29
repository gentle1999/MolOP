"""
Author: TMJ
Date: 2024-01-25 22:46:35
LastEditors: TMJ
LastEditTime: 2024-01-25 23:10:59
Description: 请填写简介
"""
import os
from typing import List, Union
import multiprocessing
from joblib import Parallel, delayed, cpu_count

import pandas as pd
from tqdm import tqdm

from molop.io.bases.file_base import BaseFileParser, BlockType
from molop.io.coords_file.gjf_parser import GJFParser
from molop.io.coords_file.sdf_parser import SDFParser, SDFBlockParser
from molop.io.coords_file.xyz_parser import XYZParser
from molop.io.qm_file.g16fchk_parser import G16FCHKParser
from molop.io.qm_file.g16irc_parser import G16IRCParser
from molop.io.qm_file.g16log_parser import G16LOGParser
from molop.io.qm_file.xtbout_parser import XTBOUTParser
from molop.logger.logger import logger

PARSERTPES = Union[
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
) -> Union[PARSERTPES, None]:
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
                elif parser in (XYZParser, GJFParser):
                    return parser(file_path, charge=charge, multiplicity=multiplicity)
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
    __path: str
    __file_names: List[str]
    __parsers = List[PARSERTPES]

    def __init__(
        self,
        dir_path: str,
        charge=None,
        multiplicity=None,
        only_extract_structure=False,
        only_last_frame=False,
    ) -> None:
        if os.path.isdir(dir_path):
            self.__path = dir_path
            self.__file_names = [
                file_name
                for file_name in os.listdir(dir_path)
                if os.path.splitext(file_name)[-1] in parsers
            ]
        else:
            raise ValueError(f"{dir_path} is not a directory")
        self.__parsers: List[PARSERTPES] = []
        self._charge = charge
        self._multiplicity = multiplicity
        self._only_extract_structure = only_extract_structure
        self._only_last_frame = only_last_frame

        if len(self.__file_names) > 100000:
            logger.warning(
                f"{len(self.__file_names)} files found in {self.__path}, maybe too many files"
            )
        arguments_list = [
            {
                "file_path": os.path.join(self.__path, file_path),
                "charge": self._charge,
                "multiplicity": self._multiplicity,
                "only_extract_structure": self._only_extract_structure,
                "only_last_frame": self._only_last_frame,
            }
            for file_path in self.__file_names
        ]

        self.__parsers = [
            parser
            for parser in tqdm(
                Parallel(
                    return_as="generator", n_jobs=cpu_count(), pre_dispatch="1.5*n_jobs"
                )(
                    delayed(singlefile_parser)(**arguments)
                    for arguments in arguments_list
                ),
                total=len(arguments_list),
                desc=f"MolOP parsing with {cpu_count()} jobs",
            )
            if len(parser) > 0
        ]

        logger.info(
            f"{self.__path}: {len(self.__file_names) - len(self.__parsers)} files failed to parse, {len(self.__parsers)} successfully parsed"
        )

    def to_GJF_file(
        self, file_path: str = None, prefix: str = None, suffix="\n\n"
    ) -> None:
        """file_path should be a directory, prefix and suffix are optional"""
        if not file_path:
            file_path = self.__path
        assert os.path.isdir(
            file_path
        ), f"file_path should be a directory, got {file_path}"
        for parser in self.__parsers:
            parser.to_GJF_file(
                os.path.join(file_path, parser.file_name), prefix=prefix, suffix=suffix
            )

    def to_XYZ_file(self, file_path: str = None) -> None:
        if not file_path:
            file_path = self.__path
        assert os.path.isdir(
            file_path
        ), f"file_path should be a directory, got {file_path}"
        for parser in self.__parsers:
            parser.to_XYZ_file(os.path.join(file_path, parser.file_name))

    def to_SDF_file(self, file_path: str = None) -> None:
        if not file_path:
            file_path = self.__path
        assert os.path.isdir(
            file_path
        ), f"file_path should be a directory, got {file_path}"
        for parser in self.__parsers:
            parser.to_SDF_file(os.path.join(file_path, parser.file_name))

    def to_chemdraw(self, file_path: str = None, frameID=-1, keep3D=False) -> None:
        if not file_path:
            file_path = self.__path
        assert os.path.isdir(
            file_path
        ), f"file_path should be a directory, got {file_path}"
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
        logger.info(
            f"{len(new_parsers)} files successfully replaced, {len(self.__parsers) - len(new_parsers)} files failed to replace"
        )
        return new_parsers

    def reset_atom_index(self, mapping_smarts: str) -> List[BlockType]:
        """only consider the last frame of each file"""
        new_parsers = []
        for parser in tqdm(self.__parsers):
            try:
                temp_parser = parser[-1].reset_atom_index(mapping_smarts=mapping_smarts)
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

    def __getitem__(self, parserID: int) -> PARSERTPES:
        return self.__parsers[parserID]

    def __len__(self) -> int:
        return len(self.__parsers)

    def __iter__(self):
        self.__index = 0
        return self

    def __next__(
        self,
    ) -> PARSERTPES:
        if self.__index >= len(self):
            raise StopIteration
        else:
            self.__index += 1
            return self.__parsers[self.__index - 1]

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self.__path})"

    def to_summary_df(self):
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
                    "zero-point": [
                        parser[-1].sum_energy["zero-point"]
                        if parser.__class__ in qm_parsers
                        else None
                        for parser in self.__parsers
                    ],
                    "H": [
                        parser[-1].sum_energy["thermal energy"]
                        if parser.__class__ in qm_parsers
                        else None
                        for parser in self.__parsers
                    ],
                    "S": [
                        parser[-1].sum_energy["thermal enthalpy"]
                        if parser.__class__ in qm_parsers
                        else None
                        for parser in self.__parsers
                    ],
                    "G": [
                        parser[-1].sum_energy["thermal gibbs free energy"]
                        if parser.__class__ in qm_parsers
                        else None
                        for parser in self.__parsers
                    ],
                    "zero-point correction": [
                        parser[-1].sum_energy["zero-point correction"]
                        if parser.__class__ in qm_parsers
                        else None
                        for parser in self.__parsers
                    ],
                    "H correction": [
                        parser[-1].sum_energy["thermal energy correction"]
                        if parser.__class__ in qm_parsers
                        else None
                        for parser in self.__parsers
                    ],
                    "S correction": [
                        parser[-1].sum_energy["thermal enthalpy correction"]
                        if parser.__class__ in qm_parsers
                        else None
                        for parser in self.__parsers
                    ],
                    "G correction": [
                        parser[-1].sum_energy["thermal gibbs free energy correction"]
                        if parser.__class__ in qm_parsers
                        else None
                        for parser in self.__parsers
                    ],
                    "total energy": [
                        parser[-1].energy if parser.__class__ in qm_parsers else None
                        for parser in self.__parsers
                    ],
                    "HOMO": [
                        parser[-1].alpha_energy["homo"]
                        if parser.__class__ in qm_parsers
                        else None
                        for parser in self.__parsers
                    ],
                    "LUMO": [
                        parser[-1].alpha_energy["lumo"]
                        if parser.__class__ in qm_parsers
                        else None
                        for parser in self.__parsers
                    ],
                    "GAP": [
                        parser[-1].alpha_energy["gap"]
                        if parser.__class__ in qm_parsers
                        else None
                        for parser in self.__parsers
                    ],
                }
            )

    def to_summary_csv(self, file_path: str = None):
        if not file_path:
            file_path = os.path.join(self.__path, "summary.csv")
        self.to_summary_df().to_csv(file_path)

    def to_summary_excel(self, file_path: str = None):
        if not file_path:
            file_path = os.path.join(self.__path, "summary.xlsx")
        self.to_summary_df().to_excel(file_path)
