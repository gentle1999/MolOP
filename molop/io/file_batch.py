"""
Author: TMJ
Date: 2024-01-25 22:46:35
LastEditors: TMJ
LastEditTime: 2024-01-25 23:10:59
Description: 请填写简介
"""
import os
from collections.abc import MutableMapping
from collections import OrderedDict
from typing import Dict, List, Union, Generator
from molop.config import molopconfig
import pandas as pd
from joblib import Parallel, cpu_count, delayed
from tqdm import tqdm

from molop.config import molopconfig
from molop.io.bases.file_base import BaseFileParser, BlockType, BaseBlockParser
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
    ".gau": (GJFParser,),
    ".com": (GJFParser,),
    ".gjc": (GJFParser,),
    ".log": (G16LOGParser,),
    ".g16": (G16LOGParser,),
    ".gal": (G16LOGParser,),
    ".xyz": (XYZParser,),
    ".sdf": (SDFParser,),
    ".mol": (SDFParser,),
    ".out": (XTBOUTParser, G16IRCParser),
    ".irc": (G16IRCParser,),
    ".fchk": (G16FCHKParser,),
    ".fck": (G16FCHKParser,),
    ".fch": (G16FCHKParser,),
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


class FileParserBatch(MutableMapping):
    __n_jobs: int
    __parsers: Dict[str, PARSERTYPES]

    def __init__(
        self,
        n_jobs: int = -1,
    ) -> None:
        self.__n_jobs = (
            min(n_jobs, cpu_count(), molopconfig.max_jobs)
            if n_jobs > 0
            else min(cpu_count(), molopconfig.max_jobs)
        )
        self.__parsers: Dict[str, PARSERTYPES] = OrderedDict()

    def __getitem__(self, key: Union[int, str]) -> PARSERTYPES:
        if isinstance(key, int):
            return list(self.__parsers.values())[key]
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
                raise FileNotFoundError(f"{file_path} is not a file.")
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
                        if len(parser) > 0
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
                        if len(parser) > 0
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
                        if len(parser) > 0
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
                        if len(parser) > 0
                    }
                )

        logger.info(
            f"{len(file_paths) - len(self.__parsers)} files failed to parse, {len(self.__parsers)} successfully parsed"
        )

    def add_file_parsers(self, file_parsers: List[PARSERTYPES]) -> None:
        for parser in file_parsers:
            if parser.file_path not in self.__parsers:
                self[parser.file_path] = parser
            else:
                logger.warning(
                    f"File {parser.file_path} already exists in the batch, skipped"
                )

    def to_GJF_file(
        self,
        file_path: str = None,
        charge: int = None,
        multiplicity: int = None,
        prefix: str = "# g16 gjf \n",
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
                os.path.join(file_path, parser.file_name),
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
            parser.to_XYZ_file(os.path.join(file_path, parser.file_name))
        logger.info(f"xyz files saved to {os.path.abspath(file_path)}")

    def to_SDF_file(self, file_path: str = None) -> None:
        if file_path is None:
            file_path = os.path.curdir
        if not os.path.isdir(file_path):
            raise NotADirectoryError(f"{file_path} is not a directory")
        for parser in self:
            parser.to_SDF_file(os.path.join(file_path, parser.file_name))
        logger.info(f"sdf files saved to {os.path.abspath(file_path)}")

    def to_chemdraw(self, file_path: str = None, frameID=-1, keep3D=True) -> None:
        if file_path is None:
            file_path = os.path.curdir
        if not os.path.isdir(file_path):
            raise NotADirectoryError(f"{file_path} is not a directory")
        for parser in self:
            parser.to_chemdraw(
                os.path.join(file_path, parser.file_name),
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
    ) -> List[BlockType]:
        """only consider the last frame of each file"""
        new_parsers = []
        if molopconfig.show_progress_bar:
            for parser in tqdm(self):
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
            for parser in self:
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

    def reset_atom_index(self, mapping_smarts: str) -> "FileParserBatch":
        """only consider the last frame of each file"""
        new_parsers = []
        if molopconfig.show_progress_bar:
            for parser in tqdm(self):
                try:
                    temp_parser = parser[-1].reset_atom_index(
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
                    temp_parser = parser[-1].reset_atom_index(
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

    def filter_TS(self) -> "FileParserBatch":
        TS_parsers = [
            parser
            for parser in self
            if parser.__class__ in qm_parsers and parser[-1].is_TS
        ]
        new_batch = FileParserBatch()
        new_batch.add_file_parsers(TS_parsers)
        return new_batch

    def filter_by_charge(self, charge: int) -> "FileParserBatch":
        TS_parsers = [parser for parser in self if parser[-1].charge == charge]
        new_batch = FileParserBatch()
        new_batch.add_file_parsers(TS_parsers)
        return new_batch

    def filter_by_multi(self, multi: int) -> "FileParserBatch":
        TS_parsers = [parser for parser in self if parser[-1].multiplicity == multi]
        new_batch = FileParserBatch()
        new_batch.add_file_parsers(TS_parsers)
        return new_batch

    def filter_by_format(self, format: str) -> "FileParserBatch":
        if not format.startswith("."):
            format = "." + format
        TS_parsers = [parser for parser in self if parser._file_format == format]
        new_batch = FileParserBatch()
        new_batch.add_file_parsers(TS_parsers)
        return new_batch

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

        return pd.DataFrame(
            {
                "parser": [parser.__class__.__name__ for parser in self],
                "file_name": [parser.file_name for parser in self],
                "file_path": [parser.file_path for parser in self],
                "file_format": [parser._file_format for parser in self],
                "charge": [parser[-1].charge for parser in self],
                "multiplicity": [parser[-1].multiplicity for parser in self],
                "SMILES": [parser[-1].to_standard_SMILES() for parser in self],
                "status": [
                    parser[-1].state if parser.__class__ in qm_parsers else None
                    for parser in self
                ],
                "ZPE": [
                    parser[-1].dimensionless_sum_energy["zero-point correction"]
                    if parser.__class__ in qm_parsers
                    else None
                    for parser in self
                ],
                "TCE": [
                    parser[-1].dimensionless_sum_energy["TCE"]
                    if parser.__class__ in qm_parsers
                    else None
                    for parser in self
                ],
                "TCH": [
                    parser[-1].dimensionless_sum_energy["TCH"]
                    if parser.__class__ in qm_parsers
                    else None
                    for parser in self
                ],
                "TCG": [
                    parser[-1].dimensionless_sum_energy["TCG"]
                    if parser.__class__ in qm_parsers
                    else None
                    for parser in self
                ],
                "ZPE-Gas": [
                    parser[-1].dimensionless_sum_energy["zero-point gas"]
                    if parser.__class__ in qm_parsers
                    else None
                    for parser in self
                ],
                "E-Gas": [
                    parser[-1].dimensionless_sum_energy["E gas"]
                    if parser.__class__ in qm_parsers
                    else None
                    for parser in self
                ],
                "H-Gas": [
                    parser[-1].dimensionless_sum_energy["H gas"]
                    if parser.__class__ in qm_parsers
                    else None
                    for parser in self
                ],
                "G-Gas": [
                    parser[-1].dimensionless_sum_energy["G gas"]
                    if parser.__class__ in qm_parsers
                    else None
                    for parser in self
                ],
                "sp": [
                    parser[-1].dimensionless_energy
                    if parser.__class__ in qm_parsers
                    else None
                    for parser in self
                ],
                "HOMO": [
                    parser[-1].dimensionless_alpha_energy["homo"]
                    if parser.__class__ in qm_parsers
                    else None
                    for parser in self
                ],
                "LUMO": [
                    parser[-1].dimensionless_alpha_energy["lumo"]
                    if parser.__class__ in qm_parsers
                    else None
                    for parser in self
                ],
                "GAP": [
                    parser[-1].dimensionless_alpha_energy["gap"]
                    if parser.__class__ in qm_parsers
                    else None
                    for parser in self
                ],
                "first freq": [
                    get_generator_value(parser[-1].dimensionless_frequencies, 0)["freq"]
                    if parser.__class__ in qm_parsers
                    else None
                    for parser in self
                ],
                "first freq tag": [
                    get_generator_value(parser[-1].dimensionless_frequencies, 0)[
                        "is imaginary"
                    ]
                    if parser.__class__ in qm_parsers
                    else None
                    for parser in self
                ],
                "second freq": [
                    get_generator_value(parser[-1].dimensionless_frequencies, 1)["freq"]
                    if parser.__class__ in qm_parsers
                    else None
                    for parser in self
                ],
                "second freq tag": [
                    get_generator_value(parser[-1].dimensionless_frequencies, 1)[
                        "is imaginary"
                    ]
                    if parser.__class__ in qm_parsers
                    else None
                    for parser in self
                ],
                "S**2": [
                    parser[-1].spin_eigenvalue
                    if parser.__class__ in qm_parsers
                    else None
                    for parser in self
                ],
                "S": [
                    parser[-1].spin_multiplicity
                    if parser.__class__ in qm_parsers
                    else None
                    for parser in self
                ],
            }
        )

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

    def replace_substituent(
        self,
        query_smi: str,
        replacement_smi: str,
        bind_idx: int = None,
        replace_all=False,
        attempt_num: int = 10,
    ) -> Generator[BaseBlockParser, None, None]:
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

        Returns:
            The new parser.
        """
        return (
            f.replace_substituent(
                query_smi, replacement_smi, bind_idx, replace_all, attempt_num
            )
            for f in self
        )

    def reset_atom_index(
        self,
        mapping_smarts: str,
    ) -> Generator[BaseBlockParser, None, None]:
        """
        Reset the atom index of the molecule according to the mapping SMARTS.

        Parameters:
            mapping_smarts str:
                The SMARTS to query the molecule substructure.
                The queried atoms will be renumbered and placed at the beginning of all atoms according to the order of the atoms in SMARTS. The relative order of the remaining atoms remains unchanged.

        Returns:
            The new parser.
        """
        return (f.reset_atom_index(mapping_smarts) for f in self)

    def standard_orient(
        self,
        anchor_list: List[int],
    ) -> Generator[BaseBlockParser, None, None]:
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
            The new parser.
        """
        return (f.standard_orient(anchor_list) for f in self)
