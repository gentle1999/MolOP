"""
Author: TMJ
Date: 2024-02-14 14:40:02
LastEditors: TMJ
LastEditTime: 2024-02-15 18:13:15
Description: 请填写简介
"""

import os
from typing import List, Tuple, Union

import pandas as pd
from joblib import Parallel, cpu_count, delayed

from molop.config import molopconfig
from molop.io.bases.molblock_base import BaseBlockParser
from molop.io.coords_file.GJFBlockParser import GJFBlockParser
from molop.io.coords_file.SDFBlockParser import SDFBlockParser
from molop.io.coords_file.XYZBlockParser import XYZBlockParser
from molop.io.qm_file.G16FCHKBlockParser import G16FCHKBlockParser
from molop.io.qm_file.G16LOGBlockParser import G16LOGBlockParser
from molop.io.qm_file.XTBOUTBlockParser import XTBOUTBlockParser
from molop.unit import atom_ureg

BLOCKTYPES = Union[
    BaseBlockParser,
    G16LOGBlockParser,
    G16FCHKBlockParser,
    XTBOUTBlockParser,
    GJFBlockParser,
    SDFBlockParser,
    XYZBlockParser,
]
QMBLOCKTYPES = Union[
    G16LOGBlockParser,
    G16FCHKBlockParser,
    XTBOUTBlockParser,
]


class BaseFileParser:
    """
    Base class for multi-frame parsers.

    Attributes:
        _file_path str:
            File path of the file to be parsed.
        __frames List[BlockType]:
            List of parsed frames.
        __index int:
            Current frame index.
        _allowed_formats Tuple[str]:
            Allowed file formats.
    """

    _file_path: str
    __frames: List[BLOCKTYPES]
    __index: int
    _allowed_formats: Tuple[str]

    def __init__(self, file_path: str) -> None:
        self._file_path = os.path.abspath(file_path)
        self.__frames = []
        self.__index: int = 0

    def _check_formats(self, file_path: str) -> None:
        """
        Check if the file exists and is a file, and if its format is allowed.

        Parameters:
            file_path (str): The path of the file to be checked.

        Raises:
            FileNotFoundError: If the file does not exist.
            IsADirectoryError: If the file is a directory.
            ValueError: If the file format is not allowed.
        """
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"{file_path} not found.")
        elif not os.path.isfile(file_path):
            raise IsADirectoryError(f"{file_path} is not a file.")
        else:
            _, file_format = os.path.splitext(file_path)
            if file_format not in self._allowed_formats:
                raise ValueError(
                    f"File format {file_format} not in {self._allowed_formats}"
                )
        self._file_format = file_format

    def _post_parse(self) -> None:
        """
        Post-parsing operations.
        """
        pass

    def __iter__(self):
        self.__index = 0
        return self

    def __next__(
        self,
    ) -> BLOCKTYPES:
        if self.__index >= len(self):
            raise StopIteration
        else:
            self.__index += 1
            return self.__frames[self.__index - 1]

    def __getitem__(self, frameID: int) -> BLOCKTYPES:
        return self.__frames[frameID]

    def __len__(self) -> int:
        return len(self.__frames)

    def __str__(self) -> str:
        return f"{self.__class__.__name__}({os.path.basename(self._file_path)})"

    def __hash__(self) -> int:
        return hash(str(self))

    def __repr__(self) -> str:
        return f"{str(self)}(frame num: {len(self)})"

    @property
    def file_path(self) -> str:
        """
        Get the file path of the object.

        Returns:
            str: The file path of the object.
        """
        return self._file_path

    @property
    def file_name(self) -> str:
        """
        Get the file name from the file path.

        Returns:
            str: The file name.
        """
        return os.path.basename(self._file_path)

    @property
    def pure_filename(self) -> str:
        return os.path.splitext(self.file_name)[0]

    @property
    def file_format(self) -> str:
        """
        Get the file format of the object.
        Returns:
            str: The file format of the object.
        """
        return self._file_format

    def _check_path(self, file_path: str = None, format: str = ".xyz"):
        if file_path is None:
            return os.path.splitext(self._file_path)[0] + format
        if os.path.isdir(file_path):
            return os.path.join(
                file_path,
                self.pure_filename + format,
            )
        return file_path

    @property
    def frames(
        self,
    ) -> List[BLOCKTYPES]:
        """
        Get the list of parsed frames.

        Returns:
            List[BlockType]: The list of parsed frames.
        """
        return self.__frames

    def append(
        self,
        frame: BLOCKTYPES,
    ):
        """
        Append a frame to the list of frames.

        Parameters:
            frame (BlockType): The frame to be appended.

        Raises:
            TypeError: If the type of the frame is not a subclass of BaseBlockParser.
        """
        if not issubclass(type(frame), BaseBlockParser):
            raise TypeError(f"{type(frame)} is not a subclass of {BaseBlockParser}")
        frame._frameID = len(self.__frames)
        if frame._frameID > 0:
            self.__frames[frame._frameID - 1]._next_block = frame
            frame._prev_block = self.__frames[frame._frameID - 1]
        self.__frames.append(frame)

    def to_XYZ_block(self) -> str:
        """
        Convert all frames to a string representation in XYZ format.

        Returns:
            str: The string representation of the object in XYZ format.
        """
        return "\n".join([frame.to_XYZ_block() for frame in self.__frames])

    def to_SDF_block(self) -> str:
        """
        Convert all frames to a string representation in SDF format.

        Returns:
           str: The string representation in SDF format.
        """
        return "$$$$\n".join([frame.to_SDF_block() for frame in self.__frames])

    def to_GJF_block(
        self,
        file_path: str = None,
        charge: int = None,
        multiplicity: int = None,
        prefix: str = f"#p opt b3lyp def2svp freq EmpiricalDispersion=GD3BJ NoSymm\n",
        suffix="",
        template: str = None,
        chk: bool = True,
        oldchk: bool = False,
        frameID=-1,
    ) -> str:
        """
        Convert the frame at the specified index to a Gaussian input block in the Gaussian 16 format (GJF).

        Parameters:
            charge int:
                The forced charge. If specified, will be used to overwrite the charge in the gjf file.
            multiplicity int:
                The forced multiplicity. If specified, will be used to overwrite the multiplicity in the gjf file.
            template str:
                path to read a gjf file as a template.
            prefix str:
                prefix to add to the beginning of the gjf file, priority is lower than template.
            suffix str:
                suffix to add to the end of the gjf file, priority is lower than template.
            chk bool:
                If true, add the chk keyword to the link0 section. Will use the file name as the chk file name.
            oldchk bool:
                If true, add the oldchk keyword to the link0 section. Will use the file name as the chk file name.

        Returns:
            A modified GJF block.
        """
        return self.__frames[frameID].to_GJF_block(
            file_path=file_path,
            charge=charge,
            multiplicity=multiplicity,
            prefix=prefix,
            suffix=suffix,
            template=template,
            chk=chk,
            oldchk=oldchk,
        )

    def to_XYZ_file(self, file_path: str = None):
        """
        Write the XYZ file.

        Parameters:
            file_path str:
                The file path. If not specified, will be generated in situ.
        Returns:
            The absolute path of the XYZ file.
        """
        _file_path = self._check_path(file_path, ".xyz")
        with open(_file_path, "w") as f:
            f.write(self.to_XYZ_block())
        f.close()
        return os.path.abspath(_file_path)

    def to_SDF_file(self, file_path: str = None):
        """
        Write the SDF file.

        Parameters:
            file_path str:
                The file path. If not specified, will be generated in situ.
        Returns:
            The absolute path of the SDF file.
        """
        _file_path = self._check_path(file_path, ".sdf")
        with open(_file_path, "w") as f:
            f.write(self.to_SDF_block())
        f.close()
        return os.path.abspath(_file_path)

    def to_GJF_file(
        self,
        file_path: str = None,
        charge: int = None,
        multiplicity: int = None,
        template: str = None,
        prefix: str = f"#p opt b3lyp def2svp freq EmpiricalDispersion=GD3BJ NoSymm\n",
        suffix="\n\n",
        chk: bool = True,
        oldchk: bool = False,
        frameID=-1,
    ):
        """
        Write the GJF file.

        Parameters:
            file_path str:
                The path to write the GJF file. If not specified, will be generated in situ.
            charge int:
                The forced charge. If specified, will be used to overwrite the charge in the gjf file.
            multiplicity int:
                The forced multiplicity. If specified, will be used to overwrite the multiplicity in the gjf file.
            template str:
                path to read a gjf file as a template.
            prefix str:
                prefix to add to the beginning of the gjf file, priority is lower than template.
            suffix str:
                suffix to add to the end of the gjf file, priority is lower than template.
            frameID int:
                The frame ID to write.
            chk bool:
                If true, add the chk keyword to the link0 section. Will use the file name as the chk file name.
            oldchk bool:
                If true, add the oldchk keyword to the link0 section. Will use the file name as the chk file name.
        Returns:
            The path to the GJF file.
        """
        _file_path = self._check_path(file_path, ".gjf")
        with open(_file_path, "w") as f:
            f.write(
                self.to_GJF_block(
                    file_path=file_path,
                    charge=charge,
                    multiplicity=multiplicity,
                    prefix=prefix,
                    suffix=suffix,
                    template=template,
                    chk=chk,
                    oldchk=oldchk,
                    frameID=frameID,
                )
            )
        f.close()
        return os.path.abspath(_file_path)

    def to_chemdraw(self, file_path: str = None, frameID=-1, keep3D=True):
        """
        Write the ChemDraw file.

        Parameters:
            file_path str:
                The path to write the ChemDraw file. If not specified, will be generated in situ.
            frameID int:
                The frame ID to write.
            keep3D bool:
                Whether to keep the 3D information.
        Returns:
            The path to the ChemDraw file.
        """
        return self.__frames[frameID].to_chemdraw(file_path, keep3D=keep3D)

    def geometry_analysis_df(
        self, key_atoms: List[List[int]], precision: int = 1, one_start=False
    ) -> pd.DataFrame:
        """
        Get the geometry infos among the atoms with all frames in the file

        Parameters:
            key_atoms List[List[int]]:
                A list of list of index of the atoms, starts from 0:

                    - If the length of atom_idxs is 2, the bond length with unit Angstrom between the two atoms will be returned.
                    - If the length of atom_idxs is 3, the angle with unit degree between  the three atoms will be returned.
                    - If the length of atom_idxs is 4, the dihedral angle with unit degree between the four atoms will be returned.
            precision int:
                The precision of the geometry analysis. Default is 1. e.g. 1 means 1.0001 ==> 1.0
            one_start bool:
                If true, consider atom index starts from 1, so let index value subtracts 1 for all the atoms

        Returns:
            A pandas DataFrame:
                The columns of the DataFrame is the atom_idxs, and the index are the frames.
        """
        values = {
            "-".join([str(idx) for idx in atom_idxs]): [
                round(molblock.geometry_analysis(atom_idxs, one_start), precision)
                for molblock in self.__frames
            ]
            for atom_idxs in key_atoms
        }
        return pd.DataFrame(values)

    def geometry_analysis(
        self,
        key_atoms: List[List[int]],
        file_path: str = None,
        precision: int = 1,
        one_start=False,
    ):
        """
        Get the geometry infos among the atoms with all frames in the file

        Parameters:
            key_atoms List[List[int]]:
                A list of list of index of the atoms, starts from 0

                    - If the length of atom_idxs is 2, the bond length with unit Angstrom between the two atoms will be returned.
                    - If the length of atom_idxs is 3, the angle with unit degree between  the three atoms will be returned.
                    - If the length of atom_idxs is 4, the dihedral angle with unit degree between the four atoms will be returned.
            file_path str:
                The path of the csv file to be saved. If None, the file will be saved in the same directory of the file_path.
            precision int:
                The precision of the geometry analysis. Default is 1. e.g. 1 means 1.0001 ==> 1.0
            one_start bool:
                If true, consider atom index starts from 1, so let index value subtracts 1 for all the atoms

        Returns:
            file_path
        """
        df = self.geometry_analysis_df(
            key_atoms, precision=precision, one_start=one_start
        )
        if file_path is None:
            file_path = self._file_path
        if os.path.isdir(file_path):
            raise IsADirectoryError(f"{file_path} is a directory.")
        file_path = os.path.splitext(file_path)[0] + "_geometry_analysis.csv"
        df.to_csv(file_path)
        return file_path

    def replace_substituent(
        self,
        query_smi: str,
        replacement_smi: str,
        bind_idx: int = None,
        replace_all=False,
        attempt_num: int = 10,
        frameID=-1,
    ) -> BaseBlockParser:
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
                The frame ID to replace.
        Returns:
            The new parser.
        """
        return self.__frames[frameID].replace_substituent(
            query_smi, replacement_smi, bind_idx, replace_all, attempt_num
        )

    def reset_atom_index(
        self,
        mapping_smarts: str,
        frameID: int = -1,
    ) -> BaseBlockParser:
        """
        Reset the atom index of the molecule according to the mapping SMARTS.

        Parameters:
            mapping_smarts str:
                The SMARTS to query the molecule substructure.
                The queried atoms will be renumbered and placed at the beginning of all atoms according to the order of the atoms in SMARTS. The relative order of the remaining atoms remains unchanged.
            frameID int:
                The frame ID to reset.
        Returns:
            The new parser.
        """
        return self.__frames[frameID].reset_atom_index(mapping_smarts)

    def standard_orient(
        self,
        anchor_list: List[int],
        frameID: int = -1,
    ) -> BaseBlockParser:
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
            frameID int:
                The frame ID to standardize orientation.

        Returns:
            The new parser.
        """
        return self.__frames[frameID].standard_orient(anchor_list)

    def to_summary_df(self) -> pd.DataFrame:
        """
        Returns:
            A pandas DataFrame containing the summary information of the parser.
        """
        self.recover_structures()
        return pd.concat(
            [frame.to_summary_series() for frame in self.__frames], axis=1
        ).T

    def recover_structures(self):
        return [frame.to_standard_SMILES() for frame in self.frames]


class BaseQMFileParser(BaseFileParser):
    """
    Base class for QM multi-frame parsers.

    Attributes:
        _parameter_comment str:
            The parameter comment in the file.
        _only_extract_structure bool:
            Whether to only extract the structure.
        _only_last_frame bool:
            Whether to only extract the last frame.
        _version str:
            The version of the QM software.
    """

    _only_extract_structure: bool
    _only_last_frame: bool

    def __init__(
        self,
        file_path: str,
        only_extract_structure=False,
        only_last_frame=False,
    ) -> None:
        super().__init__(file_path)
        self.parameter_comment: str = None
        self.basis: str = None
        self.functional: str = None
        self.temperature: float = 298.15

        self.solvent_model: str = None
        self.solvent: str = None

        self._only_extract_structure: bool = only_extract_structure
        self._only_last_frame = only_last_frame
        self.version: str = None

    def __repr__(self) -> str:
        return (
            f"{str(self)}\n"
            + f"functional: {self.functional}\n"
            + f"basis: {self.basis}\n"
            + f"solvent_model: {self.solvent_model}\n"
            + f"solvent: {self.solvent}\n"
            + f"frame num: {len(self)}"
        )
