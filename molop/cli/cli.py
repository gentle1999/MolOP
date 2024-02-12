"""
Author: TMJ
Date: 2024-02-04 11:04:35
LastEditors: TMJ
LastEditTime: 2024-02-12 16:49:38
Description: 请填写简介
"""
import os

import fire

from molop import AutoParser
from molop.config import molopconfig
from molop.io import FileParserBatch


class MolOPCLI:
    """
    CLI for MolOP.
    """

    def __init__(self) -> None:
        self._file_batch: FileParserBatch = None

    def auto(self):
        """
        Auto process the current directory.
        """
        self._file_batch = AutoParser(r"*.log")
        self._file_batch.to_summary_csv()
        self._file_batch.to_SDF_file()

    def read(
        self,
        file_path: str,
        charge=None,
        multiplicity=None,
        only_extract_structure=False,
        only_last_frame=False,
    ):
        """
        Read the files given and set the file batch object.

        Parameters:
            file_path str:
                use regax to match files.
            charge int:
                forced charge of the molecule, if not given, will use the charge written in the file or 0.
            multiplicity int:
                forced multiplicity of the molecule, if not given, will use the charge written in the file or 1.
            only_extract_structure bool:
                if True, only extract the structure, else extract the whole file.
            only_last_frame bool:
                if True, only extract the last frame, else extract all frames.
        """
        self._file_batch = AutoParser(
            file_path,
            charge,
            multiplicity,
            only_extract_structure=only_extract_structure,
            only_last_frame=only_last_frame,
        )
        return self

    def summary(self, file_dir: str = None):
        """
        Save the summary csv file.
        """
        self._file_batch.to_summary_csv(file_dir)
        return self

    def xyz(self, file_dir: str = None):
        """
        Save the XYZ file of all frames of each file.
        """
        self._file_batch.to_XYZ_file(file_dir)
        return self

    def sdf(self, file_dir: str = None):
        """
        Save the SDF file of all frames of each file.
        """
        self._file_batch.to_SDF_file(file_dir)
        return self

    def chemdraw(self, file_dir: str = None, frameID=-1, keep3D=True):
        """
        Save the cdxml file of specified frames of each file.
        """
        self._file_batch.to_chemdraw(file_dir, frameID, keep3D)
        return self

    def gjf(
        self,
        file_dir: str = None,
        charge: int = None,
        multiplicity: int = None,
        prefix: str = "# g16 gjf \n",
        suffix: str = "\n\n",
        template: str = None,
        frameID: int = -1,
    ):
        """
        Save the GJF file of any frame of each file.
        """
        self._file_batch.to_GJF_file(
            file_dir,
            charge=charge,
            multiplicity=multiplicity,
            prefix=prefix,
            suffix=suffix,
            template=template,
            frameID=frameID,
        )
        return self

    def smiles(self):
        """
        Print the SMILES of last frame of each files.
        """
        for _file in self._file_batch:
            print(_file[-1].to_standard_SMILES())
        return self

    def end(self):
        """
        End the command chain, stop printing help comments.
        """
        self._file_batch = None


def app():
    fire.Fire(MolOPCLI)
