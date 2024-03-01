"""
Author: TMJ
Date: 2024-02-04 11:04:35
LastEditors: TMJ
LastEditTime: 2024-02-12 16:49:38
Description: 请填写简介
"""

import fire

from molop import AutoParser
from molop.io import FileParserBatch


class MolOPCLI:
    """
    CLI for MolOP.
    """

    def __init__(self) -> None:
        self._file_batch: FileParserBatch = None
        self._temp_batch: FileParserBatch = None

    @property
    def temp_batch(self):
        """
        Get the file batch object.
        """
        if self._temp_batch is None:
            return self._file_batch
        else:
            return self._temp_batch

    def auto(self):
        """
        Auto process the current directory.
        """
        self._file_batch = AutoParser(r"*.log")
        self.temp_batch.to_summary_csv()
        self.temp_batch.to_SDF_file()
        return self

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
        self.temp_batch.to_summary_csv(file_dir)
        return self

    def xyz(self, file_dir: str = None):
        """
        Save the XYZ file of all frames of each file.
        """
        self.temp_batch.to_XYZ_file(file_dir)
        return self

    def sdf(self, file_dir: str = None):
        """
        Save the SDF file of all frames of each file.
        """
        self.temp_batch.to_SDF_file(file_dir)
        return self

    def chemdraw(self, file_dir: str = None, frameID=-1, keep3D=True):
        """
        Save the cdxml file of specified frames of each file.
        """
        self.temp_batch.to_chemdraw(file_dir, frameID, keep3D)
        return self

    def gjf(
        self,
        file_dir: str = None,
        charge: int = None,
        multiplicity: int = None,
        prefix: str = "#p opt b3lyp def2svp freq EmpiricalDispersion=GD3BJ NoSymm\n",
        suffix: str = "\n\n",
        template: str = None,
        chk: bool = True,
        oldchk: bool = False,
        frameID: int = -1,
    ):
        """
        Save the GJF file of any frame of each file.
        """
        self.temp_batch.to_GJF_file(
            file_dir,
            charge=charge,
            multiplicity=multiplicity,
            prefix=prefix,
            suffix=suffix,
            template=template,
            chk=chk,
            oldchk=oldchk,
            frameID=frameID,
        )
        return self

    def smiles(self):
        """
        Print the SMILES of last frame of each files.
        """
        for _file in self.temp_batch:
            print(_file[-1].to_standard_SMILES())
        return self

    def charge(self, charge: int):
        """
        Filter the file batch by charge.
        """
        self._temp_batch = self.temp_batch.filter_by_charge(charge)
        return self

    def multi(self, multiplicity: int):
        """
        Filter the file batch by multiplicity.
        """
        self._temp_batch = self.temp_batch.filter_by_multi(multiplicity)
        return self

    def normal(self):
        """
        Filter the file batch by normal judgement.
        """
        self._temp_batch = self.temp_batch.filter_normal()
        return self

    def error(self):
        """
        Filter the file batch by error judgement.
        """
        self._temp_batch = self.temp_batch.filter_error()
        return self

    def ts(self):
        """
        Filter the file batch by TS judgement, which means the structure has a unique imagnary frequency.
        """
        self._temp_batch = self.temp_batch.filter_TS()
        return self

    def format(self, format: str):
        """
        Filter the file batch by format. e.g. "sdf" or ".sdf" are equal.
        """
        self._temp_batch = self.temp_batch.filter_by_format(format)
        return self

    def paths(self):
        """
        Print the file paths of each file.
        """
        for _file in self.temp_batch:
            print(_file.file_path)
        return self

    def end(self):
        """
        End the command chain, stop printing help comments.
        """
        self._file_batch = None


def app():
    fire.Fire(MolOPCLI)
