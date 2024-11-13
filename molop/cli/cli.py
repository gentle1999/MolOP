"""
Author: TMJ
Date: 2024-02-04 11:04:35
LastEditors: TMJ
LastEditTime: 2024-02-12 16:49:38
Description: 请填写简介
"""

import fire

from molop import AutoParser
from molop.config import molopconfig
from molop.io import FileParserBatch


class MolOPCLI:
    """
    CLI for MolOP.
    """

    def __init__(self) -> None:
        molopconfig.log_on()
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

    def auto(self, hong=True):
        """
        Auto process the current directory.
        """
        self._file_batch = AutoParser(r"*.log")
        self.temp_batch.to_summary_csv(use_hong_style=hong)

    def read(
        self,
        file_path: str,
        charge=0,
        multiplicity=1,
        n_jobs=-1,
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
            n_jobs=n_jobs,
            only_extract_structure=only_extract_structure,
            only_last_frame=only_last_frame,
        )
        return self

    def summary(
        self,
        file_dir: str = None,
        full: bool = True,
        with_units: bool = True,
        use_hong_style=False,
    ):
        """
        Save the summary csv file.

        Parameters:
            file_dir str:
                the directory to save the summary csv file.
            full bool:
                if True, save all the information, else save only the essential information.
            with_units bool:
                if True, save the information with units, else save the information without units.
        """
        self.temp_batch.to_summary_csv(
            file_dir, full, with_units, use_hong_style=use_hong_style
        )
        return self

    def xyz(self, file_dir: str = None):
        """
        Save the XYZ file of all frames of each file.

        Parameters:
            file_dir str:
                the directory to save the XYZ file.
        """
        self.temp_batch.to_XYZ_file(file_dir)
        return self

    def sdf(self, file_dir: str = None):
        """
        Save the SDF file of all frames of each file.

        Parameters:
            file_dir str:
                the directory to save the SDF file.
        """
        self.temp_batch.to_SDF_file(file_dir)
        return self

    def chemdraw(self, file_dir: str = None, frameID=-1, keep3D=True):
        """
        Save the cdxml file of specified frames of each file.

        Parameters:
            file_dir str:
                the directory to save the cdxml file.
            frameID int:
                the frame ID to save, -1 means the last frame.
            keep3D bool:
                if True, keep the 3D information, else remove the 3D information.
        """
        self.temp_batch.to_chemdraw(file_dir, frameID, keep3D)
        return self

    def gjf(
        self,
        file_dir: str = None,
        charge: int = None,
        multiplicity: int = None,
        prefix: str = "",
        suffix: str = "",
        template: str = None,
        chk: bool = True,
        oldchk: bool = False,
        frameID: int = -1,
    ):
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

    def smiles(self, frameID: int = -1):
        """
        Print the SMILES of the specified frame of each files.

        Parameters:
            frameID int:
                the frame ID to print, -1 means the last frame.
        """
        for _file in self.temp_batch:
            try:
                print(_file[frameID].to_canonical_SMILES())
            except AttributeError:
                print(f"No SMILES found in {_file.file_path} frame {frameID}")
        return self

    def charge(self, charge: int):
        """
        Filter the file batch by charge.

        Parameters:
            charge int:
                the charge to filter.
        """
        self._temp_batch = self.temp_batch.filter_by_charge(charge)
        return self

    def multi(self, multiplicity: int):
        """
        Filter the file batch by multiplicity.

        Parameters:
            multiplicity int:
                the multiplicity to filter.
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
        molopconfig.log_off()
        return self

    def log_on(self):
        """
        Turn on the log.
        """
        molopconfig.log_on()
        return self


def app():
    fire.Fire(MolOPCLI)
