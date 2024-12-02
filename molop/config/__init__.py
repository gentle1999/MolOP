"""
Author: TMJ
Date: 2024-10-19 09:57:26
LastEditors: TMJ
LastEditTime: 2024-12-02 21:28:20
Description: 请填写简介
"""

import logging
import warnings

from openbabel import pybel
from pint.errors import UnitStrippedWarning
from rdkit import RDLogger
from rdkit.Chem.rdFingerprintGenerator import (
    GetAtomPairGenerator,
    GetMorganGenerator,
    GetRDKitFPGenerator,
    GetTopologicalTorsionGenerator,
)

from molop.logger.logger import moloplogger, stream_handler

warnings.simplefilter("ignore", UnitStrippedWarning)


class MolOPConfig:
    def __init__(self) -> None:
        self.show_progress_bar = True
        RDLogger.DisableLog("rdApp.*")
        pybel.ob.obErrorLog.StopLogging()
        self.max_jobs = 16
        self.morgan_fpgen = GetMorganGenerator(
            radius=3, fpSize=1024, includeChirality=True
        )
        self.atompair_fpgen = GetAtomPairGenerator(fpSize=1024, includeChirality=True)
        self.rdkit_fpgen = GetRDKitFPGenerator(fpSize=1024)
        self.topological_torsion_fpgen = GetTopologicalTorsionGenerator(fpSize=1024)
        self.strict_structure_recovery = False
        self.max_structure_recovery_time = 10
        self.allow_spin_change = False
        self.force_unit_transform = False
        self.parallel_max_size = 8 * 1024**2

    def quiet(self):
        """
        Turn off progress bar and log messages to stdout.
        """
        self.show_progress_bar = False
        moloplogger.removeHandler(stream_handler)

    def verbose(self):
        """
        Turn on progress bar and log messages to stdout.
        """
        self.show_progress_bar = True
        moloplogger.addHandler(stream_handler)

    def log_off(self):
        """
        Turn off log messages to file.
        """
        pass
        # moloplogger.removeHandler(file_handler)

    def log_on(self):
        """
        Turn on log messages to file.
        """
        pass
        # moloplogger.addHandler(file_handler)

    def set_log_level(self, level: str):
        moloplogger.setLevel(level)


molopconfig = MolOPConfig()
