"""
Author: TMJ
Date: 2024-02-14 14:40:02
LastEditors: TMJ
LastEditTime: 2024-08-31 15:24:50
Description: 请填写简介
"""

import logging
import warnings

from openbabel import pybel
from pint.errors import UnitStrippedWarning
from rdkit import RDLogger
from rdkit.Chem.rdFingerprintGenerator import (
    GetMorganGenerator,
    GetAtomPairGenerator,
    GetRDKitFPGenerator,
    GetTopologicalTorsionGenerator,
)

from molop.logger.logger import moloplogger, file_handler, stream_handler

warnings.simplefilter("ignore", UnitStrippedWarning)


class MolOPConfig:
    def __init__(self) -> None:
        self.show_progress_bar = True
        RDLogger.DisableLog("rdApp.*")
        # ob_log_handler = pybel.ob.OBMessageHandler()
        # ob_log_handler.SetOutputLevel(0)
        pybel.ob.obErrorLog.StopLogging()
        self.max_jobs = 16
        self.morgan_fpgen = GetMorganGenerator(
            radius=3, fpSize=1024, includeChirality=True
        )
        self.atompair_fpgen = GetAtomPairGenerator(fpSize=1024, includeChirality=True)
        self.rdkit_fpgen = GetRDKitFPGenerator(fpSize=1024)
        self.topological_torsion_fpgen = GetTopologicalTorsionGenerator(fpSize=1024)
        self.max_structure_recovery_time = 10
        self.allow_spin_change = False

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
        moloplogger.removeHandler(file_handler)

    def log_on(self):
        """
        Turn on log messages to file.
        """
        moloplogger.addHandler(file_handler)

    def set_log_level(self, level: str):
        moloplogger.setLevel(level)


molopconfig = MolOPConfig()
