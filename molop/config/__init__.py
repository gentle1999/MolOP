"""
Author: TMJ
Date: 2024-02-14 14:40:02
LastEditors: TMJ
LastEditTime: 2024-08-04 20:36:07
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

from molop.logger.logger import logger

warnings.simplefilter("ignore", UnitStrippedWarning)


class MolOPConfig:
    def __init__(self) -> None:
        self.show_progress_bar = True
        self.show_log = True
        self._sh = logging.StreamHandler()
        logger.addHandler(self._sh)
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
        self.show_progress_bar = False
        self.show_log = False
        logger.removeHandler(self._sh)

    def verbose(self):
        self.show_progress_bar = True
        self.show_log = True
        logger.addHandler(self._sh)


molopconfig = MolOPConfig()
